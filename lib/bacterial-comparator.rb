# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maximilien1er@gmail.com
# review:
# date:    	15-02-24
# version: 	0.0.1
# licence:

require 'bio'
require 'fileutils'
require 'parallel'
require 'helper'

require 'bacterial-annotator/sequence-annotation'

class BacterialComparator

  attr_reader :genomes_list, :stats

  # Initialize BacterialAnnotator
  # options[:input], options[:refgenome], ROOT, options[:outdir], options)
  def initialize options, root

    @root = root
    @outdir = File.expand_path(File.dirname(options[:outdir])) + "/#{File.basename(options[:outdir])}"
    Dir.mkdir(@outdir) if ! Dir.exists? @outdir
    @genomes_list = options[:genomes_list]
    @proc = options[:proc].to_i
    @phylo_nb_genes = options[:phylo_nb_genes]
    @refgenome_file = options[:refgenome]
    @refgenome = ""

    if ["fasttree","raxml"].include? options[:software]
      @software = options[:software]
    else
      @software = "fasttree"
    end

    min_cov = options[:min_cov].to_f
    min_pid = options[:pidentity].to_f
    if min_cov > 1
      min_cov = min_cov/100
    end
    if min_pid > 1
      min_pid = min_pid/100
    end

    @min_cov = min_cov
    @min_pid = min_pid

    @aln_opt = options[:align].downcase
    @run_phylo = 0
    if options[:phylogeny] == 1
      @bootstrap = options[:bootstrap]
      @run_phylo = 1
    end

    @ref_prot = get_ref_prot
    @synteny = read_prot_synteny
    @stats = extract_syntenic_fasta min_cov, min_pid

  end


  def run_comparison

    run_mafft_aln

    run_phylo if @run_phylo != 0

  end

  def read_prot_synteny

    puts "# Reading genome synteny files START.."
    # print(@genomes_list)

    start_time = Time.now
    synteny = {}

    # If genome is genbank (.gbk or .gb) then do syntheny first
    new_genomes_dir = @outdir+"/new-genomes/"
    Dir.mkdir(new_genomes_dir) if ! Dir.exists? new_genomes_dir

    @genomes_list.each_with_index do |g,i|

      genome_synteny = []

      if File.directory?(g)
        @genomes_list[i] = g
      else

        if ! File.exists? g
          puts "Downloading genbank file because it doesn't exists"
          id = g
          if id =~ /.gbk/i
            id = id.gsub(".gbk","")
          else
            g = id+".gbk"
          end
          puts "#{id} #{g}"
          Helper.download_genbank id, g
        end

        genome_dir = new_genomes_dir + "/" + File.basename(g).gsub(".gbk","")
        Dir.mkdir(genome_dir) if ! Dir.exists? genome_dir
        genome_to_annotate = SequenceAnnotation.new(@root,
                                                    genome_dir,
                                                    g,
                                                    "refGbk")
        query_prot_file = Dir["#{genome_dir}/*.pep"][0]
        query_gene_file = Dir["#{genome_dir}/*.dna"][0]

        File.symlink(query_prot_file, File.dirname(query_prot_file)+"/Proteins.fa")
        File.symlink(query_gene_file, File.dirname(query_gene_file)+"/Genes.fa")

        run_synteny_prot @root, genome_dir, @ref_prot_file, query_prot_file

        g = genome_dir
        @genomes_list[i] = genome_dir

      end


      file = File.open("#{g}/Prot-Synteny.tsv", "r")
      l = file.gets             # skip header
      while l = file.gets
        # AAK98805.1      spr0001 453     1.0     100.0   ABAC01000005_14 453     1.0	1|0
        lA = l.chomp.split("\t")
        synteny[lA[0]] =  [] if ! synteny.has_key? lA[0]
        synteny[lA[0]] << {ref_cov: lA[3].to_f, pId: lA[4].to_f, query_prot: lA[5], query_cov: lA[7].to_f}
        genome_synteny << lA[0]
      end

      @ref_prot.each do |ref_prot|
        if ! genome_synteny.include? ref_prot
          synteny[lA[0]] << {ref_cov: "-", pId: "-", query_prot: "-", query_cov: "-"}
        end
      end
      file.close
    end
    end_time = Time.now
    c_time = Helper.sec2str(end_time-start_time)

    puts "# Reading genome synteny files [DONE] (in #{c_time})"

    synteny

  end

  def get_ref_prot

    ref_prot = []
    if File.exist?("#{@genomes_list[0]}/*.pep")
      pep_file = Dir["#{@genomes_list[0]}/*.pep"]
      @ref_prot_file = pep_file[0]
      flatfile = Bio::FlatFile.auto("#{@ref_prot_file}")
      flatfile.each_entry do |entry|
        ref_prot << entry.definition.split(" ")[0]
      end
      flatfile.close
      dna_file = Dir["#{@genomes_list[0]}/*.dna"]
      @ref_dna_file = dna_file[0]
    else
      if @refgenome_file == ""
        abort "You need to provide a reference genome to add a non-annotated genome"
      elsif @refgenome == ""
        new_genomes_dir = @outdir+"/new-genomes/"
        Dir.mkdir(new_genomes_dir) if ! Dir.exists? new_genomes_dir
        refgenome_dir = new_genomes_dir + "/" + File.basename(@refgenome_file).gsub(".gbk","")
        Dir.mkdir(refgenome_dir) if ! Dir.exists? refgenome_dir
        @refgenome = SequenceAnnotation.new(@root,
                                            refgenome_dir,
                                            @refgenome_file,
                                            "refGbk")
        pep_file = Dir["#{refgenome_dir}/*.pep"]
        @ref_prot_file = pep_file[0]
        flatfile = Bio::FlatFile.auto("#{@ref_prot_file}")
        flatfile.each_entry do |entry|
          ref_prot << entry.definition.split(" ")[0]
        end
        flatfile.close
        dna_file = Dir["#{refgenome_dir}/*.dna"]
        @ref_dna_file = dna_file[0]
      end
    end

    ref_prot

  end

  # load all id => sequences from multifasta
  def load_genome_cds file

    proteins = {}
    flatfile = Bio::FlatFile.auto(file)
    flatfile.each_entry do |entry|
      name = entry.definition.split(" ")[0]
      bioseq = Bio::Sequence.auto(entry.seq)
      out = bioseq.output_fasta("#{name}",60)
      proteins[name] = out
    end

    proteins

  end


  def build_multifasta synteny_list

    pep_out_dir = "#{@outdir}/align-genes-pep"

    ref_proteins = load_genome_cds(@ref_prot_file)
    synteny_list.each do |k,v|
      pep_out = File.open(pep_out_dir+"/#{k}.pep", "w")
      pep_out.write(ref_proteins[k])
      pep_out.close
    end

    @genomes_list.each_with_index do |g,i|

      genome_proteins = load_genome_cds("#{g}/Proteins.fa")
      synteny_list.each do |k,v|
        pep_out = File.open(pep_out_dir+"/#{k}.pep", "a")
        # puts "#{v[i][:query_prot]}"
        if ! genome_proteins.has_key? v[i][:query_prot]
          puts "#{g} doesn't have the query prot"
        end
        pep_out.write(genome_proteins[v[i][:query_prot]])
        pep_out.close
      end

    end

    dna_out_dir = "#{@outdir}/align-genes-dna"
    ref_genes = load_genome_cds(@ref_dna_file)
    synteny_list.each do |k,v|
      dna_out = File.open(dna_out_dir+"/#{k}.dna", "w")
      dna_out.write(ref_genes[k])
      dna_out.close
    end

    @genomes_list.each_with_index do |g,i|

      genome_genes = load_genome_cds("#{g}/Genes.fa")
      synteny_list.each do |k,v|
        dna_out = File.open(dna_out_dir+"/#{k}.dna", "a")
        dna_out.write(genome_genes[v[i][:query_prot]])
        dna_out.close
      end

    end

  end

  # extract and dump multifasta for syntenic genes and proteins
  def extract_syntenic_fasta min_cov, min_pid

    puts "# Extracting Proteins and Genes multifasta  START.."
    start_time = Time.now

    nb_of_syntenic = 0
    stats = {}
    stats[:syntenic] = []
    fout = File.open("#{@outdir}/cds-synteny.tsv", "w")
    genomes_name = []
    @genomes_list.each do |g|
      genomes_name.push(File.basename(g))
    end
    fout.write("Gene\t"+genomes_name.join("\t")+"\n")

    to_build_multifasta = []

    @synteny.each do |k,v|

      is_syntenic = 1
      v.each do |v_|
        if v_[:query_cov] == "-"
          is_syntenic = 0
          break
        elsif v_[:query_cov] > min_cov and
             v_[:ref_cov] > min_cov and
             v_[:pId] > min_pid
        # synteny -> great !
        else
          is_syntenic = 0
          break
        end
      end

      fout.write("#{k}")
      if is_syntenic == 1
        nb_of_syntenic += 1
        # build_multifasta k, v
        to_build_multifasta << [k,v]
        v.each do |x|
          fout.write("\t#{x[:query_prot]}|#{x[:pId]}|#{x[:query_cov]}|#{x[:ref_cov]}")
          stats[:syntenic] << k
        end
        fout.write("\n")
      else
        v.each do |x|
          if x[:pId] == 0.0
            fout.write("\t-")
          else
            fout.write("\t#{x[:query_prot]}|#{x[:pId]}|#{x[:query_cov]}|#{x[:ref_cov]}")
          end
        end
        fout.write("\n")
      end

    end

    fout.close

    pep_out_dir = "#{@outdir}/align-genes-pep"
    dna_out_dir = "#{@outdir}/align-genes-dna"
    Dir.mkdir(pep_out_dir) if ! Dir.exists? pep_out_dir
    Dir.mkdir(dna_out_dir) if ! Dir.exists? dna_out_dir

    synteny_list = to_build_multifasta.each_slice((to_build_multifasta.length/@proc)+1).to_a

    Parallel.map(synteny_list, in_processes: @proc) { |list|
      build_multifasta list
    }

    end_time = Time.now
    c_time = Helper.sec2str(end_time-start_time)
    puts "# Extracting Proteins and Genes multifasta  [DONE] (in #{c_time})"

    stats[:nb_of_syntenic] = nb_of_syntenic
    #puts "   Syntenic genes : " + nb_of_syntenic.to_s + " / " + @ref_prot.length.to_s

  end

  def mafft_align f

    trying = 0
    begin
      cmd = system("#{@root}/mafft.linux --quiet #{f} > #{f}.aln")
      if File.size("#{f}.aln") == 0
        puts "File size of 0.. --#{f}--"
        puts "Command used : #{@root}/mafft.linux --quiet #{f} > #{f}.aln"
        fail
      else
        status = "OK"
        status = "FAILED" if cmd != true
        # puts "Alignment #{f} : #{status}"
      end
    rescue
      if trying < 3
        trying += 1
        retry
      end
      status = "FAILED"
      puts "Alignment #{f} : #{status}"
    end

  end


  def mafft_align_all_pep

    puts "# Sequence alignments - individual proteins a.a. (MAFFT)  START.."
    start_time = Time.now

    ori_dir = Dir.pwd
    Dir.chdir("#{@outdir}/align-genes-pep/")

    is_done = 1
    if Dir["*.pep"].length == Dir["*.aln"].length
      Dir["*.aln"].each do |a|
        if File.size(a) == 0
          is_done = 0
        end
      end
    else
      is_done = 0
    end

    if is_done==0
      Parallel.map(Dir["*.pep"], in_processes: @proc) { |f|
        mafft_align f
      }
    end

    end_time = Time.now
    c_time = Helper.sec2str(end_time-start_time)
    puts "# Sequence alignments - individual proteins a.a. (MAFFT)  [DONE] (in #{c_time})"

    # FIXME ugly hack to find out the reference genome
    Dir.chdir(ori_dir)
    ref_id = @ref_prot_file.split('/')[-1].gsub(".pep","")

    concat_alignments "#{@outdir}/align-genes-pep.all.fasta", ref_id

    Dir.chdir(ori_dir)

  end

  def mafft_align_all_dna

    puts "# Sequence alignments - individual genes dna (MAFFT)  START.."
    start_time = Time.now

    ori_dir = Dir.pwd
    Dir.chdir("#{@outdir}/align-genes-dna/")

    is_done = 1
    if Dir["*.dna"].length == Dir["*.aln"].length
      Dir["*.aln"].each do |a|
        if File.size(a) == 0
          is_done = 0
        end
      end
    else
      is_done = 0
    end

    if is_done == 0
      Parallel.map(Dir["*.dna"], in_processes: @proc) { |f|
        mafft_align f
      }
    end

    # ugly hack to find out the reference genome FIXME
    Dir.chdir(ori_dir)
    ref_id = @ref_prot_file.split('/')[-1].gsub(".pep","")

    end_time = Time.now
    c_time = Helper.sec2str(end_time-start_time)
    puts "# Sequence alignments - individual genes dna (MAFFT)  [DONE] (in #{c_time})"

    concat_alignments "#{@outdir}/align-genes-dna.all.fasta", ref_id

    Dir.chdir(ori_dir)

  end


  def concat_alignments outfile, ref_id

    if File.exists?("#{outfile}") and File.size("#{outfile}") > 0
      puts "..Alignment concatenated file already exists, skipping."
      return
    end

    fout = File.open("#{outfile}", "w")

    seq = ""
    aln_dir = outfile.split(".")[0..-3].join(".")

    Dir["#{aln_dir}/*.aln"].each do |f|
      flat = Bio::FlatFile.auto(f)
      ref_seq = flat.entries[0]
      flat.close
      seq += ref_seq.seq
    end

    bioseq = Bio::Sequence.auto(seq)
    out = bioseq.output_fasta("#{ref_id}", 60)
    fout.write(out)

    for i in 1..@genomes_list.length
      seq = ""
      Dir["#{aln_dir}/*.aln"].each do |f|
        flat = Bio::FlatFile.auto(f)
        j=0
        flat.each_entry do |entry|
          if j<i
            j+=1
            next
          elsif i == j
            seq += entry.seq
            j+=1
          else
            break
          end
        end
        flat.close
      end
      bioseq = Bio::Sequence.auto(seq)
      # get the file name without path prefix and extension
      genome_name = genomes_list[i-1].split("/")[-1].split(".")[0]
      out = bioseq.output_fasta(genome_name,60)
      fout.write(out)
    end

    fout.close

  end

  def run_mafft_aln

    if @aln_opt == "both"
      mafft_align_all_pep
      mafft_align_all_dna
    elsif @aln_opt == "prot"
      mafft_align_all_pep
    elsif @aln_opt == "dna"
      mafft_align_all_dna
    end

  end

  def raxml_tree_dna bt
    puts "# Genes DNA tree creation (RAXML)  START.."
    start_time = Time.now
    ori_dir = Dir.pwd
    Dir.chdir(@outdir)
    Dir.mkdir("tree-genes-dna") if ! Dir.exists?("tree-genes-dna")
    current_dir = Dir.pwd
    tree_dir = "#{current_dir}/tree-genes-dna"
    cmd = system("#{@root}/raxml.linux -T #{@proc} -f d -N #{bt} -s align-genes-dna.all.fasta  -m GTRGAMMA -p 123454321 -n DnaTree -w #{tree_dir} > /dev/null 2>&1")
    cmd = system("cat #{tree_dir}/RAxML_result.DnaTree.RUN.* >> #{tree_dir}/RAxML_result.BS")
    cmd = system("#{@root}/raxml.linux -T #{@proc} -f b -z #{tree_dir}/RAxML_result.BS -t #{tree_dir}/RAxML_bestTree.DnaTree -m GTRGAMMA -n DNA_BS_TREE -w #{tree_dir} > /dev/null 2>&1")
    cmd = system("ln -s #{tree_dir}/RAxML_bipartitionsBranchLabels.DNA_BS_TREE #{tree_dir}/../tree-genes-dna.nwk")
    Dir.chdir(ori_dir)
    end_time = Time.now
    c_time = Helper.sec2str(end_time-start_time)
    puts "# Genes DNA tree creation (RAXML)  [DONE] (in #{c_time})"
  end

  def raxml_tree_pep bt
    puts "# Proteins AA tree creation (RAXML)  START.."
    start_time = Time.now
    ori_dir = Dir.pwd
    Dir.chdir(@outdir)
    Dir.mkdir("tree-genes-pep") if ! Dir.exists?("tree-genes-pep")
    current_dir = Dir.pwd
    tree_dir = "#{current_dir}/tree-genes-pep"
    cmd = system("#{@root}/raxml.linux -T #{@proc} -f d -N #{bt} -s align-genes-pep.all.fasta  -m PROTGAMMAAUTO -p 123454321 -n PepTree -w #{tree_dir} > /dev/null 2>&1")
    cmd = system("cat #{tree_dir}/RAxML_result.PepTree.RUN.* >> #{tree_dir}/RAxML_result.BS")
    cmd = system("#{@root}/raxml.linux -T #{@proc} -f b -z #{tree_dir}/RAxML_result.BS -t #{tree_dir}/RAxML_bestTree.PepTree -m PROTGAMMAAUTO -n PEP_BS_TREE -w #{tree_dir} > /dev/null 2>&1")
    cmd = system("ln -s #{tree_dir}/RAxML_bipartitionsBranchLabels.PEP_BS_TREE #{tree_dir}/../tree-proteins-aa.nwk")
    Dir.chdir(ori_dir)
    end_time = Time.now
    c_time = Helper.sec2str(end_time-start_time)
    puts "# Proteins AA tree creation (RAXML)  [DONE] (in #{c_time})"
  end


  def fasttree_tree_dna bt
    puts "# Genes DNA tree creation (FastTree)  START.."
    start_time = Time.now
    ori_dir = Dir.pwd
    Dir.chdir(@outdir)
    Dir.mkdir("tree-genes-dna") if ! Dir.exists?("tree-genes-dna")
    current_dir = Dir.pwd
    cmd = system("export OMP_NUM_THREADS=#{@proc} && #{@root}/fasttree.linux -nosupport -fastest -nt -gtr align-genes-dna.all.fasta > tree-genes-dna.nwk")
    Dir.chdir(ori_dir)
    end_time = Time.now
    c_time = Helper.sec2str(end_time-start_time)
    puts "# Genes DNA tree creation (FastTree)  [DONE] (in #{c_time})"
  end


  def fasttree_tree_pep bt
    puts "# Proteins AA tree creation (FastTree)  START.."
    start_time = Time.now
    ori_dir = Dir.pwd
    Dir.chdir(@outdir)
    Dir.mkdir("tree-genes-pep") if ! Dir.exists?("tree-genes-pep")
    current_dir = Dir.pwd
    cmd = system("export OMP_NUM_THREADS=#{@proc} && #{@root}/fasttree.linux -nosupport -fastest align-genes-pep.all.fasta > tree-proteins-aa.nwk")
    Dir.chdir(ori_dir)
    end_time = Time.now
    c_time = Helper.sec2str(end_time-start_time)
    puts "# Proteins AA tree creation (FastTree)  [DONE] (in #{c_time})"
  end


  def run_phylo

    if @software == "raxml"
      if @aln_opt == "both"
        raxml_tree_dna @bootstrap
        raxml_tree_pep @bootstrap
      elsif @aln_opt == "prot"
        raxml_tree_pep @bootstrap
      elsif @aln_opt == "dna"
        raxml_tree_dna @bootstrap
      end
    elsif @software == "fasttree"
      if @aln_opt == "both"
        fasttree_tree_dna @bootstrap
        fasttree_tree_pep @bootstrap
      elsif @aln_opt == "prot"
        fasttree_tree_pep @bootstrap
      elsif @aln_opt == "dna"
        fasttree_tree_dna @bootstrap
      end
    end

  end


  def get_fasta_length fasta
    flatfile = Bio::FlatFile.auto(fasta)
    prot_lengths = {}
    flatfile.each_entry do |entry|
      prot_id = entry.definition.split(" ")[0]
      prot_length = entry.length
      prot_lengths[prot_id] = prot_length
    end
    flatfile.close
    prot_lengths
  end


  def run_synteny_prot root, outdir, ref_prot_file, query_prot_file

    puts query_prot_file
    puts ref_prot_file

    ref_synteny_prot = SequenceSynteny.new(root,
                                           outdir,
                                           query_prot_file,
                                           ref_prot_file,
                                           "Prot-Ref",
                                           @min_cov,
                                           @min_cov,
                                           "prot")

    print "# Running alignment with Reference Genome CDS (diamond).."
    start_time = Time.now
    ref_synteny_prot.run_diamond
    end_time = Time.now
    c_time = Helper.sec2str(end_time - start_time)
    print "done (#{c_time})\n"

    ref_synteny_prot.extract_hits :refgenome

    synteny_file = File.open("#{outdir}/Prot-Synteny.tsv","w")
    synteny_file.write("RefLocusTag\tRefProtID\tRefLength\tRefCoverage\tIdentity\tQueryGene\tQueryLength\tQueryCoverage\tQueryPartial\n")
    ref_annotated = {}

    ref_synteny_prot.query_sequences.each do |prot, syn_val|

      next if ! syn_val.has_key? :homology
      next if syn_val[:homology][:assert_cutoff].inject(:+) < 3
      next if ref_annotated.has_key? syn_val[:homology][:hits][0] and
        ref_annotated[syn_val[:homology][:hits][0]][:partial] == 0 and
        ref_annotated[syn_val[:homology][:hits][0]][:score] > syn_val[:homology][:score]

      ref_annotated[syn_val[:homology][:hits][0]] = {
        key: prot,
        pId: syn_val[:homology][:pId],
        score: syn_val[:homology][:score],
        cov_query: syn_val[:homology][:cov_query],
        cov_subject: syn_val[:homology][:cov_subject],
        assert_cutoff: syn_val[:homology][:assert_cutoff],
        length: syn_val[:homology][:length][0],
        partial: (syn_val[:partial] ? 1 : 0)
      }

      # ref_annotated[syn_val[:homology][:hits][0]] = {
      #   key: prot,
      #   pId: syn_val[:homology][:pId],
      #   cov_query: syn_val[:homology][:cov_query],
      #   cov_subject: syn_val[:homology][:cov_subject],
      #   assert_cutoff: syn_val[:homology][:assert_cutoff],
      #   length: syn_val[:homology][:length][0],
      #   partial: (syn_val[:partial] ? 1 : 0)
      # }

    end

    # print ref_annotated
    query_lengths = get_fasta_length query_prot_file

    @refgenome.coding_seq.each do |ref_k, ref_v|
      gene = ""
      coverage_ref = ""
      coverage_query = ""
      query_length = ""
      pId = ""

      if ref_annotated[ref_v[:protId]] != nil

        if ref_annotated[ref_v[:protId]][:pId] >= @min_pid and
          ref_annotated[ref_v[:protId]][:cov_query] >= @min_cov and
          ref_annotated[ref_v[:protId]][:cov_subject] >= @min_cov

          gene = ref_annotated[ref_v[:protId]][:key]
          coverage_ref = ref_annotated[ref_v[:protId]][:cov_subject]
          query_length = query_lengths[ref_annotated[ref_v[:protId]][:key]]
          coverage_query = ref_annotated[ref_v[:protId]][:cov_query]
          pId = ref_annotated[ref_v[:protId]][:pId]
          partial = ref_annotated[ref_v[:protId]][:partial]
        end

      end

      _locus_tag = ref_v[:locustag] || ""
      _seq_len = "NA"
      # _seq_len = ref_v[:bioseq].seq.length.to_s if ! ref_v[:bioseq].nil?
      _seq_len = ref_v[:length].to_s if ! ref_v[:length].nil?

      synteny_file.write(ref_v[:protId])
      synteny_file.write("\t"+_locus_tag)
      synteny_file.write("\t"+_seq_len)
      synteny_file.write("\t"+coverage_ref.to_s)
      synteny_file.write("\t"+pId.to_s)
      synteny_file.write("\t"+gene)
      synteny_file.write("\t"+query_length.to_s)
      synteny_file.write("\t"+coverage_query.to_s)
      synteny_file.write("\t"+partial.to_s)
      synteny_file.write("\n")

    end

    synteny_file.close

  end


end                             # end of Class
