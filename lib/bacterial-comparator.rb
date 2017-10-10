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

class BacterialComparator

  attr_reader :genomes_list, :stats

  # Initialize BacterialAnnotator
  # options[:input], options[:refgenome], ROOT, options[:outdir], options)
  def initialize options, root

    @root = root
    @outdir = options[:outdir]
    Dir.mkdir(@outdir) if ! Dir.exists? @outdir
    @genomes_list = options[:genomes_list]
    @proc = options[:proc].to_i
    @phylo_nb_genes = options[:phylo_nb_genes]
    min_cov = options[:min_cov].to_f
    min_pid = options[:pidentity].to_f
    if min_cov > 1
      min_cov = min_cov/100
    end
    if min_pid > 1
      min_pid = min_pid/100
    end

    @ref_prot = get_ref_prot
    @synteny = read_prot_synteny
    @stats = extract_syntenic_fasta min_cov, min_pid

  end

  def read_prot_synteny

    print "# Reading genome synteny files - from genome annotations.."
    start_time = Time.now
    synteny = {}
    @genomes_list.each do |g|
      genome_synteny = []
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
    print "done (#{c_time})\n"

    synteny

  end

  def get_ref_prot
    ref_prot = []
    pep_file = Dir["#{@genomes_list[0]}/*.pep"]
    flatfile = Bio::FlatFile.auto("#{pep_file[0]}")
    flatfile.each_entry do |entry|
      ref_prot << entry.definition.split(" ")[0]
    end
    flatfile.close
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

    pep_out_dir = "./#{@outdir}/align-genes-pep"

    ref_proteins = load_genome_cds(Dir["#{@genomes_list[0]}/*.pep"][0])
    synteny_list.each do |k,v|
      pep_out = File.open(pep_out_dir+"/#{k}.pep", "w")
      pep_out.write(ref_proteins[k])
      pep_out.close
    end

    @genomes_list.each_with_index do |g,i|

      genome_proteins = load_genome_cds("#{g}/Proteins.fa")
      synteny_list.each do |k,v|
        pep_out = File.open(pep_out_dir+"/#{k}.pep", "a")
        pep_out.write(genome_proteins[v[i][:query_prot]])
        pep_out.close
      end

    end

    dna_out_dir = "./#{@outdir}/align-genes-dna"
    ref_genes = load_genome_cds(Dir["#{@genomes_list[0]}/*.dna"][0])
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

    print "# Extracting Proteins and Genes multifasta.."
    start_time = Time.now

    nb_of_syntenic = 0
    stats = {}
    stats[:syntenic] = []
    fout = File.open("#{@outdir}/cds-synteny.tsv", "w")
    fout.write("Gene\t"+@genomes_list.join("\t")+"\n")

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

    pep_out_dir = "./#{@outdir}/align-genes-pep"
    dna_out_dir = "./#{@outdir}/align-genes-dna"
    Dir.mkdir(pep_out_dir) if ! Dir.exists? pep_out_dir
    Dir.mkdir(dna_out_dir) if ! Dir.exists? dna_out_dir

    synteny_list = to_build_multifasta.each_slice((to_build_multifasta.length/@proc)+1).to_a

    Parallel.map(synteny_list, in_processes: @proc) { |list|
      build_multifasta list
    }

    end_time = Time.now
    c_time = Helper.sec2str(end_time-start_time)
    print "done (#{c_time})\n"

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

    print "# Sequence alignments - conserved single proteins a.a. (MAFFT).."
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
    print "done (#{c_time})\n"

    # FIXME ugly hack to find out the reference genome
    ref_id = Dir["#{ori_dir}/#{@genomes_list[0]}/*.pep"][0].split('/')[-1].gsub(".pep","")

    concat_alignments "align-genes-pep.all.fasta", ref_id

    Dir.chdir(ori_dir)

  end

  def mafft_align_all_dna

    print "# Sequence alignments - conserved single genes dna (MAFFT).."
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

    # ugly hack to find out the reference genome
    ref_id = Dir["#{ori_dir}/#{@genomes_list[0]}/*.pep"][0].split('/')[-1].gsub(".pep","")

    end_time = Time.now
    c_time = Helper.sec2str(end_time-start_time)
    print "done (#{c_time})\n"

    concat_alignments "align-genes-dna.all.fasta", ref_id

    Dir.chdir(ori_dir)

  end


  def concat_alignments outfile, ref_id

    if File.exists?("../#{outfile}") and File.size("../#{outfile}") > 0
      puts "..Alignment concatenated file already exists, skipping."
      return
    end

    fout = File.open("../#{outfile}", "w")

    seq = ""
    Dir["*.aln"].each do |f|
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
      Dir["*.aln"].each do |f|
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

  def mafft_aln aln_opt

    if aln_opt == "both"
      mafft_align_all_pep
      mafft_align_all_dna
    elsif aln_opt == "prot"
      mafft_align_all_pep
    elsif aln_opt == "dna"
      mafft_align_all_dna
    end

  end


  def raxml_tree_dna bt
    print "# Genes DNA tree creation (RAXML).."
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
    print "done (#{c_time})\n"
  end

  def raxml_tree_pep bt
    print "# Proteins AA tree creation (RAXML).."
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
    print "done (#{c_time})\n"
  end


  def raxml_tree aln_opt, bt

    if aln_opt == "both"
      raxml_tree_dna bt
      raxml_tree_pep bt
    elsif aln_opt == "prot"
      raxml_tree_pep bt
    elsif aln_opt == "dna"
      raxml_tree_dna bt
    end

  end

end                             # end of Class
