# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maximilien1er@gmail.com
# review:  	
# date:    	15-02-24
# version: 	0.0.1
# licence:  	

require 'bio'
require 'fileutils'

require 'bacterial-annotator/sequence-fasta'
require 'bacterial-annotator/sequence-annotation'
require 'bacterial-annotator/sequence-synteny'
require 'helper'

class BacterialAnnotator

  # Initialize BacterialAnnotator
  # options, ROOT (path)
  def initialize options, root

    @root = root
    @options = options

    abort if ! @options.has_key? :input

    @minlength = @options[:minlength].to_i
    @options[:minlength] = @options[:minlength].to_i
    @options[:pidentity] = @options[:pidentity].to_f
    @options[:pidentity] = @options[:pidentity] * 100 if @options[:pidentity] <= 1.00
    @options[:pcoverage] = @options[:pcoverage].to_f
    @options[:pcoverage] = @options[:pcoverage] / 100 if @options[:pcoverage] > 1.00

    if ! @options.has_key? :name
      @options[:name] = @options[:input].gsub(/.fasta|.fa|.fna/,"")
    end

    if File.exists? (@options[:outdir])
      if ! options.has_key? :force
        abort "Output directory already exist ! Choose another one or use -f to overwrite"
      else
        puts "Overwriting output directory #{@options[:outdir]}"
        FileUtils.remove_dir(@options[:outdir], :force=>true)
      end
    end
    Dir.mkdir(@options[:outdir])

    @query_fasta = SequenceFasta.new(@root,
                                     @options[:outdir],
                                     @options[:input],
                                     @options[:meta])

    @with_refence_genome = false
    @with_db = false
    if @options.has_key? :refgenome
      @with_refence_genome = true
      @ref_genome = SequenceAnnotation.new(@root,
                                           @options[:outdir],
                                           @options[:refgenome],
                                           "refGbk")
    elsif @options[:mergem]
      @with_db = true
      @ref_genome = SequenceAnnotation.new(@root,
                                           @options[:outdir],
                                           @options[:mergem],
                                           "db")
    end

    @with_external_db = false
    @with_external_db = true if @options.has_key? :external_db

    @prot_synteny = nil
    @annotation_stats = {
      by_contigs: {},
      annotated_cds: 0,
      flagged_cds: [],
      total_cds: 0,
      foreign_contigs: [],
      synteny_contigs: [],
      short_contigs: []
    }

    @contig_foreign_cds = {}

    @contig_annotations = {}

    @contig_annotations_externaldb = {}

    @contig_annotations_cds = {}

  end                           # end of method


  # run_alignment of reference genome proteins and the query
  def run_annotation

    prepare_files_for_annotation

    # process reference genome synteny
    if @with_refence_genome        # Annotation with the Reference Genome

      @prot_synteny_refgenome = run_reference_synteny_prot

      # iterate over each contig
      #     discard short contig
      #     cumulate statistics of homolog CDS
      @query_fasta.annotation_files[:contigs].each_with_index do |contig, contig_index|

        # Skip short contigs
        if @query_fasta.annotation_files[:contigs_length][contig_index] < @minlength
          @annotation_stats[:short_contigs] << contig
          next
        end

        remaining_cds = cumulate_annotation_stats_reference contig

        if remaining_cds != []
          @contig_foreign_cds[contig] = remaining_cds
        end

      end

      # dump foreign proteins to file
      foreign_cds_file = dump_cds

      # dump reference CDS synteny to file
      dump_ref_synteny_to_file

      # run RNA annotation
      @rna_synteny = SequenceSynteny.new(@root,
                                         @options[:outdir],
                                         @query_fasta.fasta_file,
                                         @ref_genome.rna_file,
                                         "RNA-Ref",
                                         @options[:pidentity],
                                         @options[:pcoverage],
                                         "dna")

      print "# Running alignment with Reference Genome RNA (blat).."
      start_time = Time.now
      @rna_synteny.run_blat
      end_time = Time.now
      c_time = Helper.sec2str(end_time-start_time)
      print "done (#{c_time})\n"

      # # takes too long
      # print "# Running alignment with Reference Genome RNA (fasta36).."
      # start_time = Time.now
      # @rna_synteny.run_fasta36
      # end_time = Time.now
      # c_time = Helper.sec2str(end_time-start_time)
      # print "done (#{c_time})\n"

      @rna_synteny.extract_hits_dna :rna
      @contig_annotations_rna = {}
      @query_fasta.annotation_files[:contigs].each_with_index do |contig, contig_index|
        @contig_annotations_rna[contig] = @rna_synteny.get_annotation_for_contig contig
      end


    elsif @with_db

      @prot_synteny_refgenome = run_mergem_synteny_prot
      # iterate over each contig
      #     discard short contig
      #     cumulate statistics of homolog CDS
      @query_fasta.annotation_files[:contigs].each_with_index do |contig, contig_index|

        # Skip short contigs
        if @query_fasta.annotation_files[:contigs_length][contig_index] < @minlength
          @annotation_stats[:short_contigs] << contig
          next
        end

        remaining_cds = cumulate_annotation_stats_reference contig

        if remaining_cds != []
          @contig_foreign_cds[contig] = remaining_cds
        end

      end

      # dump foreign proteins to file
      foreign_cds_file = dump_cds

      # dump reference CDS synteny to file
      dump_ref_synteny_to_file


    else                        # no reference genome

      # no reference genome .. will process all the CDS as foreign for the external db
      foreign_cds_file = @query_fasta.annotation_files[:proteins]

    end

    # Finishing annotation for foreign proteins
    finish_annotation foreign_cds_file

    # Parse annotations to genbank files
    parse_genbank_files

    print "# Printing Statistics.."
    print_stats "#{@options[:outdir]}"
    print "done\n"

  end                           # end of method


  # Prepare files for the annotation
  # Will run prodigal on the query and prepare reference genome files
  def prepare_files_for_annotation
    print "# Running Prodigal on your genome.."
    start_time = Time.now
    @query_fasta.run_prodigal
    end_time = Time.now
    c_time = Helper.sec2str(end_time - start_time)
    print "done (#{c_time})\n"
  end                           # end of method


  def run_mergem_synteny_prot


    ref_synteny_prot = SequenceSynteny.new(@root,
                                           @options[:outdir],
                                           @query_fasta.annotation_files[:proteins],
                                           @ref_genome.cds_file,
                                           "Prot-Ref",
                                           @options[:pidentity],
                                           @options[:pcoverage],
                                           "prot")

    print "# Running alignment with Reference Genome CDS (diamond).."
    start_time = Time.now
    ref_synteny_prot.run_diamond
    end_time = Time.now
    c_time = Helper.sec2str(end_time - start_time)
    print "done (#{c_time})\n"

    ref_synteny_prot.extract_hits :refgenome

    ref_synteny_prot.query_sequences.each do |k,v|
      if v.has_key? :homology
        @contig_annotations_cds[v[:contig]] = [] if ! @contig_annotations_cds.has_key? v[:contig]
        @contig_annotations_cds[v[:contig]] << k
      end
    end

    ref_synteny_prot


  end



  def run_reference_synteny_prot

    ref_synteny_prot = SequenceSynteny.new(@root,
                                           @options[:outdir],
                                           @query_fasta.annotation_files[:proteins],
                                           @ref_genome.cds_file,
                                           "Prot-Ref",
                                           @options[:pidentity],
                                           @options[:pcoverage],
                                           "prot")

    print "# Running alignment with Reference Genome CDS (diamond).."
    start_time = Time.now
    ref_synteny_prot.run_diamond
    end_time = Time.now
    c_time = Helper.sec2str(end_time - start_time)
    print "done (#{c_time})\n"

    # print "# Running alignment with Reference Genome CDS (blat).."
    # start_time = Time.now
    # ref_synteny_prot.run_blat
    # end_time = Time.now
    # c_time = Helper.sec2str(end_time - start_time)
    # print "done (#{c_time})\n"

    # print "# Running alignment with Reference Genome CDS (fasta36).."
    # start_time = Time.now
    # ref_synteny_prot.run_fasta36
    # end_time = Time.now
    # c_time = Helper.sec2str(end_time - start_time)
    # print "done (#{c_time})\n"

    ref_synteny_prot.extract_hits :refgenome

    ref_synteny_prot.query_sequences.each do |k,v|
      if v.has_key? :homology
        @contig_annotations_cds[v[:contig]] = [] if ! @contig_annotations_cds.has_key? v[:contig]
        @contig_annotations_cds[v[:contig]] << k
      end
    end

    ref_synteny_prot

  end


  # Finishing the annotation of the remaining CDS
  def finish_annotation remaining_cds_file

    # only finish the annotation with an external DB
    if @options.has_key? :external_db	# from an external DB

      db_file = @options[:external_db]
      ref_cds = SequenceAnnotation.new(@root,
                                       @options[:outdir],
                                       db_file,
                                      "fasta")

      # ref_cds = extract_externaldb_prot_info db_file

      @externaldb_synteny = SequenceSynteny.new(@root,
                                                @options[:outdir],
                                                remaining_cds_file,
                                                db_file,
                                                "Prot-ExternalDB",
                                                @options[:pidentity],
                                                @options[:pcoverage],
                                                "prot")

      print "# Running BLAT alignment with External Database.."
      start_time = Time.now
      @externaldb_synteny.run_blat
      end_time = Time.now
      c_time = Helper.sec2str(end_time-start_time)
      print "done (#{c_time})\n"
      @externaldb_synteny.extract_hits :externaldb

      @externaldb_synteny.query_sequences.each do |k, v|

        contig_of_protein = k.split("_")[0..-2].join("_")

        if ! @contig_annotations_externaldb.has_key? contig_of_protein
          @contig_annotations_externaldb[contig_of_protein] = {}
        end

        next if ! v.has_key? :homology

        if ! @contig_annotations_cds.has_key? contig_of_protein
          @contig_annotations_cds[contig_of_protein] = []
        end
        @contig_annotations_cds[contig_of_protein] << k

        hit_gi = v[:homology][:hits][0]

        # note = "Protein homology (#{v[:pId]}% identity) with gi:#{hit_gi}"
        cov_query = (v[:homology][:cov_query]*100).round(2)
        cov_subject = (v[:homology][:cov_subject]*100).round(2)
        note = "Protein homology (AA identity: #{v[:homology][:pId]}%; coverage (q,s): #{cov_query}%,#{cov_subject}%) with #{ref_cds.coding_seq[hit_gi][:prot_id]}"
        inference = "similar to AA sequence:#{ref_cds.coding_seq[hit_gi][:db_source]}:#{ref_cds.coding_seq[hit_gi][:prot_id]}"

        if ref_cds.coding_seq[hit_gi][:org] != ""
          note +=  " from #{ref_cds.coding_seq[hit_gi][:org]}"
        end

        @contig_annotations_externaldb[contig_of_protein][v[:homology][:hits][0]] = {
          product: ref_cds.coding_seq[hit_gi][:product],
          feature: "cds",
          gene: nil,
          prot_id: ref_cds.coding_seq[hit_gi][:prot_id],
          locustag: nil,
          note: note,
          inference: inference
        }

      end

    end

  end                           # end of method


  # parse all genbank files
  def parse_genbank_files

    print "# Parsing annotation into genbank files.."
    start_time = Time.now
    @contig_annotations_cds.each do |contig, contig_prots|

      gbk_path = @query_fasta.annotation_files[:gbk_path]
      gbk_to_annotate = SequenceAnnotation.new(@root,
                                               "#{gbk_path}",
                                               "#{gbk_path}/#{contig}.gbk",
                                               "newGbk")

      if @with_external_db and @with_refence_genome
        gbk_to_annotate.add_annotation_ref_synteny_prot(
          (@prot_synteny_refgenome.query_sequences.merge(@externaldb_synteny.query_sequences)),
          @contig_annotations_externaldb[contig].merge(@ref_genome.coding_seq),
          (File.basename @options[:refgenome]).gsub(/.gb.*/,"")
        )
      elsif @with_external_db
        gbk_to_annotate.add_annotation_ref_synteny_prot(
          @externaldb_synteny.query_sequences,
          @contig_annotations_externaldb[contig]
        )
      elsif @with_db
        gbk_to_annotate.add_annotation_ref_synteny_prot(
          @prot_synteny_refgenome.query_sequences,
          @ref_genome.coding_seq
        )
      else
        gbk_to_annotate.add_annotation_ref_synteny_prot(
          @prot_synteny_refgenome.query_sequences,
          @ref_genome.coding_seq,
          (File.basename @options[:refgenome]).gsub(/.gb.*/,"")
        )
      end

      if @contig_annotations_rna and @contig_annotations_rna.has_key? contig
        # puts "RNA annotation"
        gbk_to_annotate.add_annotations @contig_annotations_rna[contig], "new"
      end

      gbk_to_annotate.save_genbank_to_file

    end
    end_time = Time.now
    c_time = Helper.sec2str(end_time-start_time)
    print "done (#{c_time})\n"
  end                           # end of method


  # cumulate the stats for the synteny
  # return : unannotated cds array
  # def cumulate_annotation_stats_reference contig, contig_prots_ann
  def cumulate_annotation_stats_reference contig

    remaining_cds = []
    contig_prots = @query_fasta.annotation_files[:prot_ids_by_contig][contig]

    @annotation_stats[:total_cds] += contig_prots.length if contig_prots

    # count contig as foreign if no cds homolog in reference genome
    if @contig_annotations_cds.has_key? contig and
       @contig_annotations_cds[contig].length > 0
      @annotation_stats[:synteny_contigs] << contig
    else
      @annotation_stats[:foreign_contigs] << contig
      return
    end

    contig_prots.each do |prot|

      if @contig_annotations_cds[contig].include? prot

        if @prot_synteny_refgenome.query_sequences[prot].has_key? :homology and
           @prot_synteny_refgenome.query_sequences[prot][:homology][:hits].length > 0

          assert_sum = @prot_synteny_refgenome.query_sequences[prot][:homology][:assert_cutoff].inject(:+)
          if assert_sum > 2
            @annotation_stats[:annotated_cds] += 1
          else
            flag = "#{prot}"
            flag += "\t#{@prot_synteny_refgenome.query_sequences[prot][:homology][:assert_cutoff].join(',')}"
            flag += "\t#{@prot_synteny_refgenome.query_sequences[prot][:homology][:pId]}"
            flag += "\t#{(@prot_synteny_refgenome.query_sequences[prot][:homology][:cov_query]*100).round(2)}"
            flag += "\t#{(@prot_synteny_refgenome.query_sequences[prot][:homology][:cov_subject]*100).round(2)}"
            @annotation_stats[:flagged_cds] << flag
          end

        else

          puts "No " + prot

        end

      else

        remaining_cds << prot

      end

    end

    remaining_cds

  end                           # end of method


  # print statistics to file
  def print_stats file_dir

    file = file_dir + "/Annotation-Stats.txt"
    total_nb_contigs = @annotation_stats[:foreign_contigs].length +
                       @annotation_stats[:synteny_contigs].length +
                       @annotation_stats[:short_contigs].length
    p_contigs_annotated = @annotation_stats[:synteny_contigs].length.to_f/total_nb_contigs.to_f
    p_cds_annotated = @annotation_stats[:annotated_cds].to_f/@annotation_stats[:total_cds].to_f

    File.open(file, "w") do |fopen|

      fopen.write("#Contigs annotation based on reference genomes\n")
      fopen.write("Short Contigs (< #{@minlength}) :\t\t" + @annotation_stats[:short_contigs].length.to_s + "\n")
      fopen.write("Foreign Contigs :\t\t" + @annotation_stats[:foreign_contigs].length.to_s + "\n")
      fopen.write("Annotated Contigs :\t\t" + @annotation_stats[:synteny_contigs].length.to_s + "\n")
      fopen.write("Total Contigs :\t\t\t" + total_nb_contigs.to_s + "\n")
      fopen.write("% Contigs annotated :\t\t" + (p_contigs_annotated*100).round(2).to_s + "\n")
      fopen.write("\n")

      fopen.write("#CDS annotations based on reference genomes\n")
      fopen.write("Annotated CDS :\t\t\t" + @annotation_stats[:annotated_cds].to_s + "\n")
      fopen.write("Flagged CDS :\t\t\t" + @annotation_stats[:flagged_cds].length.to_s + "\n")
      fopen.write("Total CDS :\t\t\t" + @annotation_stats[:total_cds].to_s + "\n")
      fopen.write("% CDS annotated :\t\t" + (p_cds_annotated*100).round(2).to_s + "\n")
      fopen.write("\n")

    end

    file_flagged_cds = file_dir + "/Prot-flagged.tsv"
    File.open(file_flagged_cds, "w") do |fopen|
      fopen.write("CDS locus\tAssertion-CutOff\tAA Identity\tCovQuery(%)\tCovSubject(%)\n")
      @annotation_stats[:flagged_cds].each do |fcds|
        fopen.write("#{fcds}\n")
      end
    end

  end                           # end of method


  # dump cds to file for blast
  def dump_cds

    cds_outfile = File.open("#{@options[:outdir]}/Proteins-foreign.fa","w")
    foreign_cds = []
    @contig_foreign_cds.each_value do |v|
      foreign_cds.push(*v)
    end
    inprot = false
    File.open(@query_fasta.annotation_files[:proteins]) do |fprot|
      while l=fprot.gets
        if l[0] == ">"
          inprot = false
          prot_id = l.chomp.split(" ")[0].delete(">")
          if foreign_cds.include? prot_id
            inprot = true
            cds_outfile.write(l)
          end
        elsif inprot
          cds_outfile.write(l)
        end
      end
    end
    cds_outfile.close
    return "#{@options[:outdir]}/Proteins-foreign.fa"

  end                           # end of method


  # extract the information on protein from an externaldb
  def extract_externaldb_prot_info db

    # NCBI
    # >gi|103485499|ref|YP_615060.1| chromosomal replication initiation protein [Sphingopyxis alaskensis RB2256]
    # Swissprot
    # >sp|C7C422|BLAN1_KLEPN Beta-lactamase NDM-1 OS=Klebsiella pneumoniae GN=blaNDM-1 PE=1 SV=1
    # TrEMBL
    # >tr|E5KIY2|E5KIY2_ECOLX Beta-lactamase NDM-1 OS=Escherichia coli GN=blaNDM-1 PE=1 SV=1
    # MERGEM
    # >Genome_ID|location|Protein_ID|LocusTag|Gene|Protein_Product

    ref_cds = {}

    File.open(db, "r") do |dbfile|
      while l=dbfile.gets

        if l[0] == ">"

          lA = l.chomp.split("|")
          #key_gi = lA[1]
          key_gi = l.split(" ")[0][1..-1]
          product_long = lA[-1]

          organism = ""
          product = ""
          db_source = "[DBSource]"

          if product_long.scan(/|/).count >= 5 # MERGEM
            product = product_long
            db_source = "RefSeq"
          elsif product_long.include? " [" and product_long.include? "]" # NCBI
            organism = product_long[/\[.*?\]/]
            product = product_long.split(" [")[0].strip
          elsif product_long.include? "OS=" # Swissprot / TrEMBL
            product_tmp = product.split("OS=")
            organism = product_tmp[1].split(/[A-Z][A-Z]=/)[0].strip
            product = product_tmp[0].strip
          elsif product_long.include? "[A-Z][A-Z]=" # NCBI
            product = product_long.split(/[A-Z][A-Z]=/)[0].strip
          else
            product = product_long
          end

          org = organism.gsub("[","").gsub("]","")

          product.lstrip!
          prot_id = nil

          if key_gi.count("|") == 4
            if lA[2] == "ref"
              db_source = "RefSeq"
            end
            prot_id = lA[3]
          elsif key_gi.count("|") == 2
            if lA[0].include? == "sp" or
              lA[0].include? == "tr"
              db_source = "UniProtKB"
            end
            prot_id = lA[1]
          elsif key_gi.count("|") == 5
            db_source = "RefSeq"
            prot_id = lA[2]
          end

          ref_cds[key_gi] = {product: product, org: org, prot_id: prot_id, db_source: db_source}

        end

      end

    end                         # end of file reading

    ref_cds

  end                           # end of method


  # split fasta file to multiple fasta
  def split_remaining_cds_file file

    cds_files = []
    outdir = "#{@options[:outdir]}/Protein-foreign.split"

    Dir.mkdir(outdir) if ! Dir.exists? outdir

    iter = 0
    file_nb = 0
    fout = File.open("#{outdir}/ProtForeign.#{file_nb}.fa", "w")
    cds_files << "#{outdir}/ProtForeign.#{file_nb}.fa"

    File.open(file, "r") do |fopen|
      while l=fopen.gets
        if l[0] == ">"
          if iter > 9
            fout.close
            iter = 0
            file_nb += 1
            fout = File.open("#{outdir}/ProtForeign.#{file_nb}.fa", "w")
            cds_files << "#{outdir}/ProtForeign.#{file_nb}.fa"
          end
          iter += 1
        end
        fout.write(l)
      end
    end

    fout.close

    cds_files

  end                           # end of method

  # will dump reference CDS synteny to file
  def dump_ref_synteny_to_file

    # Iterate over each Ref protein and print syntheny
    synteny_file = File.open("#{@options[:outdir]}/Prot-Synteny.tsv","w")
    synteny_file.write("RefLocusTag\tRefProtID\tRefLength\tRefCoverage\tIdentity\tQueryGene\tQueryLength\tQueryCoverage\tQueryPartial\n")
    ref_annotated = {}

    @prot_synteny_refgenome.query_sequences.each do |prot, syn_val|
      next if ! syn_val.has_key? :homology
      next if syn_val[:homology][:assert_cutoff].inject(:+) < 3
      next if ref_annotated.has_key? syn_val[:homology][:hits][0] and ref_annotated[syn_val[:homology][:hits][0]][:partial] == 0
      ref_annotated[syn_val[:homology][:hits][0]] = {
        key: prot,
        pId: syn_val[:homology][:pId],
        cov_query: syn_val[:homology][:cov_query],
        cov_subject: syn_val[:homology][:cov_subject],
        assert_cutoff: syn_val[:homology][:assert_cutoff],
        length: syn_val[:homology][:length][0],
        partial: (syn_val[:partial] ? 1 : 0)
      }
    end

    @contig_annotations.each do |contig, prot_annotations|
      prot_annotations.each do |key,prot|
        ref_annotated[prot[:protId]] = {key: key, length: prot[:length], pId: prot[:pId]} if prot != nil
      end
    end

    @ref_genome.coding_seq.each do |ref_k, ref_v|

      gene = ""
      coverage_ref = ""
      coverage_query = ""
      query_length = ""
      pId = ""
      if ref_annotated[ref_v[:protId]] != nil
        gene = ref_annotated[ref_v[:protId]][:key]
        coverage_ref = ref_annotated[ref_v[:protId]][:cov_subject]
        query_length = @query_fasta.annotation_files[:prot_ids_length][gene]
        coverage_query = ref_annotated[ref_v[:protId]][:cov_query]
        pId = ref_annotated[ref_v[:protId]][:pId]
        partial = ref_annotated[ref_v[:protId]][:partial]
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

  private :dump_cds, :split_remaining_cds_file, :dump_ref_synteny_to_file

end                             # end of class
