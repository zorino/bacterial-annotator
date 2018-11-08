# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maximilien1er@gmail.com
# date:    	15-02-24
# version: 	0.0.1
# licence:  	

require 'json'
require 'zlib'
require 'pp'

class SequenceAnnotationV2

  attr_accessor :gbk, :fasta, :annotation_stats

  # Initialize a SequenceAnnotation
  def initialize root, options

    @options = options

    @root = root
    @outdir = options[:outdir]

    @coding_seq = {}
    @rna_seq = {}

    # TODO remove that shit -
    # uniformise pour toute les types d'annotations avec des classes
    @with_refence_genome = true
    @with_external_db = false
    @with_db = false

    fasta_sequence_init options
    stats_init

  end

  # Initialize a FastaSequence for annotation
  def fasta_sequence_init options

    @fasta = SequenceFastaV2.new(@root, options)
    print "# Running Prodigal on your genome.."
    start_time = Time.now
    @fasta.run_prodigal
    end_time = Time.now
    c_time = Helper.sec2str(end_time - start_time)
    print "done (#{c_time})\n"

  end

  # Initialize statistics of annotation
  def stats_init

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

  end


  # run annotation
  # args:
  #   aln_type : dna / prot
  #   db file : fasta file for alignment
  #   genes : hash of gene information
  def run_annotation aln_type, db_file, genes, software

    # seq_synteny = SequenceSynteny.new(@root,
    #                                   @options[:outdir],
    #                                   @fasta.annotation_files[:proteins],
    #                                   ref_genome.cds_file,
    #                                   "Prot-Ref",
    #                                   @options[:pidentity],
    #                                   @options[:pcoverage],
    #                                   "prot")

    seq_synteny = SequenceSyntenyV2.new(@root,
                                        @options,
                                        @fasta.annotation_files[:proteins],
                                        db_file,
                                        "Prot-Ref",
                                        aln_type,
                                        software)

    print "# Running #{aln_type} alignment using #{software}.."
    start_time = Time.now
    case software
    when "diamond"
      seq_synteny.run_diamond
    when "fasta36"
      seq_synteny.run_fasta36
    when "blat"
      seq_synteny.run_blat
    end
    end_time = Time.now
    c_time = Helper.sec2str(end_time - start_time)
    print "done (#{c_time})\n"

    seq_synteny.extract_hits :refgenome

    seq_synteny.query_sequences.each do |k,v|
      if v.has_key? :homology
        @contig_annotations_cds[v[:contig]] = [] if ! @contig_annotations_cds.has_key? v[:contig]
        @contig_annotations_cds[v[:contig]] << k
      end
    end

    pp(seq_synteny)

    @prot_synteny_refgenome = seq_synteny

  end



  # run annotation with reference genome
  def run_annotation_ref_genome ref_genome

    @ref_genome = ref_genome

    run_reference_synteny_prot ref_genome

    # iterate over each contig
    #     discard short contig
    #     cumulate statistics of homolog CDS
    @fasta.annotation_files[:contigs].each_with_index do |contig, contig_index|

      # Skip short contigs
      if @fasta.annotation_files[:contigs_length][contig_index] < @options[:minlength]
        @annotation_stats[:short_contigs] << contig
        next
      end

      remaining_cds = cumulate_annotation_stats_reference contig

      if remaining_cds != []
        @contig_foreign_cds[contig] = remaining_cds
      end

    end

    # dump foreign proteins to file
    foreign_cds_file = dump_foreign_cds

    # Finishing annotation for foreign proteins
    # finish_annotation foreign_cds_file

    # dump reference CDS synteny to file - TODO fix this
    dump_ref_synteny_to_file

    # Parse annotations to genbank files
    parse_genbank_files

    print "# Printing Statistics.."
    print_stats
    print "done\n"

  end


  # run cds synteny with reference genome
  def run_reference_synteny_prot ref_genome

    ref_synteny_prot = SequenceSynteny.new(@root,
                                           @options[:outdir],
                                           @fasta.annotation_files[:proteins],
                                           ref_genome.cds_file,
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

    @prot_synteny_refgenome = ref_synteny_prot

  end


  # output annotation to genbank files
  def parse_genbank_files

    print "# Parsing annotation into genbank files.."
    start_time = Time.now

    @contig_annotations_cds.each do |contig, contig_prots|

      gbk_path = @fasta.annotation_files[:gbk_path]
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


  end



  # cumulate the stats for the synteny
  # return : unannotated cds array
  # def cumulate_annotation_stats_reference contig, contig_prots_ann
  def cumulate_annotation_stats_reference contig

    remaining_cds = []
    contig_prots = @fasta.annotation_files[:prot_ids_by_contig][contig]

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
            remaining_cds << prot
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


  # dump cds to file for blast
  def dump_foreign_cds

    cds_outfile = File.open("#{@outdir}/Proteins-foreign.fa","w")
    foreign_cds = []
    @contig_foreign_cds.each_value do |v|
      foreign_cds.push(*v)
    end
    inprot = false
    File.open(@fasta.annotation_files[:proteins]) do |fprot|
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
    return "#{@outdir}/Proteins-foreign.fa"

  end                           # end of method


  # will dump reference CDS synteny to file
  def dump_ref_synteny_to_file

    # Iterate over each Ref protein and print syntheny
    synteny_file = File.open("#{@outdir}/Prot-Synteny.tsv","w")
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
        query_length = @fasta.annotation_files[:prot_ids_length][gene]
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


  # print statistics to file
  def print_stats

    file = @outdir + "/Annotation-Stats.txt"
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

    file_flagged_cds = @outdir + "/Prot-flagged.tsv"
    File.open(file_flagged_cds, "w") do |fopen|
      fopen.write("CDS locus\tAssertion-CutOff\tAA Identity\tCovQuery(%)\tCovSubject(%)\tNote\n")
      @annotation_stats[:flagged_cds].each do |fcds|
        fopen.write("#{fcds}\n")
      end
    end

  end                           # end of method



  ######### TODO - STILL USED


  # New Genbank Holder to add annotation to it
  def new_gbk gbk_file

    if ! File.exists? gbk_file
      fetch_ncbi_genome(gbk_file)
      gbk_file = "#{@outdir}/#{gbk_file}.gbk"
      # gbk_file += ".gbk"
    end

    flat_gbk = Bio::FlatFile.auto(gbk_file)

    # Check if gbk is valid
    if flat_gbk.dbclass != Bio::GenBank
      abort "Aborting : The input #{gbk_file} is not a valid genbank file !"
    else
      @gbk = flat_gbk.next_entry
    end

    @bioseq = @gbk.to_biosequence

  end


  # add annotation from reference prot synteny
  def add_annotation_ref_synteny_prot synteny_prot, annotations, ref_genome=nil

    contig = @gbk.definition

    prot_iterator = 0
    @gbk.features.each_with_index do |cds, ft_index|

      next if cds.feature != "CDS"

      prot_iterator+=1
      prot_id = contig+"_"+prot_iterator.to_s

      ftArray = []
      cds.qualifiers = []

      hit = nil

      if ! synteny_prot.has_key? prot_id or
         ! synteny_prot[prot_id].has_key? :homology or
         ! annotations.has_key? synteny_prot[prot_id][:homology][:hits][0] or
         synteny_prot[prot_id][:homology][:assert_cutoff].inject(:+) < 3

        qLocusTag = Bio::Feature::Qualifier.new('locus_tag', "#{prot_id}")
        qProd = Bio::Feature::Qualifier.new('product', "hypothetical protein")

        ftArray.push(qLocusTag)
        ftArray.push(qProd)

        if synteny_prot.has_key? prot_id and
          synteny_prot[prot_id].has_key? :homology and
          synteny_prot[prot_id][:homology][:assert_cutoff] == [1,1,0]
          hit = annotations[synteny_prot[prot_id][:homology][:hits][0]]
          qNote = Bio::Feature::Qualifier.new('note', "possible pseudo gene of #{hit[:locustag]} from #{ref_genome}")
          ftArray.push(qNote)
        end

      else

        hit = annotations[synteny_prot[prot_id][:homology][:hits][0]]

        if synteny_prot.has_key? prot_id

          locus, gene, product, note, inference = nil
          locus = hit[:locustag]
          gene = hit[:gene]
          product = hit[:product]
          note = hit[:note]
          inference = hit[:inference]
          pId = synteny_prot[prot_id][:homology][:pId]
          cov_query = (synteny_prot[prot_id][:homology][:cov_query]*100).round(2)
          cov_subject = (synteny_prot[prot_id][:homology][:cov_subject]*100).round(2)
          reference_prot_id = synteny_prot[prot_id][:homology][:hits][0]

          qLocusTag = Bio::Feature::Qualifier.new('locus_tag', "#{prot_id}")
          ftArray.push(qLocusTag)

          if gene != nil
            qGene = Bio::Feature::Qualifier.new('gene', gene)
            ftArray.push(qGene)
          end

          if product != nil
            qProd = Bio::Feature::Qualifier.new('product', product)
            ftArray.push(qProd)
          end

          # check if there is a reference genome.. reference_locus shouldn't be nil in that case
          if locus != nil

            qNote = Bio::Feature::Qualifier.new('note', "corresponds to #{locus} locus (AA identity: #{pId}%; coverage(q,s): #{cov_query}%,#{cov_subject}%) from #{ref_genome}")
            ftArray.push(qNote)

            db_source = "[DBSource]"
            if reference_prot_id.include? "_"
              db_source = "RefSeq"
            else
              db_source = "INSD"
            end
            qInference = Bio::Feature::Qualifier.new('inference', "similar to AA sequence:#{db_source}:#{reference_prot_id}")
            ftArray.push(qInference)

          end

          if note != nil
            qNote = Bio::Feature::Qualifier.new('note', note)
            ftArray.push(qNote)
          end

          if inference != nil
            qInference = Bio::Feature::Qualifier.new('inference', inference)
            ftArray.push(qInference)
          end

        end

      end

      cds.qualifiers = ftArray

    end


  end


  # add annotation to a genbank file produced by prodigal
  def add_annotations annotations, mode, reference_locus=nil

    # nb_of_added_ft = 0
    i = 0

    contig = @gbk.definition

    if mode == "inplace"

      # iterate through
      @gbk.features.each_with_index do |cds, ft_index|

        next if cds.feature != "CDS"

        ftArray = []
        cds.qualifiers = []

        i += 1
        prot_id = contig+"_"+i.to_s
        hit = nil

        if annotations.has_key? prot_id
          hit = annotations[prot_id]
        else
          puts "no hit for #{prot_id}"
          next
        end

        if hit != nil

          locus, gene, product, note = nil
          locus = hit[:locustag]
          gene = hit[:gene]
          product = hit[:product]
          note = hit[:note]
          pId = hit[:pId]

          if gene != nil
            qGene = Bio::Feature::Qualifier.new('gene', gene)
            ftArray.push(qGene)
          end

          if product != nil
            qProd = Bio::Feature::Qualifier.new('product', product)
            ftArray.push(qProd)
          end

          # check if there is a reference genome.. reference_locus shouldn't be nil in that case
          if locus != nil
            qNote = Bio::Feature::Qualifier.new('note', "corresponds to #{locus} locus (#{pId}% identity) from #{reference_locus.entry_id}")
            ftArray.push(qNote)
          end

          if note != nil
            qNote = Bio::Feature::Qualifier.new('note', note)
            ftArray.push(qNote)
          end

        end
        cds.qualifiers = ftArray

      end


    elsif mode == "new"

      sorted_annotations = annotations.sort_by { |k, v| v[:query_location][0][0] }

      new_features = {}
      annotations_done = {}
      gbk_features_len = @gbk.features.length

      @gbk.features.each_with_index do |ft, ft_index|

        sorted_annotations.each do |k,v|

          next if annotations_done.has_key? k

          if v[:query_location][0][0] < ft.locations[0].from or ft_index == gbk_features_len-1

            if v[:subject_location][0][0] > v[:subject_location][0][1]
              location = "complement(#{v[:query_location][0][0]}..#{v[:query_location][0][1]})"
            else
              location = "#{v[:query_location][0][0]}..#{v[:query_location][0][1]}"
            end

            feature = Bio::Feature.new(v[:feature][0],location)
            feature.qualifiers.push(Bio::Feature::Qualifier.new('product',v[:product][0])) if ! v[:product][0].nil? or v[:product][0] != ""
            if ft_index == gbk_features_len-1
              new_features[gbk_features_len] = feature
            else
              new_features[ft_index] = feature
            end

            annotations_done[k] = 1
            break

          end

        end

      end

      new_features.each do |k,v|
        @gbk.features.insert(k,v)
      end

    end

  end


  def save_genbank_to_file

    File.open("#{@outdir}/#{@gbk.definition}.gbk", "w") do |f|
      f.write(@gbk.to_biosequence.output(:genbank))
    end

  end

  ###################
  # Private Methods #
  ###################

  # Fct: Get dna sequence
  def get_DNA cds, seq
    loc = cds.locations
    sbeg = loc[0].from.to_i
    send = loc[0].to.to_i
    fasta = Bio::Sequence::NA.new(seq.subseq(sbeg,send))
    # position = "#{sbeg}..#{send}"
    if loc[0].strand == -1
      fasta.reverse_complement!
    end
    dna = Bio::Sequence.auto(fasta)
    return dna
  end


  # Fetch genbank genome from NCBI
  def fetch_ncbi_genome refgenome_id
    Bio::NCBI.default_email = 'default@default.com'
    ncbi = Bio::NCBI::REST.new
    genbankstring = ncbi.efetch(refgenome_id, {"db"=>'nucleotide', "rettype"=>'gb'})
    File.open("#{@outdir}/#{refgenome_id}.gbk", "w") do |f|
      f.write(genbankstring)
    end
  end

  private :fetch_ncbi_genome, :get_DNA


end                             # end of Class

