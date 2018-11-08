# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maximilien1er@gmail.com
# review:  	
# date:    	18-11-07
# version: 	0.0.1
# licence:  	

require 'json'
require 'zlib'
require 'pp'


class ReferenceGenome

  attr_accessor :gbk, :genes, :coding_seq, :rna_seq, :cds_file, :rna_file

  # Initialize then genbank file
  def initialize root, options

    @options = options

    @root = root
    @outdir = options[:outdir]
    @refGenome = options[:refgenome]

    @coding_seq = {}
    @rna_seq = {}
    @genes = {}

    reference_gbk @refGenome

  end


  # Use a Genbank Reference and read annotation from it
  def reference_gbk gbk_file

    puts "# Preparing reference genome files.."
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

    write_cds_to_file
    write_rna_to_file

  end

  # Print CDS to files
  # RETURN : cds_file path
  def write_cds_to_file

    cds_file = "#{@gbk.accession}.pep"
    dna_file = "#{@gbk.accession}.dna"

    if @coding_seq.empty?
      get_cds
    end

    dna_out = File.open("#{@outdir}/#{dna_file}", "w")
    File.open("#{@outdir}/#{cds_file}", "w") do |fwrite|
      @coding_seq.each_key do |k|
        seqout = @coding_seq[k][:bioseq].output_fasta("#{k}",60)
        seqout_dna = @coding_seq[k][:bioseq_gene].output_fasta("#{k}",60)
        fwrite.write(seqout)
        dna_out.write(seqout_dna)
      end
    end
    dna_out.close

    @cds_file = "#{@outdir}/" + cds_file

  end

  # Print RNA to files
  # RETURN : rna_file path
  def write_rna_to_file

    rna_file = "#{@gbk.accession}.rna"

    if @rna_seq.empty?
      get_rna
    end

    File.open("#{@outdir}/#{rna_file}", "w") do |fwrite|
      @rna_seq.each_key do |k|
        seqout_dna = @rna_seq[k][:bioseq_gene].output_fasta("#{k}|#{@rna_seq[k][:type]}|#{@rna_seq[k][:product]}",60)
        fwrite.write(seqout_dna)
      end
    end

    @rna_file = "#{@outdir}/" + rna_file

  end

  # Prepare CDS/proteins
  def get_cds

    if @coding_seq.empty?

      # Iterate over each CDS
      @gbk.each_cds do |ft|
        ftH = ft.to_hash
        loc = ft.locations
        gene = []
        product = []
        _id = ""
        if ftH.has_key? "pseudo"
          next
        end

        gene = ftH["gene"] if !ftH["gene"].nil?
        product = ftH["product"] if !ftH["product"].nil?
        _id = ftH["protein_id"][0] if !ftH["protein_id"].nil?
        locustag = ftH["locus_tag"][0] if !ftH["locus_tag"].nil?

        dna = get_DNA(ft,@bioseq)
        pep = dna.translate
        pepBioSeq = Bio::Sequence.auto(pep)
        dnaBioSeq = Bio::Sequence.auto(dna)

        if _id.strip == ""
          _id = locustag
        end

        @coding_seq[_id] = {
          protId: _id,
          location: loc,
          locustag: locustag,
          gene: gene[0],
          product: product[0],
          bioseq: pepBioSeq,
          bioseq_gene: dnaBioSeq,
          length: pepBioSeq.length
        }

        @genes[_id] = {
          id: _id,
          type: "CDS",
          location: loc,
          locustag: locustag,
          gene: gene[0],
          product: product[0],
          bioseq: pepBioSeq,
          length: pepBioSeq.length
        }

      end

    end

    @coding_seq

  end

  # Prepare rRNA tRNA
  def get_rna

    if @rna_seq.empty?

      @rna_seq = {}
      @gbk.features do |ft|

        next if ! ft.feature.to_s.include? "RNA"

        ftH = ft.to_hash
        loc = ft.locations
        # seqBeg = loc[0].from.to_s
        # seqEnd = loc[0].to.to_s
        # strand = loc[0].strand.to_s
        if ftH.has_key? "pseudo"
          next
        end
        # gene = ftH["gene"] if !ftH["gene"].nil?
        # protId = ftH["protein_id"][0] if !ftH["protein_id"].nil?

        gene = ""
        product = ""
        locustag = ""

        gene = ftH["gene"] if !ftH["gene"].nil?
        product = ftH["product"] if !ftH["product"].nil?
        locustag = ftH["locus_tag"][0] if !ftH["locus_tag"].nil?

        dna = get_DNA(ft,@bioseq)
        dnaBioSeq = Bio::Sequence.auto(dna)

        @rna_seq[locustag] = {
          type: ft.feature.to_s,
          location: loc,
          locustag: locustag,
          product: product,
          bioseq_gene: dnaBioSeq
        }

        @genes[locustag] = {
          id: locustag,
          type: ft.feature.to_s,
          location: loc,
          locustag: locustag,
          gene: gene[0],
          product: product[0],
          bioseq: dnaBioSeq,
          length: dnaBioSeq.length
        }

      end

    end

    @rna_seq

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


end
