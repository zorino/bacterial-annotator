# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maximilien1er@gmail.com
# date:    	15-02-24
# version: 	0.0.1
# licence:  	



class SequenceAnnotation

  attr_accessor :gbk, :coding_seq, :cds_file, :rna_file

  # Initialize then genbank file
  def initialize gbk_file, outdir

    @gbk_file = gbk_file
    if ! File.exists? @gbk_file
      fetch_ncbi_genome(@gbk_file, outdir)
      @gbk_file = "#{outdir}/#{gbk_file}.gbk"
      # @gbk_file += ".gbk"
    end

    flat_gbk = Bio::FlatFile.auto(@gbk_file)

    # Check if gbk is valid
    if flat_gbk.dbclass != Bio::GenBank
      abort "Aborting : The input #{@gbk_file} is not a valid genbank file !"
    else
      @gbk = flat_gbk.next_entry
    end

    @bioseq = @gbk.to_biosequence

  end


  # Prepare CDS/proteins
  def get_cds

    if @coding_seq == nil

      @coding_seq = {}

      # Iterate over each CDS
      @gbk.each_cds do |ft|
        ftH = ft.to_hash
        loc = ft.locations
        gene = []
        product = []
        protId = ""
        if ftH.has_key? "pseudo"
          next
        end
        gene = ftH["gene"] if !ftH["gene"].nil?
        product = ftH["product"] if !ftH["product"].nil?
        protId = ftH["protein_id"][0] if !ftH["protein_id"].nil?
        locustag = ftH["locus_tag"][0] if !ftH["locus_tag"].nil?

        dna = get_DNA(ft,@bioseq)
        pep = dna.translate
        pepBioSeq = Bio::Sequence.auto(pep)
        dnaBioSeq = Bio::Sequence.auto(dna)

        if protId.strip == ""
          protId = locustag
        end

        @coding_seq[protId] = {
          protId: protId,
          location: loc,
          locustag: locustag,
          gene: gene[0],
          product: product[0],
          bioseq: pepBioSeq,
          bioseq_gene: dnaBioSeq,
          bioseq_len: pepBioSeq.length
        }
      end

    end

    @coding_seq

  end

  # Prepare rRNA tRNA
  def get_rna

    if @rna_seq == nil

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
        product = ""
        product = ftH["product"][0] if !ftH["product"].nil?
        locustag = ftH["locus_tag"][0] if !ftH["locus_tag"].nil?

        # puts "#{@accession}\t#{seqBeg}\t#{seqEnd}\t#{strand}\t#{protId}\t#{locustag}\t#{gene[0]}\t#{product[0]}"
        dna = get_DNA(ft,@bioseq)
        dnaBioSeq = Bio::Sequence.auto(dna)

        @rna_seq[locustag] = {
          type: ft.feature.to_s,
          location: loc,
          locustag: locustag,
          product: product,
          bioseq_gene: dnaBioSeq
        }

      end

    end

    @rna_seq

  end


  # Print CDS to files
  # RETURN : cds_file path
  def write_cds_to_file outdir

    cds_file = "#{@gbk.accession}.pep"
    dna_file = "#{@gbk.accession}.dna"

    if @coding_seq == nil
      get_cds
    end

    dna_out = File.open("#{outdir}/#{dna_file}", "w")
    File.open("#{outdir}/#{cds_file}", "w") do |fwrite|
      @coding_seq.each_key do |k|
        seqout = @coding_seq[k][:bioseq].output_fasta("#{k}",60)
        seqout_dna = @coding_seq[k][:bioseq_gene].output_fasta("#{k}",60)
        fwrite.write(seqout)
        dna_out.write(seqout_dna)
      end
    end
    dna_out.close

    @cds_file = "#{outdir}/" + cds_file

  end

  # Print RNA to files
  # RETURN : rna_file path
  def write_rna_to_file outdir

    rna_file = "#{@gbk.accession}.rna"

    if @rna_seq == nil
      get_rna
    end

    File.open("#{outdir}/#{rna_file}", "w") do |fwrite|
      @rna_seq.each_key do |k|
        seqout_dna = @rna_seq[k][:bioseq_gene].output_fasta("#{k}|#{@rna_seq[k][:type]}|#{@rna_seq[k][:product]}",60)
        fwrite.write(seqout_dna)
      end
    end

    @rna_file = "#{outdir}/" + rna_file

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

      next if ! synteny_prot.has_key? prot_id or
        ! synteny_prot[prot_id].has_key? :homology

      # puts "#{annotations.keys}"
      if annotations.has_key? synteny_prot[prot_id][:homology][:hits][0]
        hit = annotations[synteny_prot[prot_id][:homology][:hits][0]]
        # puts hit
      else
        puts "no hit for #{prot_id}"
        next
      end

      # hit = annotations[synteny_prot[prot_id][:homology][:hits][0]]

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

      cds.qualifiers = ftArray

    end


  end


  # add annotation to a genbank file produced by prodigal
  def add_annotations annotations, mode, reference_locus=nil

    # nb_of_added_ft = 0
    i = 0

    fdebug = File.open("debug-add-annotation.txt","w")

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

          fdebug.write(hit)
          fdebug.write("\n")

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

      @gbk.features.each_with_index do |ft, ft_index|

        sorted_annotations.each do |k,v|

          next if annotations_done.has_key? k

          if v[:query_location][0][0] < ft.locations[0].from

            if v[:subject_location][0][0] > v[:subject_location][0][1]
              location = "complement(#{v[:query_location][0][0]}..#{v[:query_location][0][1]})"
            else
              location = "#{v[:query_location][0][0]}..#{v[:query_location][0][1]}"
            end

            feature = Bio::Feature.new(v[:feature][0],location)
            feature.qualifiers.push(Bio::Feature::Qualifier.new('product',v[:product][0])) if ! v[:product][0].nil? or v[:product][0] != ""
            new_features[ft_index] = feature
            annotations_done[k] = 1
            break

          end

        end

      end

      new_features.each do |k,v|
        @gbk.features.insert(k,v)
      end

    end

    fdebug.close

  end


  def save_genbank_to_file outdir

    File.open("#{outdir}/#{@gbk.definition}.gbk", "w") do |f|
      f.write(@gbk.to_biosequence.output(:genbank))
    end

  end

  ###################
  # Private Methods #
  ###################

  # Fct: Get dna sequence
  def get_DNA (cds, seq)
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
  def fetch_ncbi_genome refgenome_id, outdir
    Bio::NCBI.default_email = 'default@default.com'
    ncbi = Bio::NCBI::REST.new
    genbankstring = ncbi.efetch(refgenome_id, {"db"=>'nucleotide', "rettype"=>'gb'})
    File.open("#{outdir}/#{refgenome_id}.gbk", "w") do |f|
      f.write(genbankstring)
    end
  end

  private :fetch_ncbi_genome, :get_DNA


end                             # end of Class

