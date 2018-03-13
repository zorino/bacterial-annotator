# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maximilien1er@gmail.com
# date:    	15-02-24
# version: 	0.0.1
# licence:  	



class SequenceAnnotation

  attr_accessor :gbk, :coding_seq, :cds_file, :rna_file

  # Initialize then genbank file
  def initialize root, outdir, file_ref, type

    @root = root
    @outdir = outdir
    @coding_seq = {}
    @rna_seq = {}

    case type
    when "refGbk"
      # reference genome use for annotation
      reference_gbk file_ref
    when "db"
      # reference database use for annotation
      reference_db file_ref
    when "fasta"
      # single fasta database for annotation (completion)
      single_fasta file_ref
    when "newGbk"
      # new genbank holder to be annotated
      new_gbk file_ref
    end

  end


  # Use a MERGEM database to get annotation from it
  def reference_db dir

    abort "Aborting: Can't find MERGEM db direcotry" if ! File.exists? dir

    @cds_file = "#{dir}/cds.dmnd"
    @rna_file = "#{dir}/rnas.fasta"

    File.open("#{dir}/cds.txt") do |f|
      while l = f.gets
        lA = l.chomp.split(" ")
        @coding_seq[lA[0].gsub(">","")] = {
          protId: lA[0].gsub(">",""),
          location: nil,
          product: lA[1],
        }
      end
    end

    File.open("#{dir}/rnas.txt") do |f|
      while l = f.gets
        lA = l.chomp.split(" ")
        @rna_seq[lA[0].gsub(">","")] = {
          protId: lA[0].gsub(">",""),
          location: nil,
          product: lA[1],
        }
      end
    end

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

  # Use a Genbank Reference and read annotation from it
  def single_fasta fasta_file

    return "" if ! File.exists? fasta_file

    File.open(fasta_file, "r") do |dbfile|

      while l=dbfile.gets

        if l[0] == ">"

          lA = l.chomp.split("|")

          if lA.length > 1      # refseq, ncbi, trembl, swissprot

            key_gi = l.split(" ")[0][1..-1]
            product_long = lA[-1]

            organism = ""
            product = ""
            db_source = "[DBSource]"

            if product_long.scan(/|/).count >= 5 # FROM BIORUBY SCRIPTS
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


          else                  # mergem


          end

          @coding_seq[key_gi] = { product: product,
                                  org: org,
                                  prot_id: prot_id,
                                  db_source: db_source }

        end

      end

    end

  end


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


  # Prepare CDS/proteins
  def get_cds

    if @coding_seq.empty?

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

    if @rna_seq.empty?

      @rna_seq = {}
      @gbk.features do |ft|

        next if ! ft.feature.to_s.include? "rRNA"

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

