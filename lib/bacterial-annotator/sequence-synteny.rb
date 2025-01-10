# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maximilien1er@gmail.com
# review:
# date:    	15-02-24
# version: 	0.0.1
# licence:

require 'json'
require 'zlib'

class SequenceSynteny

  attr_reader :query_file, :subject_file, :aln_hits, :query_sequences, :subject_sequences

  def initialize root, outdir, query_file, subject_file, name, pidentity, min_coverage, type

    @root = root
    @outdir = outdir
    @query_file = query_file
    @subject_file = subject_file
    @query_sequences = get_sequences(query_file)
    @subject_sequences = get_sequences(subject_file)

    @name = name
    @pidentity = pidentity
    @min_coverage = min_coverage
    @aln_file = nil
    @type = type

  end                           # end of initialize


  # get sequences name with length in hash
  def get_sequences raw_file

    sequences = {}

    if raw_file.include?(".dmnd")

      seq_info_file = raw_file.gsub(".dmnd",".json.gz")

      json_genes = {}
      Zlib::GzipReader.open(seq_info_file) {|gz|
        json_genes = JSON.parse(gz.read)
      }

      json_genes.each do |gene|

        sequences[gene["cluster_id"]] = {}
        sequences[gene["cluster_id"]][:length] = gene["consensus_length"].to_f
        sequences[gene["cluster_id"]][:conserved] = false
        sequences[gene["cluster_id"]][:contig] = gene["cluster_id"].split("_")[0..-2].join("_") if gene["cluster_id"].include? "_"

      end

    else

      seq_file = raw_file
      flat = Bio::FlatFile.auto("#{seq_file}")
      flat.each_entry do |s|
        s_name = s.definition.chomp.split(" ")[0]
        # puts s_name
        sequences[s_name] = {}
        properties = s.definition.chomp.split(";")
        partial = false
        if properties.length >= 2 and properties[1].include? "partial"
          partial = (properties[1].gsub("partial=","").include? '1')
        end
        sequences[s_name][:partial] = partial
        sequences[s_name][:length] = s.seq.length
        sequences[s_name][:conserved] = false

        sequences[s_name][:contig] = s_name.split("_")[0..-2].join("_") if s_name.include? "_"

      end

    end

    sequences

  end

  # run blat on proteins
  def run_blat
    base_cmd = "#{@root}/blat.linux -out=blast8 -minIdentity=#{@pidentity} > /dev/null 2>&1"
    if @type == "prot"
      system("#{base_cmd} -prot #{@subject_file} #{@query_file} #{@outdir}/#{@name}.blat8.tsv")
    else
      system("#{base_cmd} #{@subject_file} #{@query_file} #{@outdir}/#{@name}.blat8.tsv")
    end
    @aln_file = "#{@outdir}/#{@name}.blat8.tsv"
    # extract_hits
  end                           # end of method

  # run fasta36 on proteins
  def run_fasta36
    if @type == "prot"
      system("#{@root}/fasta36.linux -T 1 -b 3 -E 1e-40 -m 8 #{@query_file} #{@subject_file} > #{@outdir}/#{@name}.fasta36.tsv")
    else
      system("#{@root}/glsearch36.linux -T 1 -b 12 -E 1e-40 -m 8 #{@query_file} #{@subject_file} > #{@outdir}/#{@name}.fasta36.tsv")
    end
    @aln_file_fasta36 = "#{@outdir}/#{@name}.fasta36.tsv"
    # extract_hits
  end                           # end of method

  # run diamond on proteins
  def run_diamond
    if @type == "prot"
      if subject_file.include? ".dmnd"
        db_file = subject_file
      else
        system("#{@root}/diamond.linux makedb --db #{subject_file} --in #{subject_file} > /dev/null 2>&1")
        db_file = subject_file
      end
      system("#{@root}/diamond.linux blastp --masking none --db #{db_file} -q #{query_file} -o #{@outdir}/#{@name}.diamond.tsv -f 6 > /dev/null 2>&1")
    else
      # system("#{@root}/glsearch36.linux -b 3 -E 1e-25 -m 8 #{@subject_file} #{@query_file} > #{@outdir}/#{@name}.fasta36.tsv")
    end
    @aln_file = "#{@outdir}/#{@name}.diamond.tsv"
    # extract_hits
  end                           # end of method


  # Extract Hit from blast8 file and save it in hash
  # contig-0_1      ABJ71957.1      96.92   65      2       0       1       65      1       65      9.2e-31 131.0
  def extract_hits mode

    feature = ""
    File.open(@aln_file,"r") do |fread|
      while l = fread.gets

        lA = l.chomp!.split("\t")
        key = lA[0]

        # extraction of hit id depends on mode ..
        if mode == :refgenome
          hit = lA[1]
          feature = "cds"
        elsif mode == :externaldb
          # hit = lA[1].chomp.split("|")[3]
          hit = lA[1]
          feature = "cds"
          # decrease min. % identity for foreign proteins
          @pidentity = 60.0
        end

        # compute coverage based on sequences length
        cov_query = (lA[3].to_f/@query_sequences[key][:length]).round(2)
        cov_subject = (lA[3].to_f/@subject_sequences[hit][:length]).round(2)

        # assert cutoff on identity and coverage
        # 1 -> pass cutoff, 0 under cutoff
        assert_cutoff = [1,1,1]
        assert_cutoff[0] = 0 if lA[2].to_f < @pidentity
        assert_cutoff[1] = 0 if cov_query < @min_coverage
        assert_cutoff[2] = 0 if cov_subject < @min_coverage and @query_sequences[key][:partial] == false

        # first hit for query
        if ! @query_sequences[key].has_key? :homology
          @query_sequences[key][:conserved] = true
          # @subject_sequences[key][:conserved] = true
          @query_sequences[key][:homology] = {
            pId: lA[2].to_f.round(2),
            cov_query: cov_query,
            cov_subject: cov_subject,
            evalue: lA[10],
            score: lA[11].to_f,
            hits: [hit],
            length: [lA[3].to_i],
            query_location: [[lA[6].to_i,lA[7].to_i]],
            subject_location: [[lA[8].to_i,lA[9].to_i]],
            feature: feature,
            assert_cutoff: assert_cutoff
          }
          @subject_sequences[hit][:hits] = [key]

        # query already got at least 1 hit and new_score > last_score
        elsif lA[11].to_f > @query_sequences[key][:homology][:score]
          @query_sequences[key][:conserved] = true
          # @subject_sequences[key][:conserved] = true
          @query_sequences[key][:homology] = {
            pId: lA[2].to_f.round(2),
            cov_query: cov_query,
            cov_subject: cov_subject,
            evalue: lA[10],
            score: lA[11].to_f,
            hits: [hit],
            length: [lA[3].to_i],
            query_location: [[lA[6].to_i,lA[7].to_i]],
            subject_location: [[lA[8].to_i,lA[9].to_i]],
            feature: feature,
            assert_cutoff: assert_cutoff
          }
          @subject_sequences[hit][:hits] =  [key]

        # query already got at least 1 hit and score == last_score
        elsif lA[11].to_f == @query_sequences[key][:homology][:score]
          @query_sequences[key][:homology][:hits] << hit
          @query_sequences[key][:homology][:length] << lA[3].to_i
          @query_sequences[key][:homology][:query_location] << [lA[6].to_i,lA[7].to_i]
          @query_sequences[key][:homology][:subject_location] << [lA[8].to_i,lA[9].to_i]
          if @subject_sequences[hit].has_key? :hits
            @subject_sequences[hit][:hits] << key
          else
            @subject_sequences[hit][:hits] = [key]
          end
        end
      end
    end

  end                           # end of method


  # Extract Hit from blast8 file and save it in hash
  # prpa    PA0668.4|rRNA|23S       99.97   2891    1       0       705042  707932  1       2891    0.0e+00 5671.0
  def extract_hits_dna mode

    @aln_hits = {}
    feature = ""
    File.open(@aln_file,"r") do |fread|
      while l = fread.gets
        lA = l.chomp!.split("\t")
        key = lA[0]+"_"+lA[6]+"_"+lA[7]
        if mode == :rna
          hit_split = lA[1].chomp.split("|")
          hit = hit_split[0]
          feature = hit_split[1]
          product = hit_split[2]
        end
        next if lA[2].to_f < @pidentity
        if ! @aln_hits.has_key? key
          @aln_hits[key] = {
            pId: lA[2].to_f.round(2),
            # cov_query: (@query_sequences[key][:length]/lA[3].to_f).round(2),
            # cov_subject: (@subject_sequences[hit][:length]/lA[3].to_f).round(2),
            evalue: lA[10],
            score: lA[11].to_f,
            hits: [hit],
            product: [product],
            length: [lA[3].to_i],
            query_location: [[lA[6].to_i,lA[7].to_i]],
            subject_location: [[lA[8].to_i,lA[9].to_i]],
            feature: [feature]
          }
        elsif lA[11].to_f > @aln_hits[key][:score]
          @aln_hits[key] = {
            pId: lA[2].to_f.round(2),
            # cov_query: (@query_sequences[key][:length]/lA[3].to_f).round(2),
            # cov_subject: (@subject_sequences[hit][:length]/lA[3].to_f).round(2),
            evalue: lA[10],
            score: lA[11].to_f,
            hits: [hit],
            product: [product],
            length: [lA[3].to_i],
            query_location: [[lA[6].to_i,lA[7].to_i]],
            subject_location: [[lA[8].to_i,lA[9].to_i]],
            feature: [feature]
          }
        elsif lA[11].to_f == @aln_hits[key][:score]
          @aln_hits[key][:hits] << hit
          @aln_hits[key][:length] << lA[3].to_i
          @aln_hits[key][:query_location] << [lA[6].to_i,lA[7].to_i]
          @aln_hits[key][:subject_location] << [lA[8].to_i,lA[9].to_i]
          @aln_hits[key][:feature] << feature
          @aln_hits[key][:product] << product
        end
      end
    end

    prune_aln_hits @aln_hits

  end                           # end of method


  # Get the annotations for a contig for RerenceGenome
  def get_annotation_for_contig contig_to_annotate, prots_to_annotate=nil, ref_cds=nil

    annotations = {}

    if prots_to_annotate != nil

      # contig_to_annotate = prots_to_annotate[0].split("_")[0..-2].join("_")
      prots = []

      @aln_hits.each_key do |k|
        contig = k.split("_")[0..-2].join("_")
        if contig == contig_to_annotate
          prots << k
        end
      end

      # sorting the prot by their appearance in the contig
      prots.sort! { |a,b| a.split("_")[-1].to_i <=> b.split("_")[-1].to_i }

      i = 0
      prots_to_annotate.each do |p|

        if @aln_hits.has_key? p

          hit_index = 0

          if @aln_hits[p][:hits].length > 1
            hit_index = choose_best_hit i, prots, ref_cds
          end

          h = @aln_hits[p][:hits][hit_index]
          hit = ref_cds[h]
          annotations[p] = hit
          annotations[p][:pId] = @aln_hits[p][:pId]
          annotations[p][:length] = @aln_hits[p][:length][hit_index]
          i+=1

        else

          annotations[p] = nil

        end

      end

    elsif ! @aln_hits.empty?

      @aln_hits.each_key do |k|
        contig = k.split("_")[0..-3].join("_")
        if contig == contig_to_annotate
          annotations[k] = @aln_hits[k]
        end
      end

    end

    annotations                 # return

  end                           # end of method


  # Choose Best Hit base on neighbor hits
  def choose_best_hit i, prots, ref_cds

    hit_index = 0
    p = prots[i]
    hit_locus_tags = []

    @aln_hits[p][:hits].each do |h|
      hit_locus_tags << ref_cds[h][:locustag].downcase.split("_")[-1].gsub(/[a-z]/,"").to_i
    end

    continue=true
    offset=1

    while continue
      fwd_end = false
      bcw_end = false
      found = false

      if (i+offset) < (prots.length-1)
        fwd_p = prots[i+offset]
        next_prot_hits = @aln_hits[fwd_p][:hits]
        if next_prot_hits.length < 2
          n = ref_cds[next_prot_hits[0]][:locustag].downcase.split("_")[-1].gsub(/[a-z]/,"").to_i
          closest = 10000
          current_ltag_i = 0
          hit_locus_tags.each_with_index do |ltag,ltag_i|
            if (ltag-n).abs < closest
              current_ltag_i = ltag_i
              closest = (ltag-n).abs
            end
          end
          hit_index = current_ltag_i
          found = true
        end
      else
        fwd_end = true
      end

      if (i-offset) >= 0 and !found
        bcw_p = prots[i-offset]
        next_prot_hits = @aln_hits[bcw_p][:hits]
        if next_prot_hits.length < 2
          n = ref_cds[next_prot_hits[0]][:locustag].downcase.split("_")[-1].gsub(/[a-z]/,"").to_i
          closest = 10000
          current_ltag_i = 0
          hit_locus_tags.each_with_index do |ltag,ltag_i|
            if (ltag-n).abs < closest
              current_ltag_i = ltag_i
              closest = (ltag-n).abs
            end
          end
          hit_index = current_ltag_i
          found = true
        end
      else
        bcw_end = true
      end

      offset += 1
      continue = (!fwd_end and !bcw_end and !found)
    end

    hit_index

  end                           # end of method

  def prune_aln_hits aln_hits

    keys_to_delete = []

    aln_hits.each do |key1,val1|

      aln_hits.each do |key2,val2|

        next if key1==key2
        next if keys_to_delete.include? key1
        next if keys_to_delete.include? key2

        if val1[:query_location][0][0] >= val2[:query_location][0][0] and
          val1[:query_location][0][0] < val2[:query_location][0][1]
          overlap_len = val2[:query_location][0][1] - val1[:query_location][0][0]
          val1_len = val1[:query_location][0][1]-val1[:query_location][0][0]
          val2_len = val2[:query_location][0][1]-val2[:query_location][0][0]
          if overlap_len.to_f/val1_len > 0.2 and overlap_len.to_f/val2_len > 0.2
            if val1[:score] < val2[:score]
              keys_to_delete << key1
            else
              keys_to_delete << key2
            end
          end
        elsif val2[:query_location][0][0] >= val1[:query_location][0][0] and
             val2[:query_location][0][0] < val1[:query_location][0][1]
          overlap_len = val1[:query_location][0][1] - val2[:query_location][0][0]
          val1_len = val1[:query_location][0][1]-val1[:query_location][0][0]
          val2_len = val2[:query_location][0][1]-val2[:query_location][0][0]
          if overlap_len.to_f/val1_len > 0.2 and overlap_len.to_f/val2_len > 0.2
            if val1[:score] < val2[:score]
              keys_to_delete << key1
            else
              keys_to_delete << key2
            end
          end
        end

      end

    end

    keys_to_delete.each do |k|
      aln_hits.delete(k)
    end

  end                           # end of method


end                             # end of class
