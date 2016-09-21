# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maxime@deraspe.net
# review:  	
# date:    	15-02-24
# version: 	0.0.1
# licence:  	



class SyntenyManip

  attr_reader :query_file, :subject_file, :aln_hits

  def initialize query_file, subject_file, name, pidentity
    @query_file = query_file
    @subject_file = subject_file
    @name = name
    @pidentity = pidentity
    @aln_file = nil
  end                           # end of initialize

  # run blat on proteins
  def run_blat root, outdir
    system("#{root}/blat.linux -out=blast8 -minIdentity=#{@pidentity} -prot #{@subject_file} #{@query_file} #{outdir}/#{@name}.blat8.tsv")
    @aln_file = "#{outdir}/#{@name}.blat8.tsv"
    # extract_hits
  end                           # end of method

  # Extract Hit from blast8 file and save it in hash
  # contig-0_1      ABJ71957.1      96.92   65      2       0       1       65      1       65      9.2e-31 131.0
  def extract_hits mode

    @aln_hits = {}
    File.open(@aln_file,"r") do |fread|
      while l = fread.gets
        lA = l.chomp!.split("\t")
        key = lA[0]
        if mode == :refgenome
          hit = lA[1]
        elsif mode == :externaldb
          hit = lA[1].chomp.split("|")[1]
        end
        if ! @aln_hits.has_key? key
          next if lA[2].to_f < @pidentity
          @aln_hits[key] = {
            pId: lA[2].to_f.round(2),
            evalue: lA[10],
            score: lA[11].to_f,
            hits: [hit],
            length: [lA[3].to_i],
            query_location: [[lA[6].to_i,lA[7].to_i]],
            subject_location: [[lA[8].to_i,lA[9].to_i]]
          }
        elsif lA[11].to_f > @aln_hits[key][:score]
          @aln_hits[key] = {
            pId: lA[2].to_f.round(2),
            evalue: lA[10],
            score: lA[11].to_f,
            hits: [hit],
            length: [lA[3].to_i],
            query_location: [[lA[6].to_i,lA[7].to_i]],
            subject_location: [[lA[8].to_i,lA[9].to_i]]
          }
        elsif lA[11].to_f == @aln_hits[key][:score]
          @aln_hits[key][:hits] << hit
          @aln_hits[key][:length] << lA[3].to_i
          @aln_hits[key][:query_location] << [lA[6].to_i,lA[7].to_i]
          @aln_hits[key][:subject_location] << [lA[8].to_i,lA[9].to_i]
        end
      end
    end

  end                           # end of method



  # Get the annotations for a contig for RerenceGenome
  def get_annotation_for_contig prots_to_annotate, ref_cds

    return {} if prots_to_annotate == nil

    contig_to_annotate = prots_to_annotate[0].split("_")[0..-2].join("_")
    annotations = {}
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



end                             # end of class
