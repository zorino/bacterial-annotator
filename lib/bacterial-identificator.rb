# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maximilien1er@gmail.com
# review:  	
# date:    	17-10-12
# version: 	0.0.1
# licence:  	

require 'bio'
require 'fileutils'
require 'parallel'
require 'helper'

class BacterialIdentificator

  attr_reader :genomes_list, :stats

  # Initialize BacterialIdentificator
  # options[:input], options[:refgenome], ROOT, options[:outdir], options)
  def initialize options, root

    @root = root
    @mash_file = options[:mash_file]
    @genome_list = options[:genome_list]
    @proc = options[:proc].to_i

    @genome_hits = {}
    @genome_list.each do |g|
      @genome_hits[g] = []
    end

    Parallel.map(@genome_list, in_threads: @proc) do |g|
      @genome_hits[g] = mash_genome g
    end

    order_hits = consensus_reference # @genome_hits

    order_hits.each do |k,v|
      puts "#{k} #{v}"
    end

    reference = order_hits[0]

    find_genome_outliers reference

  end


  def mash_genome genome

    # Reference-ID, Query-ID, Mash-distance, P-value, and Matching-hashes
    # fields = ["hit","query","distance","pvalue","match"]

    results_raw = `#{@root}/mash.linux dist #{@mash_file} #{genome}`
    results = []

    results_raw.split("\n").each do |l|
      lA = l.chomp.split("\t")
      next if lA[-1].split("/")[0] == '0' # no match
      results << lA
    end

    results_sorted = results.sort {|a,b| a[2] <=> b[2]}

    File.open("#{genome}.msh_dist", "w") do |fout|
      results_sorted.each do |f|
        fout.write(f.join("\t"))
        fout.write("\n")
      end
    end

    return results_sorted

  end

  # consensus species model
  def consensus_reference

    all_hits = {}
    @genome_hits.each do |g, hits|
      puts "#{g}"
      hits.each do |h|
        score = h[4].split("/")[0].to_i
        if ! all_hits.has_key? h[0]
          all_hits[h[0]] = score
        else
          all_hits[h[0]] += score
        end
      end
    end
    return all_hits.sort_by { |k,v| v }.to_h

  end

  # find genome that are not part of reference species
  def find_genome_outliers reference

    # @genome_hits.each do |g, hits|
    #   # need some exclusion threshold logic
    # end

  end

end

