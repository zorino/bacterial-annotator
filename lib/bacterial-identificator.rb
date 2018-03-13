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
require 'json'
require 'pp'


class BacterialIdentificator

  attr_reader :genomes_list, :stats

  # Initialize BacterialIdentificator
  # options[:input], options[:refgenome], ROOT, options[:outdir], options)
  def initialize options, root

    @root = root
    @mash_file = options[:mash_file]
    @genome_list = options[:genome_list]
    @proc = options[:proc].to_i
    @output=options[:output]

  end


  def run_identification

    @genome_hits = {}
    @genome_list.each do |g|
      @genome_hits[g] = []
    end

    Parallel.map(@genome_list, in_threads: @proc) do |g|
      @genome_hits[g] = mash_genome g
    end

    print_output

  end


  def mash_genome genome

    # Reference-ID, Query-ID, Mash-distance, P-value, and Matching-hashes
    # fields = ["hit","query","distance","pvalue","match"]

    results_raw = `#{@root}/mash.linux dist #{@mash_file} #{genome}`
    results = []

    results_raw.split("\n").each do |l|
      lA = l.chomp.split("\t")
      next if lA[-1].split("/")[0] == '0' # no match
      results << (lA[0..0] + lA[2..-1])
    end

    results_sorted = results.sort {|a,b| a[1] <=> b[1]}

    return results_sorted

  end

  # consensus species model
  def consensus_reference

    all_hits = {}
    @genome_hits.each do |g, hits|
      hits.each do |h|
        score = h[3].split("/")[0].to_i
        if ! all_hits.has_key? h[0]
          all_hits[h[0]] = score
        else
          all_hits[h[0]] += score
        end
      end
    end
    return all_hits.sort_by { |k,v| v }.to_h

  end

  # print json
  def print_output

    case @output.downcase
    when "csv"
      @genome_hits.each do |g, hits|
        hits.each do |h|
          puts "#{g},#{h.join(',')}"
        end
      end
    when "json"
      new_genome_hits = {}
      @genome_hits.each do |g, hits|
        new_genome_hits[g] = []
        hits.each do |h|
          new_genome_hits[g].push(Hash[["hit","distance","e-value","score"].zip(h)])
        end
      end
      puts JSON.pretty_generate({genomes: new_genome_hits, summary: summary})
    else
      @genome_hits.each do |g, hits|
        hits.each do |h|
          out = h.join("\t")
          puts "#{g}\t#{out}"
        end
      end
    end

  end

  def summary

    genome_hit_association =  {}

    @genome_hits.each do |g, hits|
      genome_hit_association[hits[0][0]] = 0 if ! genome_hit_association.has_key? hits[0][0]
      genome_hit_association[hits[0][0]] += 1
    end

    population = {
      consensus: consensus_reference.first[0],
      genome_hits: genome_hit_association
    }

    return population

  end


end
