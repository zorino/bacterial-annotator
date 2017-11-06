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
    @db_path = options[:database]
    @genomes_list = options[:genomes_list]
    @proc = options[:proc].to_i
    p @genomes_list

    @genomes_list.each do |g|
      mash_genome g
    end

  end


  def mash_genome genome

    # Reference-ID, Query-ID, Mash-distance, P-value, and Matching-hashes
    # fields = ["hit","query","distance","pvalue","match"]

    results_raw = `#{@root}/mash.linux dist #{@db_path}/species-sequences.msh #{genome}`
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

  end




end

