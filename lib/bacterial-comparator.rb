# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maxime@deraspe.net
# review:  	
# date:    	15-02-24
# version: 	0.0.1
# licence:  	

require 'bio'
require 'fileutils'

class BacterialComparator

  attr_reader :genomes_list, :stats

  # Initialize BacterialAnnotator
  # options[:input], options[:refgenome], ROOT, options[:outdir], options)
  def initialize genomes_list, root

    @outdir = "phylogenomics"
    Dir.mkdir(@outdir) if ! Dir.exists? @outdir
    @genomes_list = genomes_list
    @ref_prot = get_ref_prot
    @synteny = read_prot_synteny

    # TODO ARGS for min coverage,  min pId.
    min_cov = 0.8
    min_pid = 0.9

    @stats = extract_syntenic_fasta min_cov, min_pid

  end

  def read_prot_synteny
    synteny = {}
    @genomes_list.each do |g|
      puts "#{g}/Prot-Synteny.tsv"
      file = File.open("#{g}/Prot-Synteny.tsv", "r")
      l = file.gets             # skip header
      while l = file.gets
        # AAK98805.1      spr0001 453     1.0     100.0   ABAC01000005_14 453     1.0
        lA = l.chomp.split("\t")
        synteny[lA[0]] =  [] if ! synteny.has_key? lA[0]
        synteny[lA[0]] << {ref_cov: lA[3].to_f, pId: lA[4].to_f, query_prot: lA[5], query_cov: lA[7].to_f}
      end
      file.close
    end
    synteny
  end

  def get_ref_prot
    ref_prot = []
    pep_file = Dir["#{@genomes_list[0]}/*.pep"]
    flatfile = Bio::FlatFile.auto("#{pep_file[0]}")
    flatfile.each_entry do |entry|
      ref_prot << entry.definition.split(" ")[0]
    end
    ref_prot
  end


  def get_sequence_from_flatfile flatfile, name

    out = ""
    flatfile.each_entry do |entry|
      if entry.definition.split(" ")[0] == name
        bioseq = Bio::Sequence.auto(entry.seq)
        out = bioseq.output_fasta("#{name}",60)
      end
    end
    out

  end


  def build_multifasta ref_prot, synteny

    pep_out_dir = "./#{@outdir}/genes-align-pep"
    Dir.mkdir(pep_out_dir) if ! Dir.exists? pep_out_dir
    dna_out_dir = "./#{@outdir}/genes-align-dna"
    Dir.mkdir(dna_out_dir) if ! Dir.exists? dna_out_dir

    pep_out = File.open(pep_out_dir+"/#{ref_prot}.pep", "w")
    dna_out = File.open(dna_out_dir+"/#{ref_prot}.dna", "w")

    # create multifasta by syntenic proteins
    pep_file = Dir["#{@genomes_list[0]}/*.pep"]
    flatfile = Bio::FlatFile.auto("#{pep_file[0]}")
    pep_out.write(get_sequence_from_flatfile flatfile, ref_prot)
    @genomes_list.each_with_index do |g,i|
      flatfile = Bio::FlatFile.auto("#{g}/Proteins.fa")
      pep_out.write(get_sequence_from_flatfile flatfile, synteny[i][:query_prot])
    end

    # create multifasta by syntenic proteins
    dna_file = Dir["#{@genomes_list[0]}/*.dna"]
    flatfile = Bio::FlatFile.auto("#{dna_file[0]}")
    dna_out.write(get_sequence_from_flatfile flatfile, ref_prot)
    @genomes_list.each_with_index do |g,i|
      flatfile = Bio::FlatFile.auto("#{g}/Genes.fa")
      dna_out.write(get_sequence_from_flatfile flatfile, synteny[i][:query_prot])
    end

  end


  def extract_syntenic_fasta min_cov, min_pid

    "# Extracting Proteins and Genes multifasta.."
    nb_of_syntenic = 0
    stats = {}
    stats[:syntenic] = []
    fout = File.open("#{@outdir}/cds-synteny.tsv", "w")
    fout.write("Gene\t"+@genomes_list.join("\t")+"\n")

    @synteny.each do |k,v|
      is_syntenic = 1
      v.each do |v_|
        if v_[:query_cov].nil?
          is_syntenic = 0
          break
        elsif v_[:query_cov] > min_cov and
             v_[:ref_cov] > min_cov and
             v_[:pId] > min_pid
        # synteny -> great !
        else
          is_syntenic = 0
          break
        end
      end

      if is_syntenic == 1
        nb_of_syntenic += 1
        build_multifasta k, v
        fout.write("#{k}")
        v.each do |x|
          fout.write("\t#{x[:query_prot]}|#{x[:query_cov]}|#{x[:ref_cov]}")
          stats[:syntenic] << k
        end
        fout.write("\n")
      end
    end

    fout.close
    stats[:nb_of_syntenic] = nb_of_syntenic
    puts "Syntenic genes : " + nb_of_syntenic.to_s + " / " + @ref_prot.length.to_s

  end


  def mafft_align_all_pep root
    puts "# MAFFT multialign all protein sequences.."
    Dir.chdir("#{@outdir}/genes-align-pep/")
    Dir["*.pep"].each do |f|
      system("#{root}/mafft.linux #{f} > #{f}.aln")
    end

  end

  def mafft_align_all_dna root
    puts "# MAFFT multialign all gene sequences.."


  end


end                             # end of Class
