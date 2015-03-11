# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maxime@deraspe.net
# review:  	
# date:    	15-02-24
# version: 	0.0.1
# licence:  	

require 'bio'

require 'bacterial-annotator/genbank-manip'
require 'bacterial-annotator/fasta-manip'
require 'bacterial-annotator/synteny-manip'


class BacterialAnnotator

  # Initialize BacterialAnnotator
  def initialize fasta, refgenome, root, outdir
    @root = root
    @outdir = outdir
    if File.exists? (@outdir)
      abort "Output directory already exist ! Choose another one"
    end
    Dir.mkdir(@outdir)
    @fasta = FastaManip.new(fasta)
    @refgenome = GenbankManip.new(refgenome, outdir)
    @prot_synteny = nil
    @annotation_stats = {by_contigs: {},
                         annotated_cds: 0,
                         total_cds: 0,
                         foreign_contigs: [],
                         synteny_contigs: [] }
    puts "Successfully loaded #{@refgenome.gbk.definition} genomes"
  end

  # Prepare files for the annotation
  # Will run prodigal on the query and prepare reference genome files
  def prepare_files_for_annotation
    puts "\nRunning Prodigal on your genome.."
    @fasta.run_prodigal @root, @outdir
    @refgenome.write_cds_to_file @outdir
    puts "\nProdigal done."
  end

  # run_alignment of reference genome proteins and the query
  def run_annotation

    @prot_synteny = SyntenyManip.new(@fasta.prodigal_files[:proteins], @refgenome.cds_file)
    puts "\nRunning BLAT alignment.."
    @prot_synteny.run_blat @root, @outdir
    @prot_synteny.extract_hits
    puts "\nParsing annotation to genbank files.."

    @fasta.prodigal_files[:contigs].each do |contig|
      contig_prots = @fasta.prodigal_files[:prot_ids_by_contig][contig]
      contig_prot_annotations = @prot_synteny.get_annotation_for_contig contig_prots, @refgenome.coding_seq
      cumulate_annotation_stats contig, contig_prot_annotations
      gbk_path = @fasta.prodigal_files[:gbk_path]
      gbk_to_annotate = GenbankManip.new("#{gbk_path}/#{contig}.gbk", "#{gbk_path}")
      gbk_to_annotate.add_annotation contig_prot_annotations, gbk_path, 0
    end

    puts "\nPrinting Statistics.."
    print_stats "#{@outdir}/Annotation-Stats.txt"

  end


  # cumulate the stats for the synteny
  def cumulate_annotation_stats contig, contig_prots_ann

    contig_prots = @fasta.prodigal_files[:prot_ids_by_contig][contig]

    @annotation_stats[:total_cds] += contig_prots.length if contig_prots
    contig_prots_ann.each_value do |v|
      @annotation_stats[:annotated_cds] += 1 if v != nil
    end

    # Annotated Contigs
    if contig_prots_ann.keys.length < 1
      @annotation_stats[:foreign_contigs] << contig
    else
      @annotation_stats[:synteny_contigs] << contig
    end

  end


  # print statistics to file
  def print_stats file

    total_nb_contigs = @annotation_stats[:foreign_contigs].length + @annotation_stats[:synteny_contigs].length
    p_contigs_annotated = @annotation_stats[:synteny_contigs].length.to_f/total_nb_contigs.to_f
    p_cds_annotated = @annotation_stats[:annotated_cds].to_f/@annotation_stats[:total_cds].to_f

    File.open(file, "w") do |fopen| 
      fopen.write("#Contigs annotation based on reference genomes\n")
      fopen.write("Foreign Contigs :\t" + @annotation_stats[:foreign_contigs].length.to_s + "\n")
      fopen.write("Annotated Contigs :\t" + @annotation_stats[:synteny_contigs].length.to_s + "\n")
      fopen.write("Total Contigs :\t" + total_nb_contigs.to_s + "\n") 
      fopen.write("% Contigs annotated :\t" + (p_contigs_annotated*100).to_s + "\n")
      fopen.write("\n")

      fopen.write("#CDS annotations based on reference genomes\n")
      fopen.write("Annotated CDS :\t" + @annotation_stats[:annotated_cds].to_s + "\n")
      fopen.write("Total CDS :\t" + @annotation_stats[:total_cds].to_s + "\n")
      fopen.write("% CDS annotated :\t" + (p_cds_annotated*100).to_s + "\n")
      fopen.write("\n")


    end

  end


end
