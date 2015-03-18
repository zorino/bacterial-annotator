# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maxime@deraspe.net
# review:  	
# date:    	15-02-24
# version: 	0.0.1
# licence:  	



class FastaManip

  attr_reader :fasta_flat, :fasta_file, :prodigal_files

  # Initialize fasta holder
  def initialize fasta_file

    @fasta_file = fasta_file
    @fasta_flat = Bio::FlatFile.auto(@fasta_file)
    @prodigal_files = nil
    @single_fasta = nil

    if @fasta_flat.dbclass != Bio::FastaFormat
      abort "Aborting : The input sequence is not a fasta file !"
    end

  end

  # Run prodigal on the genome to annotate
  def run_prodigal root, outdir
    @prodigal_files = {}
    Dir.mkdir "#{outdir}" if ! Dir.exists? "#{outdir}"
    system("#{root}/prodigal.linux -i #{@fasta_file} -a #{outdir}/Proteins.fa -d #{outdir}/Genes.fa -o #{outdir}/Genbanks.gbk")
    @prodigal_files = {multiGBK: "#{outdir}/Genbanks.gbk",
                       contigs: [],
                       contigs_length: [],
                       genes: "#{outdir}/Genes.fa",
                       proteins: "#{outdir}/Proteins.fa",
                       prot_ids_by_contig: {},
                       fasta_path: "#{outdir}/single-fasta/",
                       gbk_path: "#{outdir}/single-genbank/"}
    split_fasta outdir
    split_genbank outdir, "#{outdir}/Genbanks.gbk"
    extract_cds_names
    @prodigal_files
  end


  # Split Multi Genbanks file
  # RETURN : array of fasta files
  def split_fasta outdir
    @single_fasta = {}
    Dir.mkdir("#{outdir}/single-fasta") if ! Dir.exists?("#{outdir}/single-fasta")
    @fasta_flat.each_entry do |seq|
      file_name = seq.definition.chomp.split(" ")[0]
      @prodigal_files[:contigs] << "#{file_name}"
      @prodigal_files[:contigs_length] << seq.seq.length
      File.open("#{outdir}/single-fasta/#{file_name}.fasta", "w") do |fwrite| 
        fwrite.write(seq)
      end
      @single_fasta[file_name] = seq
    end
  end


  # Split Multi Genbanks file
  # RETURN : array of genbank files
  def split_genbank outdir, multigbk

    Dir.mkdir("#{outdir}/single-genbank")if ! Dir.exists?("#{outdir}/single-genbank")
    File.open(multigbk,"r") do |f|
      fopen = nil
      while l = f.gets
        if l[0..9] == "DEFINITION"
          file_name = l.chomp.split(";")[2].gsub("seqhdr","").delete("\"").delete("=").split(" ")[0]
          outseq, seq_length = print_sequence_for_gbk @single_fasta[file_name]
          spacer = " " * (20-seq_length.to_s.length)
          date = DateTime.now
          month = Date::ABBR_MONTHNAMES[date.month]
          day = "%02d" % date.day
          year = date.year
          locus = "LOCUS       #{file_name}#{spacer}#{seq_length.to_s} bp    DNA     linear   BCT #{day}-#{month}-#{year}\n"
          fopen = File.open("#{outdir}/single-genbank/#{file_name}.gbk", "w")
          fopen.write(locus)
          fopen.write(l)
        elsif l[0..1] == "//"
          fopen.write(outseq)
          fopen.close
        else
          fopen.write(l)
        end
      end
    end

  end


  # Utility function to print the sequence to the end of a gbk file
  def print_sequence_for_gbk seq

    outseq = "ORIGIN\n"
    # puts "ORIGIN"

    ntNum = 0
    sequence = seq.seq.downcase

    nt_left = true
    it = 0

    while nt_left

      if sequence.length > it+60
        nt_to_add = sequence[it..(it+59)]
        # printf "%9s ", (ntNum - l.size + 2)
        outseq += "%9s " % (it+1)
        outseq += nt_to_add.scan(/.{1,10}/).join(" ")
        outseq += "\n"
        it += 60
      else
        nt_to_add = sequence[it..sequence.length-1]
        outseq += "%9s " % (it+1)
        outseq += nt_to_add.scan(/.{1,10}/).join(" ")
        outseq += "\n"
        outseq += "//"
        nt_left = false
      end

    end
    
    return outseq, sequence.length

  end


  # extract protein and gene names
  def extract_cds_names

    prot_ids = {}
    flatfile = Bio::FlatFile.auto(@prodigal_files[:proteins])
    flatfile.each_entry do |entry|
      prot_id = entry.definition.split(" ")[0]
      contig = prot_id.split("_")[0..-2].join("_")
      if !prot_ids.has_key? contig
        prot_ids[contig] = []
      end
      prot_ids[contig] << prot_id
    end

    prot_ids.each do |k,prot_array|
      prot_array.sort! { |a,b| a.split("_")[-1].to_i <=> b.split("_")[-1].to_i }
    end

    @prodigal_files[:prot_ids_by_contig] = prot_ids

  end


  private :extract_cds_names # :split_fasta, :split_genbank

end
