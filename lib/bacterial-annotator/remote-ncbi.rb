# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maxime@deraspe.net
# review:  	
# date:    	15-02-24
# version: 	0.0.1
# licence:  	

require 'mechanize'
require 'open-uri'
require 'bio'

class RemoteNCBI

  attr_reader :aln_hits, :db, :xmloutput

  # initialize stuff for a remote ncbi run
  def initialize db, seq_file, outfile

    if ! ["swissprot", "refseq_protein", "nr"].include? db
      @db = "bad database"
    else
      @db = db
    end

    url = 'http://blast.ncbi.nlm.nih.gov/Blast.cgi'\
          '?PROGRAM=blastp&BLAST_PROGRAMS=blastp'\
          '&PAGE_TYPE=BlastSearch&SHOW_DEFAULTS=on'\
          '&LINK_LOC=blasthome'

    @seq_file = seq_file
    @outfile = outfile
    @resultURI = submit_blast url

    @xmloutput = ""
    @valid = validate_output

  end                           # end of method


  # submit blast to ncbi
  def submit_blast ncbiURL

    seq_fasta = File.read(@seq_file)
    a = Mechanize.new { |agent|
      agent.user_agent_alias = 'Linux Firefox'
      agent.ignore_bad_chunking = true
    }

    resultURI = ""

    a.get(ncbiURL) do |page|

      search = page.form_with(:name => 'searchForm') { |form|
        form.textareas[0].value = File.read(@seq_file)
        form.field_with(:name => 'DATABASE').value = @db
        form.field_with(:name => 'MAX_NUM_SEQ').value = 10     
      }.submit

      toBreak = 0
      requestID = ""
      search.parser.css('td').each do |td|
        if toBreak == 1
          requestID = td.text.gsub(" ","")
          break
        end
        if td.text == "Request ID"
          toBreak = 1
        end
      end

      resultURI = URI.parse("http://blast.ncbi.nlm.nih.gov/Blast.cgi?RESULTS_FILE=on&RID=#{requestID}&FORMAT_TYPE=XML&FORMAT_OBJECT=Alignment&CMD=Get")

      puts URI.parse("http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&RID=#{requestID}")

    end

    resultURI

  end                           # end of method


  # validate the xml blast results
  def validate_output

    xmloutput = ""
    valid = true
    finish = false
    while valid and ! finish
      response = Net::HTTP.get_response(@resultURI)
      body = response.body.split("\n")

      if body[0] =~ /<?xml version=/
        xmloutput = body.join("\n")
        valid = true
        finish = true
      else
        valid = false
        body.each do |l|
          if l =~ /Status=/
            status = l.strip.gsub("Status=", "")
            if status == "WAITING"
              valid = true
            end
          end
          break if valid
        end
      end
      sleep 3
    end

    if finish
      File.open("#{@outfile}", "w") do |f|
        f.write(xmloutput)
      end
      return finish
    end
    valid

  end                           # end of method

  # extract blast results from 
  def extract_blast_results

    if !@valid
      @aln_hits = nil
      return
    end

    flat = Bio::FlatFile.auto("#{@outfile}")
    @aln_hits = {}

    flat.each_entry do |report|
      
      report.iterations.each do |query_it|
        prot_id = query_it.query_def.split(" ")[0]
        query_it.hits.each do |hit|
          if ! @aln_hits.has_key? prot_id
            p_identity = hit.identity.to_f/hit.target_len.to_f*100
            if p_identity > 70
              # cleaning product definition
              product = hit.definition.gsub("MULTISPECIES: ","").
                        gsub(/ \[.*\]/,"").
                        gsub("RecName: Full=","").
                        split("; AltName")[0].
                        split("; Flags:")[0].
                        split(" ; Short=")[0]
              gi = hit.hit_id.to_s.split("|")[1]
              organism = ""
              if ! hit.definition[/\[.*\]/].nil?
                organism = hit.definition[/\[.*\]/].gsub("[","").gsub("]","")
              end
              @aln_hits[prot_id] = {
                pId: hit.identity.to_f/hit.target_len.to_f*100,
                length: hit.target_len.to_i,
                evalue: hit.evalue,
                score: hit.bit_score.to_f,
                hits: [{gi: gi, product: product, org: organism}]
              }
            end
          end
        end
      end
    end

  end                           # end of method


end                             # end of class
