# encoding: utf-8
require 'mechanize'
require 'open-uri'
require 'bio'

Blastn = "http://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&BLAST_PROGRAMS=blastn&PAGE_TYPE=BlastSearch&SHOW_DEFAULTS=on&LINK_LOC=blasthome"
Blastp = "http://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastp&BLAST_PROGRAMS=blastp&PAGE_TYPE=BlastSearch&SHOW_DEFAULTS=on&LINK_LOC=blasthome"

def createBlastPage seq, type, db

  blastURL = "http://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastp&BLAST_PROGRAMS=blastp&PAGE_TYPE=BlastSearch&SHOW_DEFAULTS=on&LINK_LOC=blasthome"
  
  a = Mechanize.new { |agent|
    agent.user_agent_alias = 'Linux Firefox'
    agent.ignore_bad_chunking = true
  }

  resultURI = ""

  a.get(blastURL) do |page|

    search = page.form_with(:name => 'searchForm') { |form|
      form.textareas[0].value = seq
      form.field_with(:name => 'DATABASE').value = db
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

    # resultURL = "http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&RID=" + "#{requestID}"
    resultURI = URI.parse("http://blast.ncbi.nlm.nih.gov/Blast.cgi?RESULTS_FILE=on&RID=#{requestID}&FORMAT_TYPE=XML&FORMAT_OBJECT=Alignment&CMD=Get")

    sleep 1

  end

  resultURI

end



def launch_remote_blast db

  seq = "
>contig-0_20 # 16195 # 16854 # -1 # ;gc_cont=0.517
MRILLIEDDMLIGDGIKTGLSKMGFSVDWFTQGRQGKEALYSAPYDAVILDLTLPGMDGR
DILREWREKGQREPVLILTARDALAERVEGLRLGADDYLCKPFALIEVAARLEALMRRTN
GQASNELRHGNVMLDPGKRIATLAGEPLTLKPKEFALLELLMRNAGRVLSRKLIEEKLYT
WDEEVTSNAVEVHVHHLRRKLGSDFIRTVHGIGYTLGEK*
>contig-0_21 # 17006 # 17398 # 1 # ;gc_cont=0.519
MKKFAAVIAVMALCSAPVMAAEQGGFSGPSATQSQAGGFQGPNGSVTTVESAKSLRDDTW
VTLRGNIVERISDDLYVFKDASGTINVDIDHKRWNGVTVTPKDTVEIQGEVDKDWNSVEI
DVKQIRKVNP*
>contig-0_22 # 17451 # 17933 # 1 # ;gc_cont=0.532
MTNLTLDVNIIDFPSIPVAMLPHRCSPELLNYSVAKFIMWRKETGLSPVNQSQTFGVAWD
DPATTAPEAFRFDICGSVSEPIPDNRYGVSNGELTGGRYAVARHVGELDDISHTVWGIIR
HWLPASGEKMRKAPILFHYTNLAEGVTEQRLETDVYVPLA*
>contig-0_23 # 18138 # 18434 # 1 # ;gc_cont=0.401
MEKRTPHTRLSQVKKLVNAGQVRTTRSALLNADELGLDFDGMCNVIIGLSESDFYKSMTT
YSDHTIWQDVYRPRLVTGQVYLKITVIHDVLIVSFKEK*
>contig-0_24 # 18436 # 18831 # 1 # ;gc_cont=0.422
MKCPVCHQGEMVSGIKDIPYTFRGRKTVLKGIHGLYCVHCEESIMNKEESDAFMAQVKAF
RASVNAETVAPEFIVKVRKKLSLTQKEASEIFGGGVNAFSRYEKGNAQPHPSTIKLLRVL
DKHPELLNEIR*
>contig-0_25 # 18964 # 20571 # 1 # ;gc_cont=0.523
MYTRNLLWLVSLVSAAPLYAADVPANTPLAPQQVFRYNNHSDPGTLDPQKVEENTAAQIV
LDLFEGLVWMDGEGQVQPAQAERWEILDGGKRYIFHLRSGLQWSDGQPLTAEDFVLGWQR
AVDPKTASPFAGYLAQAHINNAAAIVAGKADVTSLGVKATDDRTLEVTLEQPVPWFTTML
AWPTLFPVPHHVIAKHGDSWSKPENMVYNGAFVLDQWVVNEKITARKNPKYRDAQHTVLQ
QVEYLALDNSVTGYNRYRAGEVDLTWVPAQQIPAIEKSLPGELRIIPRLNSEYYNFNLEK
PPFNDVRVRRALYLTVDRQLIAQKVLGLRTPATTLTPPEVKGFSATTFDELQKPMSERVA
MAKALLKQAGYDASHPLRFELFYNKYDLHEKTAIALSSEWKKWLGAQVTLRTMEWKTYLD
ARRAGDFMLSRQSWDATYNDASSFLNTLKSDSEENVGHWKNAQYDALLNQATQITDATKR
NALYQQAEVIINQQAPLIPIYYQPLIKLLKPYVGGFPLHNPQDYVYSKELYIKAH*
"

  url = createBlastPage seq, "blastp", db

  # response = Net::HTTP.get_response(url)
  # puts response.methods
  valid = true

  xmloutput = ""

  while valid
    response = Net::HTTP.get_response(url)
    body = response.body.split("\n")

    if body[0] =~ /<?xml version=/
      xmloutput = body.join("\n")
      valid=false
    else
      valid = false
      body.each do |l|
        if l =~ /Status=/
          status = l.strip.gsub("Status=", "")
          puts status
          if status == "WAITING"
            valid = true
            puts url
          end
        end
        break if valid
      end
    end

    sleep 1

  end

  File.open("xmloutput.xml", "w") do |f|
    f.write(xmloutput)
  end
  
end



## LAUNCH BLAST
# xmloutput = launch_remote_blast 'refseq_protein'

flat = Bio::FlatFile.auto("xmloutput.nr.xml")
if ARGV[0]
  flat = Bio::FlatFile.auto(ARGV[0])
end


@aln_hits = {}

flat.each_entry do |report| 

  report.iterations.each do |query_it|
    prot_id = query_it.query_def.split(" ")[0]
    query_it.hits.each do |hit|
      if ! @aln_hits.has_key? prot_id
        product = hit.definition.gsub("MULTISPECIES: ","").
                  gsub(/ \[.*\]/,"").gsub("RecName: Full=","").
                  split("; AltName")[0].split("; Flags:")[0]
        @aln_hits[prot_id] = {
          pId: hit.identity.to_f/hit.target_len.to_f*100,
          length: hit.target_len.to_i,
          evalue: hit.evalue,
          score: hit.bit_score.to_f,
          hits: [product]
        }
      end

    end
    
  end
  
end
# report = Bio::Blast::Report.new(xmloutput)


p @aln_hits
