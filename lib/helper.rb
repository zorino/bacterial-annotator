# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maximilien1er@gmail.com
# review:  	
# date:    	2017-10-05
# version: 	0.0.1
# licence:  	


def run_ncbi_fetch db, id, type, gbk_file
  Bio::NCBI.default_email = 'default@default.com'
  ncbi = Bio::NCBI::REST.new
  genbank = ncbi.efetch(id, {"db"=>db, "rettype"=>type})
  f = File.new(gbk_file, "w")
  f.print(genbank)
  f.close
  return genbank
end


class Helper

  def self.sec2str secs
    [[60, :seconds], [60, :minutes], [24, :hours], [1000, :days]].map{ |count, name|
      if secs > 0
        secs, n = secs.divmod(count)
        "#{n.to_i} #{name}"
      end
    }.compact.reverse.join(' ')
  end

  def self.download_genbank id, gbk_file

    db = "nucleotide"
    type = "gb"

    puts "Downloading genbank #{id}.."
    genbank = run_ncbi_fetch db, id, type, gbk_file

    if type == "gb" and (genbank =~ /^WGS/)
      Dir.mkdir(id) if ! Dir.exists? id
      genbank.split("\n").each do |l|
        if l[0..2] == "WGS"
          puts "# This is a WGS ! ..fetching sub genbank !"
          ids = l.split(/\s+/)[1].split("-")
          if ids.length == 1
            run_ncbi_fetch db, ids[0], type, "#{id}/#{ids[0]}.gbk"
          else
            for i in ids[0][4..-1]..ids[1][4..-1]
              tmp_id = ids[0][0..3]+i.to_s
              puts "  ..fetching #{tmp_id}"
              run_ncbi_fetch db, tmp_id, type, "#{id}/#{tmp_id}.gbk"
              sleep 3
            end
          end
        end
      end
    end

  end

end
