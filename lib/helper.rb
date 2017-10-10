# -*- coding: utf-8 -*-
# author:  	maxime déraspe
# email:	maximilien1er@gmail.com
# review:  	
# date:    	2017-10-05
# version: 	0.0.1
# licence:  	



class Helper

  def self.sec2str secs
    [[60, :seconds], [60, :minutes], [24, :hours], [1000, :days]].map{ |count, name|
      if secs > 0
        secs, n = secs.divmod(count)
        "#{n.to_i} #{name}"
      end
    }.compact.reverse.join(' ')
  end

end
