Gem::Specification.new do |s|
  s.name        = 'bacterial-annotator'
  s.version     = '0.9.1'
  s.date        = '2022-09-05'
  s.summary     = "Bacterial Annotator"
  s.description = "GEM to annotate bacterial genome sequence based on a reference genome and complete the annotation with an external database."
  s.authors     = ["Maxime Deraspe"]
  s.email       = 'maximilien1er@gmail.com'
  s.files       = ["lib/bacterial-comparator.rb",
                   "lib/bacterial-annotator.rb",
                   "lib/bacterial-identificator.rb",
                   "lib/helper.rb",
                   "lib/bacterial-annotator/sequence-annotation.rb",
                   "lib/bacterial-annotator/sequence-fasta.rb",
                   "lib/bacterial-annotator/sequence-synteny.rb"]
  s.homepage    = 'http://rubygems.org/gems/bacterial-annotator'
  s.license       = 'GPL-3.0'
  s.require_path = 'lib'
  s.bindir = "bin"
  s.executables = [
    "bacterial-annotator",
    "ba_prodigal",
    "ba_blat",
    "ba_mafft",
    "ba_raxml",
    "ba_diamond",
    "ba_fasta36",
    "ba_cdhit",
    "ba_fasttree",
    "ba_mash",
    "ba_trnascan"
  ]
  s.add_runtime_dependency 'bio', '~> 1.4', '>= 1.4.3'
  s.add_runtime_dependency 'mechanize', '~>2.7', '>=2.7.3'
  s.add_runtime_dependency 'parallel', '~>1.9', '>=1.9.0'
end
