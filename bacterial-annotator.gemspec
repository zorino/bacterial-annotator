Gem::Specification.new do |s|
  s.name        = 'bacterial-annotator'
  s.version     = '0.3.6'
  s.date        = '2017-01-19'
  s.summary     = "Bacterial Annotator"
  s.description = "GEM to annotate bacterial genome sequence based on a reference genome and complete the annotation with an external database or a remote database."
  s.authors     = ["Maxime Deraspe"]
  s.email       = 'maxime@deraspe.net'
  s.files       = ["lib/bacterial-comparator.rb",
                   "lib/bacterial-annotator.rb",
                   "lib/bacterial-annotator/genbank-manip.rb",
                   "lib/bacterial-annotator/fasta-manip.rb",
                   "lib/bacterial-annotator/synteny-manip.rb",
                   "lib/bacterial-annotator/remote-ncbi.rb"]
  s.homepage    = 'http://rubygems.org/gems/bacterial-annotator'
  s.license       = 'GPL-3.0'
  s.require_path = 'lib'
  s.bindir = "bin"
  s.executables = [
    "bacterial-annotator",
    "ba_prodigal",
    "ba_blat",
    "ba_mafft",
    "ba_raxml"
  ]
  s.default_executable = "bacterial-annotator"
  s.add_runtime_dependency 'bio', '~> 1.4', '>= 1.4.3'
  s.add_runtime_dependency 'mechanize', '~>2.7', '>=2.7.3'
  s.add_runtime_dependency 'parallel', '~>1.9', '>=1.9.0'
end
