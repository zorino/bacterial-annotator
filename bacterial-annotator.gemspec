Gem::Specification.new do |s|
  s.name        = 'bacterial-annotator'
  s.version     = '0.0.1'
  s.date        = '2015-02-24'
  s.summary     = "Bacterial Annotator"
  s.description = "Annotate bacterial genomes from a draft or complete genome based on a reference genome."
  s.authors     = ["Maxime Deraspe"]
  s.email       = 'maxime@deraspe.net'
  s.files       = ["lib/bacterial-annotator.rb",
                   "lib/bacterial-annotator/genbank-manip.rb",
                   "lib/bacterial-annotator/fasta-manip.rb",
                   "lib/bacterial-annotator/synteny-manip.rb"]
  s.homepage    = 'http://rubygems.org/gems/bacterial-annotator'
  s.license       = 'GPLv3'
  s.require_path = 'lib'
  s.bindir = "bin"
  s.executables = [
    "bacterial-annotator",
    "ba_prodigal",
    "ba_blat"
    # "prodigal.linux",
    # "blat.linux",
    # "multifasta-manip",
    # "genbank-manip"
  ]
  s.default_executable = "bacterial-annotator"
  s.add_runtime_dependency 'bio', '~> 1.4', '>= 1.4.3'
end
