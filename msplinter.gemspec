# coding: utf-8
lib = File.expand_path('../lib', __FILE__)
$LOAD_PATH.unshift(lib) unless $LOAD_PATH.include?(lib)
require 'msplinter/version'

Gem::Specification.new do |spec|
  spec.name          = "msplinter"
  spec.version       = Msplinter::VERSION
  spec.authors       = ["John T. Prince"]
  spec.email         = ["jtprince@gmail.com"]
  spec.description   = %q{Predicts how molecules will fragment in a mass spectrometer.  Currently focused on lipid fragmentation under CID, HCD or PQD.}
  spec.summary       = %q{Predicts how molecules will fragment in a mass spectrometer}
  spec.homepage      = "https://github.com/princelab/msplinter"
  spec.license       = "MIT"

  spec.files         = `git ls-files`.split($/)
  spec.executables   = spec.files.grep(%r{^bin/}) { |f| File.basename(f) }
  spec.test_files    = spec.files.grep(%r{^(test|spec|features)/})
  spec.require_paths = ["lib"]

  [
    ["rubabel", "~> 0.4.1"],
  ].each do |args|
    spec.add_dependency(*args)
  end
  [
    ["bundler", "~> 1.3"],
    ["rake"],
    ["rspec", "~> 2.13.0"], 
    ["rdoc", "~> 3.12"], 
    ["simplecov"],
  ].each do |args|
    spec.add_development_dependency(*args)
  end

end
