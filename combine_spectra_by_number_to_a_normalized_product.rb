#!/usr/bin/env ruby

require 'mspire/mzml'
require 'mspire/peaklist'
require 'mspire/spectrum'
require 'optparse'
require 'pry'

options = {}
opt_parse = OptionParser.new do |opts|
  opts.on_tail("-h", "--help", "Show this message") do 
    puts opts
    exit
  end
  opts.banner = "Usage: #{__FILE__} file.mzML -l i,j,k "
  opts.separator "Takes a spectral file, a list of scan numbers, and kicks them out as a single normalized spectrum"
  opts.separator "Output: file_extracted.mzML"

  opts.on("-lMANDATORY", "--list=MANDATORY x,y,z", Array, "List of scan numbers to extract", "(assumes that scan number is +1 spectral index)") do |scans|
    options[:selected_scans] = scans.map(&:to_i)
  end


end

opt_parse.parse!(ARGV)

if ARGV.size != 1
  puts opt_parse
  exit
end

file = ARGV.first
output_file = file.sub(File.extname(file), '_extracted.mzML')

# pull the intensity list for overall normalization prior to iterating through these
peaklists = Mspire::Mzml.open(file) do |mzml|
  options[:selected_scans].map do |i|
    sp = mzml[i-1]
    normalizer = 100.0 / sp.intensities.max
    sp.intensities.map! {|i| i * normalizer }
    sp.to_peaklist
  end
end
product_peaklist = Mspire::Peaklist.merge(peaklists)# add non-default options here


# Write a MZML
mzml = Mspire::Mzml.new do |mzml|
  mzml.id = 'merged spectra'
  mzml.cvs = Mspire::Mzml::CV::DEFAULT_CVS
  mzml.file_description = Mspire::Mzml::FileDescription.new  do |fd|
    fd.file_content = Mspire::Mzml::FileContent.new
    fd.source_files << Mspire::Mzml::SourceFile.new
  end
  default_instrument_config = Mspire::Mzml::InstrumentConfiguration.new("IC").describe!('MS:1000031')
  mzml.instrument_configurations << default_instrument_config
  software = Mspire::Mzml::Software.new
  mzml.software_list << software
  normalize_processing = Mspire::Mzml::DataProcessing.new("ms1_normalization") do |dp|
    # 'MS:1001484' -> intensity normalization 
    dp.processing_methods << Mspire::Mzml::ProcessingMethod.new(software).describe!('MS:1001484')
  end
  default_data_processing = Mspire::Mzml::DataProcessing.new("combined spectra")
  mzml.data_processing_list << normalize_processing
  mzml.data_processing_list << default_data_processing
  mzml.run = Mspire::Mzml::Run.new("little_run", default_instrument_config) do |run|
    spectrum_list = Mspire::Mzml::SpectrumList.new(default_data_processing, [Mspire::Mzml::Spectrum.from_arrays("combinedspectrum", product_peaklist.transpose)] )
    run.spectrum_list = spectrum_list
  end
  mzml
end

puts "Writing output file: #{output_file}"
mzml.write(output_file)
