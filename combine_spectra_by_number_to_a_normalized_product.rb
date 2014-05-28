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

  opts.on("-lMANDATORY", "--list=MANDATORY x,y,z", Array, "List of scan numbers to extract") do |scans|
    options[:selected_scans] = scans.map(&:to_i)
  end


end

opt_parse.parse!(ARGV)

if ARGV.size != 1
  puts opt_parse
  exit
end

file = ARGV.first

peaklists = []

# pull the intensity list for overall normalization prior
Mspire::Mzml.open(file) do |mzml|
  options[:selected_scans].map do |i|
    peaklists << mzml[i-1].to_peaklist
  end
end
product_peaklist = Mspire::Peaklist.merge(peaklists)# add non-default options here

#
## WRITE A SPECTRUM for MZML
#spec1 = Mspire::Mzml::Spectrum.new('scan=1') do |spec|
#  # profile and ms_level 1
#  spec.describe_many!(['MS:1000128', ['MS:1000511', 1]])
#  spec.data_arrays = [
#    Mspire::Mzml::DataArray[1,2,3].describe!('MS:1000514'),  
#    Mspire::Mzml::DataArray[4,5,6].describe!('MS:1000515')   
#  ]
#  spec.scan_list = Mspire::Mzml::ScanList.new do |sl|
#    scan = Mspire::Mzml::Scan.new do |scan|
#      # retention time of 42 seconds
#      scan.describe! 'MS:1000016', 40.0, 'UO:0000010'
#    end
#    sl << scan
#  end
#end


# Write a MZML
mzml = Mspire::Mzml.new do |mzml|
  mzml.id = 'ms1_and_ms2'
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
    spectrum_list = Mspire::Mzml::SpectrumList.new(default_data_processing, Mspire::Mzml::Spectrum.from_arrays("combinedspectrum", product_peaklist.transpose) )
    run.spectrum_list = spectrum_list
  end
end

mzml.write("writtenxml.mzML")
