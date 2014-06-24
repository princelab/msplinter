require 'mspire/mzml'
m = Rubabel["CCCCCCCCCCCCC/C=C/[C@@H]([C@@H](CO)NC(=O)CCCCC)O"]
m.adducts << "[Li+]"
mols = m.fragment
masses = mols.flatten.map(&:mass_with_adduct).flatten.sort
intensities = Array.new(masses.size) { 1000 }
spectra = Mspire::Mzml::Spectrum.from_arrays(0, [masses, intensities] )
outfile = "C6.mzml"
puts "writing outfile: #{outfile}"
require 'mspire/mzml'
file = Mspire::Mzml.new do |mzml|
  mzml.id = 'ms1_and_ms2'
  mzml.cvs = Mspire::Mzml::CV::DEFAULT_CVS
  mzml.file_description = Mspire::Mzml::FileDescription.new  do |fd|
    fd.file_content = Mspire::Mzml::FileContent.new
    fd.source_files << Mspire::Mzml::SourceFile.new
  end
  default_instrument_config = Mspire::Mzml::InstrumentConfiguration.new("IC").describe!('MS:1000031')
  mzml.instrument_configurations << default_instrument_config
  mspire = Mspire::Mzml::Software.new
  mzml.software_list.push(mspire).uniq_by(&:id)
  spectra_generation = Mspire::Mzml::DataProcessing.new("msplinter fragmentation prediction") do |dp|
    # "MS:1001458" -> spectrum_generation_information
    dp.processing_methods << Mspire::Mzml::ProcessingMethod.new(mspire).describe!('MS:1001458')
  end

  mzml.data_processing_list << spectra_generation

  mzml.run = Mspire::Mzml::Run.new("artificial_run", default_instrument_config) do |run|
    spectrum_list = Mspire::Mzml::SpectrumList.new(spectra_generation, [spectra])
    run.spectrum_list = spectrum_list
  end
end
file.write(outfile)

