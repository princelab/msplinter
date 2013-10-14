require 'rubabel/molecule/fragmentable'

module Rubabel
  class Molecule
    # define output of a molecule to MGF file format.  
    # Input: filename for the mgf file to be written
    # options: 
      # ms_level: the level of fragmentation desired
      # mode: :positive, :neutral, or :negative, signifying the pH at which the fragments will exist
      # rules: the fragmentation rules desired
      # uniq: if you only want unique fragments
      # errors: check normal documentation for this usage
      # standard_intensity_value: the default value for the intensity value for each fragment
    def to_mgf(file, mode: :positive, ms_level: 2, rules: RULES, uniq: false, errors: :remove, standard_intensity_value: 2000)
      case mode
      when :positive
        self.correct_for_ph!(1.0)
      when :neutral
        self.correct_for_ph!(7.0)
      when :negative
        self.correct_for_ph!(13.0)
      end
      preliminary_fragments = self.fragment(ms_level, rules, uniq, errors).flatten
      fragments = preliminary_fragments.select {|mol| mol.charge > 0 }
      fragments = fragments.map(&:mass).uniq.sort
      File.open(file, 'w') do |io|
        io.puts "BEGIN IONS"
        io.puts "TITLE=#{File.basename(file).sub(".mgf",'')}"
        io.puts "PEPMASS=#{self.mass}"
        fragments.zip(Array.new(fragments.size, standard_intensity_value)) {|r| io.puts r.join(" ")}
        io.puts "END IONS"
      end
    end

  end
end
