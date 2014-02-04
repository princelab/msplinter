require 'rubabel/molecule/fragmentable'

module Rubabel
  class Molecule  # OR should this be Fragment
    def to_mechanism_plot(file, mode: :positive, ms_level: 2, rules: RULES, uniq: false, errors: :remove)
      case mode
      when :positive
        self.correct_for_ph!(1.0)
      when :neutral
        self.correct_for_ph!(7.0)
      when :negative
        self.correct_for_ph!(13.0)
      end
      preliminary_fragments = self.fragment(ms_level, rules, uniq, errors)
      preliminary_fragments.flatten.map do |mol|
      
      end
  end
end


  # written svgs, atoms_count, linkages... 

  # needed methods
  #1 extract necessary parts from each SVG
  #2 resize extracted_SVG?
  #3 place SVG in canvas
  #4 add labels/arrows?
