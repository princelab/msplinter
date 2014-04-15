# Add rules to these two variables as needed
# Rule_names << :rule_names
# Rearrangements << :any_rule_names_that_are_just_rearrangements
require_relative "../fragmentable"

#        self.write( 'root.svg', :add_atom_index => true)
module Rubabel
  class Molecule
    module Fragmentable
      ::Rule_names << def jasms_2002_2_f1(only_uniqs=true)
        fragment_sets = []
        fragment = lambda do |*arr|
          # duplicate
          (nmol, (nitrogen, carbonyl, alpha_carbon, alcohol_oxygen)) = self.dup_molecule(arr)
          # fragmentation chemistry
          alcohol_oxygen.get_bond(alpha_carbon).bond_order = 2
          nmol.delete_bond(carbonyl, alpha_carbon)
          nmol.split
        end

        # Call the block on any search strings which make sense in the context of the code
        self.matches("NC(=O)C(O)C" , only_uniqs).each do |nitrogen, carbonyl, carbonyl_oxygen, alpha_carbon, alcohol_oxygen, carbon|
          fragment_sets << fragment.call(nitrogen, carbonyl, alpha_carbon, alcohol_oxygen)
        end
        fragment_sets.flatten
      end
      ::Rule_names << def jasms_2002_2_f1_aprime(only_uniqs=true)
        fragment_sets = []
        fragment = lambda do |*arr|
          # duplicate
          (nmol, (cyclized_oxygen, linked_to_nitrogen_carbon, attacked_carbon, lost_oxygen)) = self.dup_molecule(arr)
          # fragmentation chemistry
          nmol.add_bond!(cyclized_oxygen, attacked_carbon)
          nmol.delete_bond(lost_oxygen, attacked_carbon)
          nmol.split
        end
        self.matches("OCC(N)CO" , only_uniqs).each do |cyclized_oxygen, carbon, linked_to_nitrogen_carbon, nitrogen,attacked_carbon, lost_oxygen|
          fragment_sets << fragment.call(cyclized_oxygen, linked_to_nitrogen_carbon, attacked_carbon, lost_oxygen)
        end
        fragment_sets.flatten
      end
      ::Rule_names << def jasms_2002_2_f1_bprime(only_uniqs=true)
        fragment_sets = []
        fragment = lambda do |*arr|
          # duplicate
          (nmol, (lost_oxygen, attacked_carbon, nitrogen)) = self.dup_molecule(arr)
          # fragmentation chemistry
          nmol.delete_bond(lost_oxygen, attacked_carbon) #TODO ORDER MATTERS? WTHeck?
          nmol.add_bond!(nitrogen, attacked_carbon)
          nmol.split
        end
        self.matches("O[CH2]C([NH1]C=O)CO" , only_uniqs).each do |lost_oxygen, attacked_carbon, linked_to_nitrogen_carbon, nitrogen, *rest|
          fragment_sets << fragment.call(lost_oxygen, attacked_carbon, nitrogen)
        end
        fragment_sets.flatten
      end
    end # Fragmentable
  end # Molecule
end #Rubabel


if $0 == __FILE__
  require 'bundler/setup'
  require 'rubabel'
  require 'pry'
  mol = Rubabel["LMGP04010962", :lmid]
  mol = Rubabel["CCCCCCCCCCCCCCCC(=O)OC(COP(=O)([O-])[O-])COC(=O)CCCCCCC/C=C\\CCCCCCCC"]
  p mol.rule_name
end
