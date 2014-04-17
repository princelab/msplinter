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
      ::Rearrangements << def jasms_2002_2_f2_bprime_water_loss(only_uniqs=true)
        fragment_sets = []
        fragment = lambda do |*arr|
          # duplicate
          (nmol, (linked_carbon, alcohol_carbon, oxygen)) = self.dup_molecule(arr)
          # fragmentation chemistry
          nmol.delete_bond(alcohol_carbon, oxygen)
          alcohol_carbon.get_bond(linked_carbon).bond_order = 2
          nmol.split
        end
        self.matches("[N;R;r3][CH2;R;r3][CH1;R;r3]C(O)C=C" , only_uniqs).each do |nitrogen, ring_carbon, linked_carbon, alcohol_carbon, oxygen, *rest|
          fragment_sets << fragment.call(linked_carbon, alcohol_carbon, oxygen)
        end
        fragment_sets.flatten
      end
      ::Rearrangements << def jasms_2002_2_f2_b_methanol_loss(only_uniqs=true)
        fragment_sets = []
        fragment = lambda do |*arr|
          # duplicate
          (nmol, (c1, lost_carbon, oxygen, c2)) = self.dup_molecule(arr)
          # fragmentation chemistry
          nmol.delete_bond(c1, lost_carbon)
          nmol.delete_bond(oxygen, c2)
          c1.get_bond(c2).bond_order += 1
          nmol.split
        end
        self.matches("[C;R;r4][CH2;R;r4][O;R;r4][C;R;r4]", only_uniqs).each do |c1, lost_carbon, oxygen, c2|
          fragment_sets << fragment.call(c1, lost_carbon, oxygen, c2)
        end
        fragment_sets.flatten
      end
      ::Rule_names << def jasms_2002_2_e1(only_uniqs=true)
        fragment_sets = []
        fragment = lambda do |*arr|
          # duplicate
          (nmol, (nitrogen, carbonyl_carbon, carbonyl_oxygen, alcohol_carbon, alcohol_oxygen)) = self.dup_molecule(arr)
          # fragmentation chemistry
          nmol.delete_bond(nitrogen, carbonyl_carbon)
          nmol.delete_bond(carbonyl_carbon, alcohol_carbon)
          alcohol_carbon.get_bond(alcohol_oxygen).bond_order = 2
          carbonyl_oxygen.get_bond(carbonyl_carbon).bond_order = 3
          carbonyl_oxygen.remove_a_proton!
          nmol.split
        end
        self.matches("NC(=O)C(O)", only_uniqs).each do |nitrogen, carbonyl_carbon, carbonyl_oxygen, alcohol_carbon, alcohol_oxygen|
          fragment_sets << fragment.call(nitrogen, carbonyl_carbon, carbonyl_oxygen, alcohol_carbon, alcohol_oxygen)
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
