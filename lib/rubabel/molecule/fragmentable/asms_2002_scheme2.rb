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
      ::Rearrangements << def jasms_2002_2_f2_b_methoxy_loss(only_uniqs=true)
        fragment_sets = []
        fragment = lambda do |*arr|
          # duplicate
          (nmol, (c1, lost_carbon, oxygen, c2)) = self.dup_molecule(arr)
          # fragmentation chemistry
          nmol.delete_bond(c1, lost_carbon)
          nmol.delete_bond(oxygen, c2)
          c1.get_bond(c2).bond_order += 1
          lost_carbon.remove_a_proton!
          lost_carbon.remove_a_proton!
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
      ::Rule_names << def jasms_2002_2_e1_b(only_uniqs=true)
        fragment_sets = []
        fragment = lambda do |*arr|
          # duplicate
          (nmol, (alcohol_oxygen, nitrogen, attacked_carbon, leaving_oxygen)) = self.dup_molecule(arr)
          # fragmentation chemistry
          nmol.delete_bond(leaving_oxygen, attacked_carbon)
          nmol.add_bond!(alcohol_oxygen, attacked_carbon)
          nitrogen.remove_a_proton!
          nmol.split
        end
        self.matches("OCC(N)C(O)", only_uniqs).each do |alcohol_oxygen, c1, c2, nitrogen, attacked_carbon, leaving_oxygen|
          fragment_sets << fragment.call(alcohol_oxygen, nitrogen, attacked_carbon, leaving_oxygen)
        end
        fragment_sets.flatten
      end
      ::Rule_names << def jasms_2002_2_e1_bprime(only_uniqs=true)
        fragment_sets = []
        fragment = lambda do |*arr|
          # duplicate
          (nmol, (alcohol_oxygen, attacked_carbon, carbon_linker, nitrogen)) = self.dup_molecule(arr)
          # fragmentation chemistry
          nmol.delete_bond(alcohol_oxygen, attacked_carbon)
          nmol.add_bond!(nitrogen, attacked_carbon)
          nitrogen.remove_a_proton!
          nmol.split
        end
        self.matches("OCC(N)C(O)C=C", only_uniqs).each do |alcohol_oxygen, attacked_carbon, carbon_linker, nitrogen, *rest|
          fragment_sets << fragment.call(alcohol_oxygen, attacked_carbon, carbon_linker, nitrogen)
        end
        fragment_sets.flatten
      end
      ::Rule_names << def jasms_2002_2_e2_b_formaldehyde_loss(only_uniqs=true)
        fragment_sets = []
        fragment = lambda do |*arr|
          # duplicate
          (nmol, (nitrogen, c1, c2, ester_oxygen, c3)) = self.dup_molecule(arr)
          # fragmentation chemistry
          c3.get_bond(c1).bond_order = 2
          nmol.delete_bond(c1, c2)
          nmol.delete_bond(ester_oxygen, c3)
          c2.get_bond(ester_oxygen).bond_order = 2
          nitrogen.remove_a_proton!
          nmol.split
        end
        self.matches("NC1COC1C=C", only_uniqs).each do |nitrogen, c1, c2, ester_oxygen, c3, *rest|
          fragment_sets << fragment.call(nitrogen, c1, c2, ester_oxygen, c3)
        end
        fragment_sets.flatten
      end
      ::Rule_names << def jasms_2002_2_e2_bprime_water_loss(only_uniqs=true)
        fragment_sets = []
        fragment = lambda do |*arr|
          # duplicate
          (nmol, (alcohol_oxygen, c1, c2, nitrogen)) = self.dup_molecule(arr)
          # fragmentation chemistry
          nmol.delete_bond(alcohol_oxygen, c1)
          c1.get_bond(c2).bond_order = 2
          nitrogen.remove_a_proton!
          nmol.split
        end
        self.matches("OCC1CN1", only_uniqs).each do |alcohol_oxygen, c1, c2, c3, nitrogen|
          fragment_sets << fragment.call(alcohol_oxygen, c1, c2, nitrogen)
        end
        fragment_sets.flatten
      end
      ::Rule_names << def jasms_2002_2_e2_bprime_heterocyclic_loss(only_uniqs=true)
        fragment_sets = []
        fragment = lambda do |*arr|
          # duplicate
          (nmol, (nitrogen, ring_carbon, alcohol_carbon, alcohol_oxygen)) = self.dup_molecule(arr)
          # fragmentation chemistry
          nmol.delete_bond(ring_carbon, alcohol_carbon)
          alcohol_carbon.get_bond(alcohol_oxygen).bond_order = 2
          nitrogen.remove_a_proton!
          nmol.split
        end
        self.matches("N1CC1C(O)C=C", only_uniqs).each do |nitrogen, c1, ring_carbon, alcohol_carbon, alcohol_oxygen, *rest|
          fragment_sets << fragment.call(nitrogen, ring_carbon, alcohol_carbon, alcohol_oxygen)
        end
        fragment_sets.flatten
      end
      ::Rule_names << def jasms_2002_2_c2a(only_uniqs=true)
        fragment_sets = []
        fragment = lambda do |*arr|
          # duplicate
          (nmol, (nitrogen, carbonyl_carbon, carbonyl_oxygen, alcohol_carbon, alcohol_oxygen)) = self.dup_molecule(arr)
          # fragmentation chemistry
          nmol.delete_bond(nitrogen, carbonyl_carbon)
          nmol.delete_bond(carbonyl_carbon, alcohol_carbon)
          alcohol_carbon.get_bond(alcohol_oxygen).bond_order = 2
          nmol.split
        end
        self.matches("[N;R;r3]C(=O)C(O)C", only_uniqs).each do |nitrogen, carbonyl_carbon, carbonyl_oxygen, alcohol_carbon, alcohol_oxygen, *rest|
          fragment_sets << fragment.call(nitrogen, carbonyl_carbon, carbonyl_oxygen, alcohol_carbon, alcohol_oxygen)
        end
        fragment_sets.flatten
      end
      ::Rule_names << def jasms_2002_2_a1(only_uniqs=true)
        fragment_sets = []
        fragment = lambda do |*arr|
          # duplicate
          (nmol, (base_c, alcohol_carbon, alcohol_oxygen, carbonyl_carbon, carbonyl_oxygen, nitrogen, ring_carbon)) = self.dup_molecule(arr)
          # fragmentation chemistry
          nitrogen.get_bond(carbonyl_carbon).bond_order = 2
          nmol.delete_bond(carbonyl_carbon, alcohol_carbon)
          alcohol_carbon.get_bond(alcohol_oxygen).bond_order = 2
          nmol.delete_bond(nitrogen, ring_carbon)
          nmol.split
        end
        self.matches("CC(O)C(=O)N[C;R;r4]", only_uniqs).each do |base_c, alcohol_carbon, alcohol_oxygen, carbonyl_carbon, carbonyl_oxygen, nitrogen, ring_carbon|
          fragment_sets << fragment.call(base_c, alcohol_carbon, alcohol_oxygen, carbonyl_carbon, carbonyl_oxygen, nitrogen, ring_carbon)
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
