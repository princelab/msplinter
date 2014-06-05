# Add rules to these two variables as needed
# Rule_names << :rule_names
# Rearrangements << :any_rule_names_that_are_just_rearrangements
require_relative "../fragmentable"

#        self.write( 'root.svg', :add_atom_index => true)
module Rubabel
  class Molecule
    module Fragmentable
      ::Rule_names << def jasms_2002_3_b1(only_uniqs: true, fragment_adduct_state: :as_published)
        fragment_sets = []
        fragment = lambda do |*arr|
          # duplicate
          (nmol, (alcohol_oxygen, alcohol_carbon, link_carbon, nitrogen, carbonyl, carbonyl_oxygen)) = self.dup_molecule(arr)
          # fragmentation chemistry
          nmol.add_bond!(alcohol_carbon, nitrogen)
          nmol.delete_bond(alcohol_carbon, alcohol_oxygen)
          nmol.split
        end
        self.matches("O[CH2]CNC=O" , only_uniqs).each do |alcohol_oxygen, alcohol_carbon, link_carbon, nitrogen, carbonyl, carbonyl_oxygen|
          fragment_sets << fragment.call(alcohol_oxygen, alcohol_carbon, link_carbon, nitrogen, carbonyl, carbonyl_oxygen)
        end
        if self.adducts.empty?
          fragment_sets.flatten
        else
          resp = case fragment_adduct_state
          when :force_adducts
            fragment_sets.flatten.map {|frag| frag.adducts.push(*self.adducts)}
            fragment_sets
          when :as_published
            # CONFIGURE THIS HERE BY Turning OFF the adduct line
            # ADDUCT LINE
            fragment_sets.flatten.map {|frag| frag.adducts.push(*self.adducts)}
            # This line must always be here
            fragment_sets.flatten
          when :no_adducts
            fragment_sets.flatten
          when :all
            dups = fragment_sets.flatten.map(&:dup) 
            fragment_sets.flatten.map {|frag| frag.adducts.push(*self.adducts)}
            dups + fragment_sets.flatten
          end
      end

      ::Rule_names << def jasms_2002_3_b2(only_uniqs: true, fragment_adduct_state: :as_published)
        fragment_sets = []
        fragment = lambda do |*arr|
          # duplicate
          (nmol, (alcc1, alco1, alcc2, alco2)) = self.dup_molecule(arr)
          # fragmentation chemistry
          alco1.get_bond(alcc1).bond_order = 2
          alco2.get_bond(alcc2).bond_order = 2
          nmol.delete_bond(alcc1, alcc2)
          nmol.split
        end
        self.matches("C(O)C(O)C(CO)NC(=O)", only_uniqs).each do |alcc1, alco1, alcc2, alco2, link_carbon, nitrogen, carbonyl, carbonyl_oxygen|
          fragment_sets << fragment.call(alcc1, alco1, alcc2, alco2)
        end
        if self.adducts.empty?
          fragment_sets.flatten
        else
          resp = case fragment_adduct_state
          when :force_adducts
            fragment_sets.flatten.map {|frag| frag.adducts.push(*self.adducts)}
            fragment_sets
          when :as_published
            # CONFIGURE THIS HERE BY Turning OFF the adduct line
            # ADDUCT LINE
            fragment_sets.flatten.map {|frag| frag.adducts.push(*self.adducts)}
            # This line must always be here
            fragment_sets.flatten
          when :no_adducts
            fragment_sets.flatten
          when :all
            dups = fragment_sets.flatten.map(&:dup) 
            fragment_sets.flatten.map {|frag| frag.adducts.push(*self.adducts)}
            dups + fragment_sets.flatten
          end
      end
      ::Rule_names << def jasms_2002_3_b2_water_loss(only_uniqs: true, fragment_adduct_state: :as_published)
        fragment_sets = []
        fragment = lambda do |*arr|
          # duplicate
          (nmol, (alcohol_oxygen, alcohol_carbon, nitrogen)) = self.dup_molecule(arr)
          # fragmentation chemistry
          nmol.delete_bond(alcohol_carbon, alcohol_oxygen)
          nmol.add_bond!(alcohol_carbon, nitrogen)
          nmol.split
        end
        self.matches("OCC(C=O)NC(=O)" , only_uniqs).each do |alcohol_oxygen, alcohol_carbon, link_carbon, carbonyl, carbonyl_oxygen, nitrogen, *rest|
          fragment_sets << fragment.call(alcohol_oxygen, alcohol_carbon, nitrogen)
        end
        if self.adducts.empty?
          fragment_sets.flatten
        else
          resp = case fragment_adduct_state
          when :force_adducts
            fragment_sets.flatten.map {|frag| frag.adducts.push(*self.adducts)}
            fragment_sets
          when :as_published
            # CONFIGURE THIS HERE BY Turning OFF the adduct line
            # ADDUCT LINE
            fragment_sets.flatten.map {|frag| frag.adducts.push(*self.adducts)}
            # This line must always be here
            fragment_sets.flatten
          when :no_adducts
            fragment_sets.flatten
          when :all
            dups = fragment_sets.flatten.map(&:dup) 
            fragment_sets.flatten.map {|frag| frag.adducts.push(*self.adducts)}
            dups + fragment_sets.flatten
          end
      end
      ::Rule_names << def jasms_2002_3_b2_formaldehyde_loss(only_uniqs: true, fragment_adduct_state: :as_published)
        fragment_sets = []
        fragment = lambda do |*arr|
          # duplicate
          (nmol, (alcohol_oxygen, alcohol_carbon, link_carbon)) = self.dup_molecule(arr)
          # fragmentation chemistry
          nmol.delete_bond(alcohol_carbon, link_carbon)
          alcohol_carbon.get_bond(alcohol_oxygen).bond_order = 2
          nmol.split
        end
        self.matches("OCC(C=O)NC(=O)" , only_uniqs).each do |alcohol_oxygen, alcohol_carbon, link_carbon, *rest|
          fragment_sets << fragment.call(alcohol_oxygen, alcohol_carbon, link_carbon)
        end
        if self.adducts.empty?
          fragment_sets.flatten
        else
          resp = case fragment_adduct_state
          when :force_adducts
            fragment_sets.flatten.map {|frag| frag.adducts.push(*self.adducts)}
            fragment_sets
          when :as_published
            # CONFIGURE THIS HERE BY Turning OFF the adduct line
            # ADDUCT LINE
            fragment_sets.flatten.map {|frag| frag.adducts.push(*self.adducts)}
            # This line must always be here
            fragment_sets.flatten
          when :no_adducts
            fragment_sets.flatten
          when :all
            dups = fragment_sets.flatten.map(&:dup) 
            fragment_sets.flatten.map {|frag| frag.adducts.push(*self.adducts)}
            dups + fragment_sets.flatten
          end
      end
      ::Rule_names << def jasms_2002_4_g1(only_uniqs: true, fragment_adduct_state: :as_published)
        fragment_sets = []
        fragment = lambda do |*arr|
          # duplicate
          (nmol, (alcohol_oxygen, alcohol_carbon, c1)) = self.dup_molecule(arr)
          # fragmentation chemistry
          nmol.delete_bond(alcohol_carbon, c1)
          alcohol_carbon.get_bond(alcohol_oxygen).bond_order = 2
          nmol.split
        end
        self.matches("OCC=CC(O)C(N)", only_uniqs).each do |alcohol_oxygen, alcohol_carbon, c1, c2, *rest|
          fragment_sets << fragment.call(alcohol_oxygen, alcohol_carbon, c1)
        end
        if self.adducts.empty?
          fragment_sets.flatten
        else
          resp = case fragment_adduct_state
          when :force_adducts
            fragment_sets.flatten.map {|frag| frag.adducts.push(*self.adducts)}
            fragment_sets
          when :as_published
            # CONFIGURE THIS HERE BY Turning OFF the adduct line
            # ADDUCT LINE
            fragment_sets.flatten.map {|frag| frag.adducts.push(*self.adducts)}
            # This line must always be here
            fragment_sets.flatten
          when :no_adducts
            fragment_sets.flatten
          when :all
            dups = fragment_sets.flatten.map(&:dup) 
            fragment_sets.flatten.map {|frag| frag.adducts.push(*self.adducts)}
            dups + fragment_sets.flatten
          end
      end
      ::Rule_names << def jasms_2002_4_g2a_prime(only_uniqs: true, fragment_adduct_state: :as_published)
        fragment_sets = []
        fragment = lambda do |*arr|
          # duplicate
          (nmol, (alcohol_oxygen, alcohol_carbon, nitrogen)) = self.dup_molecule(arr)
          # fragmentation chemistry
          nmol.add_bond!(alcohol_carbon, nitrogen)
          nmol.delete_bond(alcohol_carbon, alcohol_oxygen)
          nmol.split
        end
        self.matches("OCC(C(O)C=C)N", only_uniqs).each do |alcohol_oxygen, alcohol_carbon, c2, c3,o2,c4,c5, nitrogen|
          fragment_sets << fragment.call(alcohol_oxygen, alcohol_carbon, nitrogen)
        end
        if self.adducts.empty?
          fragment_sets.flatten
        else
          resp = case fragment_adduct_state
          when :force_adducts
            fragment_sets.flatten.map {|frag| frag.adducts.push(*self.adducts)}
            fragment_sets
          when :as_published
            # CONFIGURE THIS HERE BY Turning OFF the adduct line
            # ADDUCT LINE
            fragment_sets.flatten.map {|frag| frag.adducts.push(*self.adducts)}
            # This line must always be here
            fragment_sets.flatten
          when :no_adducts
            fragment_sets.flatten
          when :all
            dups = fragment_sets.flatten.map(&:dup) 
            fragment_sets.flatten.map {|frag| frag.adducts.push(*self.adducts)}
            dups + fragment_sets.flatten
          end
      end
      ::Rule_names << def jasms_2002_4_g2a(only_uniqs: true, fragment_adduct_state: :as_published)
        fragment_sets = []
        fragment = lambda do |*arr|
          # duplicate
          (nmol, (alcohol_oxygen, alcohol_carbon, acc1, aco1)) = self.dup_molecule(arr)
          # fragmentation chemistry
          nmol.delete_bond(acc1, aco1)
          nmol.add_bond!(alcohol_oxygen, acc1)
          nmol.split
        end
        self.matches("C=CC(O)C(CO)N", only_uniqs).each do |c1, c2, alcohol_carbon, alcohol_oxygen, link_carbon, acc1, aco1 , *rest|
          fragment_sets << fragment.call(alcohol_oxygen, alcohol_carbon, acc1, aco1) 
        end
        if self.adducts.empty?
          fragment_sets.flatten
        else
          resp = case fragment_adduct_state
          when :force_adducts
            fragment_sets.flatten.map {|frag| frag.adducts.push(*self.adducts)}
            fragment_sets
          when :as_published
            # CONFIGURE THIS HERE BY Turning OFF the adduct line
            # ADDUCT LINE
            fragment_sets.flatten.map {|frag| frag.adducts.push(*self.adducts)}
            # This line must always be here
            fragment_sets.flatten
          when :no_adducts
            fragment_sets.flatten
          when :all
            dups = fragment_sets.flatten.map(&:dup) 
            fragment_sets.flatten.map {|frag| frag.adducts.push(*self.adducts)}
            dups + fragment_sets.flatten
          end
      end
      ::Rule_names << def jasms_2002_4_g2a_formaldehyde_loss(only_uniqs: true, fragment_adduct_state: :as_published)
        fragment_sets = []
        fragment = lambda do |*arr|
          # duplicate
          (nmol, (ring_carbon, ring_oxygen, c3, link_carbon)) = self.dup_molecule(arr)
          # fragmentation chemistry
          nmol.delete_bond(ring_carbon, ring_oxygen)
          nmol.delete_bond(c3, link_carbon)
          ring_oxygen.get_bond(c3).bond_order = 2
          link_carbon.get_bond(ring_carbon).bond_order =2 
          nmol.split
        end
        self.matches("C=C[C;R;r4][O;R;r4]C[C;R;r4]NC(=O)", only_uniqs).each do |c1, c2, ring_carbon, ring_oxygen, c3, link_carbon, *rest|
          fragment_sets << fragment.call(ring_carbon, ring_oxygen, c3, link_carbon)
        end
        if self.adducts.empty?
          fragment_sets.flatten
        else
          resp = case fragment_adduct_state
          when :force_adducts
            fragment_sets.flatten.map {|frag| frag.adducts.push(*self.adducts)}
            fragment_sets
          when :as_published
            # CONFIGURE THIS HERE BY Turning OFF the adduct line
            # ADDUCT LINE
            fragment_sets.flatten.map {|frag| frag.adducts.push(*self.adducts)}
            # This line must always be here
            fragment_sets.flatten
          when :no_adducts
            fragment_sets.flatten
          when :all
            dups = fragment_sets.flatten.map(&:dup) 
            fragment_sets.flatten.map {|frag| frag.adducts.push(*self.adducts)}
            dups + fragment_sets.flatten
          end
      end
      ::Rule_names << def jasms_2002_4_g2a_prime_water_loss(only_uniqs: true, fragment_adduct_state: :as_published)
        fragment_sets = []
        fragment = lambda do |*arr|
          # duplicate
          (nmol, (alcohol_oxygen, alcohol_carbon, ring_carbon)) = self.dup_molecule(arr)
          # fragmentation chemistry
          nmol.delete_bond(alcohol_oxygen, alcohol_carbon)
          alcohol_carbon.get_bond(ring_carbon).bond_order = 2
          nmol.split
        end
        self.matches("OC[C;R;r3][C;R;r3][N;R;r3]", only_uniqs).each do |alcohol_oxygen, alcohol_carbon, ring_carbon, *rest|
          fragment_sets << fragment.call(alcohol_oxygen, alcohol_carbon, ring_carbon)
        end
        if self.adducts.empty?
          fragment_sets.flatten
        else
          resp = case fragment_adduct_state
          when :force_adducts
            fragment_sets.flatten.map {|frag| frag.adducts.push(*self.adducts)}
            fragment_sets
          when :as_published
            # CONFIGURE THIS HERE BY Turning OFF the adduct line
            # ADDUCT LINE
            fragment_sets.flatten.map {|frag| frag.adducts.push(*self.adducts)}
            # This line must always be here
            fragment_sets.flatten
          when :no_adducts
            fragment_sets.flatten
          when :all
            dups = fragment_sets.flatten.map(&:dup) 
            fragment_sets.flatten.map {|frag| frag.adducts.push(*self.adducts)}
            dups + fragment_sets.flatten
          end
      end
      ::Rule_names << def jasms_2002_5_e1_water_loss(only_uniqs: true, fragment_adduct_state: :as_published)
        fragment_sets = []
        fragment = lambda do |*arr|
          # duplicate
          (nmol, (alcohol_carbon, alcohol_oxygen, attacked_carbon, lost_oxygen)) = self.dup_molecule(arr)
          # fragmentation chemistry
          nmol.delete_bond(lost_oxygen, attacked_carbon)
          nmol.add_bond!(alcohol_oxygen, attacked_carbon)
          nmol.split
        end
        self.matches("C=CC(O)C(N)CO", only_uniqs).each do |c1, c2, alcohol_carbon, alcohol_oxygen, c3, nitrogen, attacked_carbon, lost_oxygen|
          fragment_sets << fragment.call(alcohol_carbon, alcohol_oxygen, attacked_carbon, lost_oxygen)
        end
        if self.adducts.empty?
          fragment_sets.flatten
        else
          resp = case fragment_adduct_state
          when :force_adducts
            fragment_sets.flatten.map {|frag| frag.adducts.push(*self.adducts)}
            fragment_sets
          when :as_published
            # CONFIGURE THIS HERE BY Turning OFF the adduct line
            # ADDUCT LINE
            fragment_sets.flatten.map {|frag| frag.adducts.push(*self.adducts)}
            # This line must always be here
            fragment_sets.flatten
          when :no_adducts
            fragment_sets.flatten
          when :all
            dups = fragment_sets.flatten.map(&:dup) 
            fragment_sets.flatten.map {|frag| frag.adducts.push(*self.adducts)}
            dups + fragment_sets.flatten
          end
      end
      ::Rule_names << def jasms_2002_5_e2(only_uniqs: true, fragment_adduct_state: :as_published)
        fragment_sets = []
        fragment = lambda do |*arr|
          # duplicate
          (nmol, (nitrogen, carbonyl, carbonyl_oxygen, c2)) = self.dup_molecule(arr)
          # fragmentation chemistry
          nmol.delete_bond(nitrogen, carbonyl)
          carbonyl.get_bond(c2).bond_order = 2
          nmol.split
        end
        self.matches("N(C(=O)C)C1COC1C=CC(O)", only_uniqs).each do |nitrogen, carbonyl, carbonyl_oxygen, c2|
          fragment_sets << fragment.call(nitrogen, carbonyl, carbonyl_oxygen, c2)
        end
        if self.adducts.empty?
          fragment_sets.flatten
        else
          resp = case fragment_adduct_state
          when :force_adducts
            fragment_sets.flatten.map {|frag| frag.adducts.push(*self.adducts)}
            fragment_sets
          when :as_published
            # CONFIGURE THIS HERE BY Turning OFF the adduct line
            # ADDUCT LINE
            fragment_sets.flatten.map {|frag| frag.adducts.push(*self.adducts)}
            # This line must always be here
            fragment_sets.flatten
          when :no_adducts
            fragment_sets.flatten
          when :all
            dups = fragment_sets.flatten.map(&:dup) 
            fragment_sets.flatten.map {|frag| frag.adducts.push(*self.adducts)}
            dups + fragment_sets.flatten
          end
      end
      ::Rule_names << def jasms_2002_5_e2_water_loss(only_uniqs: true, fragment_adduct_state: :as_published)
        fragment_sets = []
        fragment = lambda do |*arr|
          # duplicate
          (nmol, (c1, alcohol_carbon, alcohol_oxygen, c3, c4, nitrogen)) = self.dup_molecule(arr)
          # fragmentation chemistry
          nmol.delete_bond(alcohol_oxygen, alcohol_carbon)
          c1.get_bond(alcohol_carbon).bond_order = 2
          nitrogen.remove_a_proton!
          nmol.split
        end
        self.matches("CC(O)C=CC1OCC1N", only_uniqs).each do |c1, alcohol_carbon, alcohol_oxygen, c3, c4, *rest|
          fragment_sets << fragment.call(c1, alcohol_carbon, alcohol_oxygen, c3, c4, rest.last)
        end
        if self.adducts.empty?
          fragment_sets.flatten
        else
          resp = case fragment_adduct_state
          when :force_adducts
            fragment_sets.flatten.map {|frag| frag.adducts.push(*self.adducts)}
            fragment_sets
          when :as_published
            # CONFIGURE THIS HERE BY Turning OFF the adduct line
            # ADDUCT LINE
            fragment_sets.flatten.map {|frag| frag.adducts.push(*self.adducts)}
            # This line must always be here
            fragment_sets.flatten
          when :no_adducts
            fragment_sets.flatten
          when :all
            dups = fragment_sets.flatten.map(&:dup) 
            fragment_sets.flatten.map {|frag| frag.adducts.push(*self.adducts)}
            dups + fragment_sets.flatten
          end
      end
      ::Rule_names << def jasms_2002_5_e5_formaldehyde_loss(only_uniqs: true, fragment_adduct_state: :as_published)
        fragment_sets = []
        fragment = lambda do |*arr|
          # duplicate
          (nmol, (nitrogen, ring_carbon, lost_carbon, lost_oxygen, ring_carbon2)) = self.dup_molecule(arr)
          # fragmentation chemistry
          nmol.delete_bond(lost_oxygen, ring_carbon2)
          nmol.delete_bond(lost_carbon, ring_carbon)
          lost_oxygen.get_bond(lost_carbon).bond_order = 2
          ring_carbon.get_bond(ring_carbon2).bond_order = 2
          nitrogen.remove_a_proton!
          nmol.split
        end
        self.matches("N[C;R;r4][C;R;r4][O;R;r4][C;R;r4]C=CC=C", only_uniqs).each do |nitrogen, ring_carbon, lost_carbon, lost_oxygen, ring_carbon2, *rest|
          fragment_sets << fragment.call(nitrogen, ring_carbon, lost_carbon, lost_oxygen, ring_carbon2)
        end
        if self.adducts.empty?
          fragment_sets.flatten
        else
          resp = case fragment_adduct_state
          when :force_adducts
            fragment_sets.flatten.map {|frag| frag.adducts.push(*self.adducts)}
            fragment_sets
          when :as_published
            # CONFIGURE THIS HERE BY Turning OFF the adduct line
            # ADDUCT LINE
            fragment_sets.flatten.map {|frag| frag.adducts.push(*self.adducts)}
            # This line must always be here
            fragment_sets.flatten
          when :no_adducts
            fragment_sets.flatten
          when :all
            dups = fragment_sets.flatten.map(&:dup) 
            fragment_sets.flatten.map {|frag| frag.adducts.push(*self.adducts)}
            dups + fragment_sets.flatten
          end
	  end
      end
    end # Fragmentable
  end # Molecule
end

if $0 == __FILE__
  require 'bundler/setup'
  require 'rubabel'
  require 'pry'
  mol = Rubabel["LMGP04010962", :lmid]
  mol = Rubabel["CCCCCCCCCCCCCCCC(=O)OC(COP(=O)([O-])[O-])COC(=O)CCCCCCC/C=C\\CCCCCCCC"]
  p mol.rule_name
end
