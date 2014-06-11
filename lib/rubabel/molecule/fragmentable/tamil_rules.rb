# Add rules to these two variables as needed
# Rule_names << :rule_names
# Rearrangements << :any_rule_names_that_are_just_rearrangements
require_relative "../fragmentable"

module Rubabel
  class Molecule
    module Fragmentable
      ::Rule_names << def break_phosphates_with_double_bond(only_uniqs: true, fragment_adduct_state: :as_published)
        fragment_sets = []
        only_uniqs = true
        fragment = lambda do |electrophile, center, center_nbr|
          (nmol, (nele, ncarb, ncarb_nbr)) = self.dup_molecule([electrophile, center, center_nbr])
          nmol.delete_bond(nele, ncarb)
          ncarb_nbr.get_bond(ncarb) + 1
          nmol.split
        end
        # Search strings
        self.matches("C-O-P-[Oh1]", only_uniqs) do |carbon, alc_oxy, phosphate, beta_carb_oxy|
            frag_set = fragment.call(alc_oxy, phosphate, beta_carb_oxy)
            frag_set.map! &:convert_dative_bonds!
            fragment_sets << frag_set
          end
        fragment_sets
      end
      ::Rule_names << def break_oxygens_with_double_bond(only_uniqs: true, fragment_adduct_state: :as_published)
        fragment_sets = []
        only_uniqs = true
        fragment = lambda do |electrophile, center, center_nbr|
          (nmol, (nele, ncarb, ncarb_nbr)) = self.dup_molecule([electrophile, center, center_nbr])
          nmol.delete_bond(nele, ncarb)
          ncarb_nbr.get_bond(ncarb) + 1
          nmol.split
        end
        # Search strings
        self.matches("[Ch1][C,O]-O", only_uniqs) do |beta_c, center, oxygen|
          fragment_sets << fragment.call(oxygen, center, beta_c)
        end
        fragment_sets
      end
      ::Rule_names << def carbon_oxygen_esteal
        fragment_sets = []
        only_uniqs = true
        fragment = lambda do |carbon, oxygen|
          nmol = self.dup
          ncarbon = nmol.atom(carbon.id)
          noxygen = nmol.atom(oxygen.id)
          nmol.delete_bond(ncarbon, noxygen)
          ncarbon.remove_a_hydride!
          noxygen.remove_a_proton! 
          nmol.split
        end
        self.each_match("C-O", only_uniqs) do |carbon, oxygen|
          fragment_sets << fragment.call(carbon, oxygen)
        end
        fragment_sets
      end
      ::Rule_names << def phosphate_oxygen_esteal(only_uniqs: true, fragment_adduct_state: :as_published)
        fragment_sets = []
        only_uniqs = true
        fragment = lambda do |carbon, oxygen|
          nmol = self.dup
          ncarbon = nmol.atom(carbon.id)
          noxygen = nmol.atom(oxygen.id)
          nmol.delete_bond(ncarbon, noxygen)
          ncarbon.remove_a_hydride!
          noxygen.remove_a_proton! 
          nmol.split
        end
        self.each_match("P-O-C", only_uniqs) do |phosphate, oxygen, carbon|
          frag_set = fragment.call(phosphate, oxygen)
          frag_set.map! &:convert_dative_bonds!
          fragment_sets << frag_set
        end

        fragment_sets
      end
      ::Rule_names << def carbonyl_oxygen_dump(only_uniqs: true, fragment_adduct_state: :as_published)
        fragment_sets = []
        only_uniqs = true
        fragment = lambda do |carbon, oxygen, carbon_nbr|
          appendage = oxygen.atoms.find {|a| a.el != :C }
          if oxygen.charge != 0
            ocharge = oxygen.charge
          end
          nmol = self.dup
          new_oxygen = nmol.atom(oxygen.id)
          new_carbon = nmol.atom(carbon.id)
          new_carbon_nbr = nmol.atom(carbon_nbr.id)
          new_appendage = nmol.atom(appendage.id) if appendage
          nmol.delete_bond(new_carbon.get_bond(new_carbon_nbr))
          if new_appendage
            nmol.delete_bond(new_oxygen.get_bond(new_appendage)) 
            nmol.add_bond!(new_carbon_nbr, new_appendage)
          end
          if ocharge
            new_carbon_nbr.charge += ocharge
            new_oxygen.charge -= ocharge
          end
          new_carbon.get_bond(new_oxygen).bond_order = 2
          nmol.split
        end
        self.each_match("C[O;h1,O]", only_uniqs) do |carbon, oxygen|
          carbon.atoms.select {|a| a.el == :C }.each do |carbon_nbr|
            fragment_sets << fragment.call(carbon, oxygen, carbon_nbr)
          end
        end
        fragment_sets
      end


    end # Fragmentable
  end # Molecule
end #Rubabel


if $0 == __FILE__
  require 'bundler/setup'
  require 'rubabel'
  require 'pry'
  mol = Rubabel["LMGP04010962", :lmid]
  p mol.fragment 
end
