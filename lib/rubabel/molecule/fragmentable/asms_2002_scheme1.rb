# Add rules to these two variables as needed
# Rule_names << :rule_names
# Rearrangements << :any_rule_names_that_are_just_rearrangements
require_relative "../fragmentable"

module Rubabel
  class Molecule
    module Fragmentable
      ::Rule_names << def asms_2002_scheme1_a_c1(only_uniqs = true)
        # Create a fragmentation block method
        fragment_sets = []
        fragment = lambda do |tbcc, o, lc|
          # Duplication and mapping identity to new atoms
          nmol = self.dup
          to_be_carbonyl_carbon = nmol.atom(tbcc.id)
          oxygen = nmol.atom(o.id)
          link_carbon = nmol.atom(lc.id)
          # manipulate bonds
          nmol.delete_bond(to_be_carbonyl_carbon, link_carbon)
          to_be_carbonyl_carbon.get_bond(oxygen).bond_order = 2
          nmol.split
        end

        # Call the block on any search strings which make sense in the context of the code
        self.matches("C=CC([OH1])CN", only_uniqs).each do |carbon2_1, carbon2_2, to_be_carbonyl_carbon, oxygen, link_carbon, nitrogen|
          fragment_sets << fragment.call(to_be_carbonyl_carbon, oxygen,link_carbon)
        end
        fragment_sets.flatten
      end

      ::Rule_names << def asms_2002_scheme1_b_d1(only_uniqs = true)
        fragment_sets = []
        fragment = lambda do |n, tcc, o|
          #duplications and mapping
          nmol = self.dup
          nitrogen = nmol.atom(n.id)
          to_cyclize_carbon = nmol.atom(tcc.id)
          oxygen = nmol.atom(o.id)
          # manipulate bonds
          nmol.delete_bond(nitrogen, to_cyclize_carbon)
          nmol.add_bond!(to_cyclize_carbon, oxygen)
          nmol.split
        end

        # call the block search strings
        self.matches("N[Ch1]C[OH1]", only_uniqs).each do |nitrogen, to_cyclize_carbon, carbon, oxygen|
          fragment_sets << fragment.call(nitrogen, to_cyclize_carbon, oxygen)
        end
        fragment_sets.flatten
      end
      ::Rearrangements << def asms_2002_scheme1_b_d1_water_loss(only_uniqs = true)
        fragment_sets = []
        fragment = lambda do |lo, lc, lh|
          #duplications and mapping
          nmol = self.dup
          leaving_oxygen = nmol.atom(lo.id)
          link_carbon = nmol.atom(lc.id)
          losing_hydrogen_atom = nmol.atom(lh.id)
          # manipulate bonds
          bond_order = leaving_oxygen.get_bond(link_carbon).bond_order
          nmol.delete_bond(leaving_oxygen, link_carbon)
          losing_hydrogen_atom.get_bond(link_carbon).bond_order += bond_order
          nmol.split
        end

        # call the block search strings
        self.matches("[Oh1]C(C=C)[Ch1][OX2h0]", only_uniqs).each do |leaving_oxygen, link_carbon,carbon1, carbon2, losing_hydrogen_carbon, oxygen2|
          fragment_sets << fragment.call(leaving_oxygen, link_carbon, losing_hydrogen_carbon)
        end
        self.matches("O=C[Nh2]", only_uniqs).each do |leaving_oxygen, link_carbon, losing_hydrogen_nitrogen|
          fragment_sets << fragment.call(leaving_oxygen, link_carbon, losing_hydrogen_nitrogen)
        end
        fragment_sets.flatten
      end

      ::Rearrangements << def asms_2002_scheme1_b_d1_formaldehyde_loss(only_uniqs = true)
        fragment_sets = []
        fragment = lambda do |tco, leavingc, leftc|
          #duplications and mapping
          nmol = self.dup
          to_carbonyl_oxygen = nmol.atom(tco.id)
          leaving_carbon = nmol.atom(leavingc.id)
          left_carbon = nmol.atom(leftc.id)
          # manipulate bonds
          nmol.delete_bond(leaving_carbon, left_carbon)
          to_carbonyl_oxygen.get_bond(leaving_carbon).bond_order += 1
          nmol.split
        end

        # call the block search strings
     #     binding.pry if self.matches("[Oh1]C[CX4]([Oh0])[CX4]C=C", only_uniqs).size > 0
        self.matches("[Oh1]C[CX4]([Oh0])[CX4]C=C", only_uniqs).each do |to_carbonyl_oxygen, leaving_carbon,left_carbon, rest|
          fragment_sets << fragment.call(to_carbonyl_oxygen, leaving_carbon, left_carbon)
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
