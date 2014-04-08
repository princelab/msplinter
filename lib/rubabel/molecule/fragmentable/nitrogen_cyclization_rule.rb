# Add rules to these two variables as needed
# Rule_names << :rule_names
# Rearrangements << :any_rule_names_that_are_just_rearrangements
require_relative "../fragmentable"

module Rubabel
  class Molecule
    module Fragmentable
      ::Rule_names << def nitrogen_cyclization
        #TODO remove this next line
        only_uniqs = true
        fragment_sets = []
        # Create a fragmentation block method
        fragment = lambda do |n, tbc, lo|
          # Splits the molecule by cyclization of the nitrogen with loss of water
          (nmol, (nitrogen, to_bond_carbon, lost_oxygen)) = dup_molecule([n, tbc, lo])
          original_n_charge = nitrogen.charge
          nmol.add_bond!(nitrogen, to_bond_carbon)
          nmol.delete_bond(lost_oxygen.get_bond(to_bond_carbon))
          nitrogen.remove_a_proton!
          nitrogen.charge= original_n_charge
          nmol.split
        end

        # Call the block on any search strings which make sense in the context of the code
        self.matches("[CX3](=O)[NX3h1]C[CX4][OX2H1]", only_uniqs).each do |arr|
          puts "C=ONCCOH"
          carbonyl_carbon = arr.first
          carbonyl_oxygen = arr[1]
          nitrogen = arr[2]
          link_carbon = arr[3]
          to_bond_carbon = arr[4]
          leaving_oxygen = arr[5]
          fragment_sets << fragment.call(nitrogen, to_bond_carbon, leaving_oxygen)
        end
        self.matches("OCC[NH2,NH1,NH3]", only_uniqs).each do |arr| # works for first case, not second
          leaving_oxygen = arr.first  
          to_bond_carbon = arr[1]
          nitrogen = arr[3]
          fragment_sets << fragment.call(nitrogen, to_bond_carbon, leaving_oxygen)
        end
        fragment_sets
      end
      ::Rearrangements << ::Rule_names.last

    end # Fragmentable
  end # Molecule
end #Rubabel


if $0 == __FILE__
  require 'bundler/setup'
  require 'rubabel'
  mol = Rubabel["LMGP04010962", :lmid]
  mol = Rubabel["CCCCCCCCCCCCCCCC(=O)OC(COP(=O)([O-])[O-])COC(=O)CCCCCCC/C=C\\CCCCCCCC"]
  mol = Rubabel["OCC(N)C(O)C=CCCCCCCCCCCCCC"]
  p mol.nitrogen_cyclization
end
