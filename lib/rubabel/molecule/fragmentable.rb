require 'set'
require 'msplinter/fragment_set'

module Rubabel
  class Molecule
    module Fragmentable

      #RULES = Set[:cod, :codoo, :oxe, :oxepd, :oxh]
      RULES = Set[:cod, :codoo, :oxe, :oxepd, :oxh, :oxhpd]

      DEFAULT_OPTIONS = {
        rules: RULES,
        errors: :remove,
        # return only the set of unique fragments
        uniq: false, 
      }

      # molecules and fragments should all have hydrogens added (add_h!)
      # before calling this method
      # 
      # For instance, water loss with double bond formation is not allowable
      # for NCC(O)CC => CCC=C[NH2+], presumably because of the lone pair and
      # double bond resonance.
      def allowable_fragmentation?(frags)
        self.num_atoms(true) == frags.reduce(0) {|cnt,fr| cnt + fr.num_atoms(true) }
      end

      # splits the molecule between the carbon and carbon_nbr, adds a double
      # bond between the carbon and oxygen, and moves whatever was on the
      # oxygen (e.g., an OH or a charge) to the carbon_nbr. Returns two new
      # molecules.
      def carbonyl_oxygen_dump(carbon, oxygen, carbon_nbr)
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
        Msplinter::FragmentSet.new nmol.split
      end

      # breaks the bond and gives the electrons to the oxygen
      def carbon_oxygen_esteal(carbon, oxygen)
        nmol = self.dup
        ncarbon = nmol.atom(carbon.id)
        noxygen = nmol.atom(oxygen.id)

        is_carboxyl = noxygen.carboxyl_oxygen?
        
        nmol.delete_bond(ncarbon, noxygen)
        ncarbon.remove_a_hydride!
        noxygen.remove_a_proton! 
        Msplinter::FragmentSet.new nmol.split
      end

      # returns the duplicated molecule and the equivalent atoms
      def dup_molecule(atoms=[])
        nmol = self.dup
        [nmol, atoms.map {|old_atom| nmol.atom(old_atom.id) }]
      end

      # returns molecules created from splitting between the electrophile and
      # the center and where the bond order is increased between the center
      # and center_nbr
      def break_with_double_bond(electrophile, center, center_nbr)
        (nmol, (nele, ncarb, ncarb_nbr)) = self.dup_molecule([electrophile, center, center_nbr])
        nmol.delete_bond(nele, ncarb)
        ncarb_nbr.get_bond(ncarb) + 1
        Msplinter::FragmentSet.new nmol.split
      end

      def matches_rules?(*to_match, &block)
        if @rules.any? {|rule| to_match.include?(rule) }
          reply = block.call
        end
        reply
      end

      # carbon oxygen dump. returns an array of FragmentSets
      def cod
        frags = []
        matches_rules? :cod, :codoo do
          self.each_match("C[O;h1,O]", @only_uniqs) do |carbon, oxygen|
            carbon.atoms.select {|a| a.el == :C }.each do |carbon_nbr|
              frags << carbonyl_oxygen_dump(carbon, oxygen, carbon_nbr)
            end
          end
        end
        frags
      end

      def oxe
        frags = []
        matches_rules? :oxe do
          self.each_match("C-O", @only_uniqs) do |carbon, oxygen|
            frags << carbon_oxygen_esteal(carbon, oxygen)
          end
        end
        frags
      end

      def oxh
        # right now implemented so that a beta hydrogen has to be availabe for
        # extraction
        frags = []
        matches_rules? :oxh do
          self.each_match("C[C,O]-O", @only_uniqs) do |beta_c, center, oxygen|
            next unless beta_c.hydrogen_count > 0
            frags << break_with_double_bond(oxygen, center, beta_c)
          end
        end
        frags
      end

      def oxhpd
        frags = []
        matches_rules? :oxhpd do
          self.each_match("C-O-P-O", @only_uniqs) do |carbon, alc_oxy, phosphate, beta_carb_oxy|
            next unless beta_carb_oxy.hydrogen_count > 0
            frag_set = break_with_double_bond(alc_oxy, phosphate, beta_carb_oxy)
            frag_set.map! &:convert_dative_bonds!
            frags << frag_set
          end
        end
        frags
      end

      def oxepd
        frags = []
        matches_rules? :oxepd do
          self.each_match("P-O-C", @only_uniqs) do |phosphate, oxygen, carbon|
            frag_set = carbon_oxygen_esteal(phosphate, oxygen)
            frag_set.map! &:convert_dative_bonds!
            frags << frag_set
          end
        end
        frags
      end


      # an empty array is returned if there are no fragments generated.
      # Hydrogens are added at a pH of 7.4, unless they have already been
      # added.
      #
      #     :rules => an array of rules (e.g., [:cod, :oxe]), 
      #        but may also be an array of arrays where each array specifies 
      #        single, double, triple, etc fragmentation
      #        (e.g. [[:cod, :oxe], [:oxe]])
      #     :uniq => false
      #     :errors => :remove | :fix | :ignore  (default is :remove)
      def fragment(opts={})
        @only_uniqs = true
        opts = DEFAULT_OPTIONS.merge(opts)
        rule_sets = opts[:rules]
        rule_sets = [rule_sets] if rule_sets.first.is_a?(Symbol)
        rule_sets.flatten.each do |rule|
          raise ArgumentError, "bad rule: #{rule}" unless RULES.include?(rule)
        end

        had_hydrogens = self.h_added?
        self.correct_for_ph!(7.4) unless had_hydrogens
        self.remove_h!

        fragment_sets = []
        to_frag = [self]
        rule_sets.each do |rules|
          @rules = rules
          tmp_frags = []
          to_frag.each do |mol|
            [:cod, :codoo, :oxe, :oxh, :oxhpd, :oxepd].each do |methd|
              tmp_frags.push self.send(methd)
            end
          end
          fragment_sets << tmp_frags
          to_frag = tmp_frags.flatten(1)
        end

        case opts[:errors]
        when :remove
          fragment_sets.each do |frag_set|
            frag_set.select! {|set| allowable_fragmentation?(set) }
          end
        when :fix
          raise NotImplementedError
        when :ignore  # do nothing
        end

        self.remove_h!
        if opts[:uniq]
          # TODO: impelent properly
          raise NotImplementedError
          #fragment_sets = fragment_sets.uniq_by(&:csmiles)
        end

        fragment_sets
      end
    end
    include Fragmentable
  end
end # Rubabel