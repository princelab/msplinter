require 'set'
Rule_names = []
Rearrangements = []

Dir.glob(File.join(File.dirname(File.absolute_path(__FILE__)),"fragmentable", "*.rb")).map {|rfile| require rfile }

module Rubabel
  class Molecule
    module Fragmentable

      #RULES = Set[:cod, :codoo, :oxe, :oxepd, :oxh]
      #RULES = Set[:cod, :codoo, :oxe, :oxepd, :oxh, :oxhpd, :paoc, :nc] # :paoc is Phosphate Attack On Carbonyl Carbon
      @@rules = Set[*Rule_names]
      @@rearrangements = Set[*Rearrangements]
      def rules
        @@rules
      end
      RULES = @@rules
      REARRANGEMENTS = @@rearrangements

      DEFAULT_OPTIONS = {
        rules: @@rules - @@rearrangements, 
        errors: :remove,
        # return only the set of unique fragments
        uniq: false,
        fragment_adduct_state: :as_published
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

      # returns the duplicated molecule and the equivalent atoms
      def dup_molecule(atoms=[])
        nmol = self.dup
        [nmol, atoms.map {|old_atom| nmol.atom(old_atom.id) }]
      end

      # An API for determining if a fragmentation rule should occur given configurable
      # parameters.  This takes into account adducts, charge state, instrument type... etc
      # Rubabel::FragmentationConditions?
      def instrument_configuration(rules, rearrangements)
      end

      # An API to select rules based upon the chemistry of the given molecule... ?  Maybe... 
      def select_rules_by_molecule
      end

      # an empty array is returned if there are no fragments generated.
      # Hydrogens are added at a pH of 7.4, unless they have already been
      # added.
      #
      #     :rules => queryable by :include? set of rules
      #     :uniq => false
      #     :errors => :remove | :fix | :ignore  (default is :remove)
      #     :fragment_adduct_state => :force_adducts, :as_published (default), :no_adducts, :all (doubles fragment list size)
      def fragment(rules: @@rules- @@rearrangements, errors: :ignore, uniq: false, rearrange: true, fragment_adduct_state: DEFAULT_OPTIONS[:fragment_adduct_state], add_rule_names: false)
        only_uniqs = true # Option is currently ignored.
        # opts = DEFAULT_OPTIONS.merge(opts)
        rules.each do |rule| 
          raise ArgumentError, "bad rule: #{rule}\nThe allowed rules are #{@@rules.entries.join(", ")}" unless @@rules.include?(rule)
        end

        had_hydrogens = self.h_added?
        self.correct_for_ph!(7.4) unless had_hydrogens
        self.remove_h!

        fragment_sets = {}
        fragments = []
        # call on API for selection rules from given experimental conditions
        #
        #
        # Fxn call to select rules... 

        # Call all rules for now... 
        rules.map do |rule|
          rule_fragments = self.send(rule, only_uniqs: only_uniqs, fragment_adduct_state: fragment_adduct_state)
          #binding.pry if rule_fragments.nil?
          rule_fragments = rule_fragments.flatten.compact
          if self.adducts? 
          #  rule_fragments.map {|m| m.adducts = self.adducts }  # TODO THIS IS BAD, it doesn't allow for rules which don't produce adduct added fragments.  I was wrong and need to refactor the rules to account for this.
          end
          fragment_sets[rule] = rule_fragments.flatten.compact
          fragments << rule_fragments
        end


        if rearrange
          (::Rearrangements+::Rule_names).map do |rule|
            fragments.flatten.map do |fragment|
              rule_fragments = fragment.send(rule)
              if fragment_sets[rule]
                fragment_sets[rule] = (fragment_sets[rule] + rule_fragments).flatten
              else
                fragment_sets[rule] = rule_fragments
              end
              fragments << rule_fragments
            end
          end
        end

        # Error handling
        case :errors
        when :remove
          # This can't work like it used to, as I've already flattened away the meaning of the set
          puts "This is like ignore right now"
          #fragments.select! {|set| allowable_fragmentation?(set) }
        when :fix
          raise NotImplementedError
        when :ignore  # do nothing
        end

        self.remove_h!
        if uniq
          # TODO: implement properly
          raise NotImplementedError
          #fragment_sets = fragment_sets.uniq_by(&:csmiles)
        end

        if add_rule_names
          fragment_sets # Hash of fragments by rule
        else
          fragments.select {|a| not a.empty? } #Array of all fragments
        end
      end
      def rearrange(rules: @@rearrangements, errors: :remove, uniq: false)
        only_uniqs = uniq
        had_hydrogens = self.h_added?
        self.correct_for_ph!(7.4) unless had_hydrogens
        self.remove_h!
        fragment_sets = Hash.new {|h,k| h[k] = [] }
        fragments = []
        rules.map do |rule|
          rule_fragments = self.send(rule)
          fragment_sets[rule] << rule_fragments.flatten.compact
          fragments << rule_fragments
        end
        fragment_sets.each {|k, v| fragment_sets[k] = v.flatten }
        fragments
      end #rearrange  
      def rearrange_but_return_sets(rules: @@rearrangements, errors: :remove, uniq: false)
        only_uniqs = uniq
        had_hydrogens = self.h_added?
        self.correct_for_ph!(7.4) unless had_hydrogens
        self.remove_h!
        fragment_sets = Hash.new {|h,k| h[k] = [] }
        fragments = []
        rules.map do |rule|
          rule_fragments = self.send(rule)
          fragment_sets[rule] << rule_fragments.flatten.compact
          fragments << rule_fragments
        end
        fragment_sets
      end #rearrange
    end
    include Fragmentable
  end
end # Rubabel
