module Rubabel
  class FragmentationTree
    class Node
      attr_accessor :parent, :root, :molecule, :child
      MinimumSizeForFragmentation = 40.0
      DefaultFragmentOpts = {}
      def initialize(molecule, ms = [])
        @molecule = molecule
        @ms = ms
        @ms_level = 0
        if @ms.empty?
          @root = self
          @parent = self
        end
      end
      def check_for_double_fragmentations(ms_level)
        @ms[ms_level].map do |fragment|
          
        end
      end
      def ms2(opts= {})
        DefaultFragmentOpts.merge(opts)
        @ms << @molecule.fragment(opts)
        @ms_level = 1
        @ms[0]
      end
      def ms3(opts = {})
        ms2 if @ms.empty?
        @ms_level = 2
        @ms[1] = @ms[0].map {|i| i.map {|j| j.fragment(opts) } }.flatten
      end
      def msn(target_level=2) # perform MS(n) at specified level, changes the ms_level to indicate that level
        levels_to_traverse = target_level - @ms.size + 1
        @ms << @molecule.fragment(DefaultFragmentOpts) if @ms.empty? and @ms_level == 0
        while levels_to_traverse > 0
          @ms[@ms_level+1] = @ms[@ms_level].map {|i| i.map {|j| j.fragment(DefaultFragmentOpts) if j.mol_wt > MinimumSizeForFragmentation } }.flatten
          levels_to_traverse -= 1
        end
        @ms_level = target_level - 1
        @ms
      end
      def tree # Method for traversing all sizes and filling in the possible fragment list
        # Fill this in
      end
    end
  end
end
