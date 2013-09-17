require 'rubabel/molecule/fragmentable'
require 'rubabel/fragmentation_tree/node'

module Rubabel
  class FragmentationTree
    include Enumerable # What does this buy me?

    attr_reader :root
    def initialize(*args)
      @root = Node.new(*args)
    end
    def push(*args)
      tmp = Node.new(*args)
      @root.child = tmp
      tmp.parent = @root
    end
    alias :<< :push 
  end
end
