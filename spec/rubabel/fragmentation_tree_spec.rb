require 'spec_helper'

require 'rubabel'
require 'rubabel/fragmentation_tree.rb'

require 'pry'

$VERBOSE = nil
FragmentationTree = Rubabel::FragmentationTree
$molecule_test_string = "LMGP04010962"
describe Rubabel::FragmentationTree do 
  before :each do 
    @a =  Rubabel::FragmentationTree::Node.new Rubabel[$molecule_test_string, :lmid]
    @tree = FragmentationTree.new @a
  end
  it "Initializes as its own root molecule" do 
    @a.root.should == @a
  end
  it "holds a parent" do 
    @a.parent.should == @a
  end
  it "handles ms2" do 
    resp = @a.ms2(rules: [:cod])
    resp.should == Rubabel[$molecule_test_string, :lmid].fragment(rules: [:cod])
  end
  it "handles ms3 for everything" do 
    resp = @a.ms3(rules: [:cod, :oxe, :codoo, :oxepd, :oxh, :oxhpd])
    resp.include?(Rubabel["[O-]C(=O)CCCCCCCCCCCCCCC"]).should be_true
    resp.include?(Rubabel["[O-]C(=O)CCCCCCC/C=C\\CCCCCCCC"]).should be_true
    # DOESN"T YET WORK
    binding.pry
    #resp.include?(Rubabel["[O-]P1(=O)OCC(OC(=O)CCCCCCCCCCCCCCCC)CO1"]).should be_true
    #resp.include?(Rubabel["CCCCCCCC/C=C\\CCCCCCCC(=O)OCC1CO[P]([O-])(=O)O1"]).should be_true
  end
  it "it traverses the stack, in order to perform all possible fragmentations" do 
    resp = @a.msn
    @a.molecule.mol_wt
    resp.flatten.compact.map(&:mol_wt).sort.uniq
  end
  it "Doesn't create the false product from this lipid test case (loss of a carbonyl from a glycerophosplipid in the middle of the chain" do
    @a.ms2.include?(Rubabel["[C@@H](P(=O)([O-])OC[C@@H](O)CO)(OC(=O)CCCCCCCCCCCCCCC)COC(=O)CCCCCCC/C=C\\CCCCCCCC"]).should be_false
    @a.ms2.include?(Rubabel["[C@@H](P(=O)([O-])OC[C@@H](O)CO)(OC(=O)CCCCCCCCCCCCCCC)COC(=O)CCCCCCC/C=C\\CCCCCCCC"]).should be_false
  end
end

