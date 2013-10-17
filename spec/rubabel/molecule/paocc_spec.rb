require 'spec_helper'

require 'rubabel'
require 'rubabel/molecule/fragmentable'
require 'pry'

$VERBOSE = nil

describe Rubabel::Molecule::Fragmentable do
  describe "Phosphate Attack On ester Carbon :paoc" do 
    it "produces cyclized product" do 
      @mol = Rubabel["CCCCCCCCCCCCCCCC(=O)O[C@@H](COP(=O)([O-])[O-])COC(=O)CCCCCCC/C=C\\CCCCCCCC"]
      @mol.write_file("paoc/root.svg")
      resp = @mol.fragment(rules: [:paoc])
      binding.pry
      p resp
      resp.size.>(0).should be_true
      # output image files for products
      resp.flatten.each_with_index{|mol,i| mol.write_file("paoc/#{mol.to_s.gsub("/", '-').gsub("\\", '_')}.svg") }

      # Test for the proper products
      resp.include?(Rubabel["[O-]C(=O)CCCCCCC/C=C/CCCCCCCC"]).should be_false
      resp.include?(Rubabel["[O-]C(=O)CCCCCCC/C=C\\CCCCCCCC"]).should be_true
      resp.flatten.include?(Rubabel["[O-]C(=O)CCCCCCCCCCCCCCC"]).should be_true
      resp.flatten.include?(Rubabel["[O-]P1(=O)OCC(COC(=O)CCCCCCC/C=C\\CCCCCCCC)O1"]).should be_true

      # Print out the mol_wt
      resp.flatten.each_with_index{|mol,i| p mol; p mol.mol_wt }
      #resp.flatten.include?(Rubabel["[O-]P1(=O)OCC(COC(=O)CCCCCCCC=CCCCCCCCC)O1"]).should be_true
      #failure "Why does it destroy the E/Z specificity in the product?"
    end
  end
end
