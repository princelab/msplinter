require 'spec_helper'

require 'rubabel'
require 'rubabel/molecule/fragmentable'
require 'pry'

$VERBOSE = nil

PRINT_IMAGES = false
PRINT_MASSES = false

describe Rubabel::Molecule::Fragmentable do
  describe "Nitrogen Cyclization :nc" do 
    DEBUGGING_FOLDER = "nitrogen"
    it "produces cyclized product" do 
      @mol = Rubabel["LMSP02030004", :lmid]
      #@mol.write_file("#{DEBUGGING_FOLDER}/ceramide18-0_24-0.svg", add_atom_index: true)
      resp = @mol.fragment(rules: [:nc])
      resp.size.>(0).should be_true
      resp = resp.flatten.uniq{|a| a.csmiles}
      # output image files for products
      if PRINT_IMAGES
        resp.each_with_index{|mol,i| mol.write_file("#{DEBUGGING_FOLDER}/#{mol.to_s.gsub("/", '-').gsub("\\", '_')}.svg", add_atom_index: true) }
      end

      # Test for the proper products
      resp.include?(Rubabel["CCCCCCCCCCCCCCCCCCCCCCCC(=O)N1[C@H](C1[C@@H](CCCCCCCCCCCCCC)O)CO"]).should be_true
      resp.include?(Rubabel["CCCCCCCCCCCCCCCCCCCCCCCC(=O)N1C[C@H]1[C@@H]([C@@H](CCCCCCCCCCCCCC)O)O"]).should be_true
      resp.include?(Rubabel["O"]).should be_true
      #TODO write a function to allow for quick product SMARTS searching.should_smarts("r3")
      
      # Print out the mol_wt
      if PRINT_MASSES
        resp.each_with_index{|mol,i| p mol; p mol.mol_wt }
      end
    end
    it "method #2 gives cyclized product" do 
      @mol = Rubabel["OCC(N)C(O)C=CCCCCCCCCCCCCC"]
      p @mol.mass
      @mol.write_file("#{DEBUGGING_FOLDER}/test.svg", add_atom_index: true)
      resp = @mol.fragment(rules: [:nc])
      resp.size.>(0).should be_true
      resp = resp.flatten.uniq{|a| a.csmiles}
      # output image files for products
      if PRINT_IMAGES
        resp.each_with_index{|mol,i| mol.write_file("#{DEBUGGING_FOLDER}/#{mol.to_s.gsub("/", '-').gsub("\\", '_')}.svg", add_atom_index: true) }
      end
      p resp

      # Test for the proper products
      resp.include?(Rubabel["CCCCCCCCCCCCC/C=C\\C(C1NC1)O"]).should be_true
      resp.include?(Rubabel["CCCCCCCCCCCCC/C=C\\C1NC1CO"]).should be_true
      resp.include?(Rubabel["O"]).should be_true

      # Print out the mol_wt
      if PRINT_MASSES
        resp.each_with_index{|mol,i| p mol; p mol.mol_wt }
      end
    end
  end
end
