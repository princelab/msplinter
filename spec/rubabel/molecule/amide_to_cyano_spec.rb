require 'spec_helper'

require 'rubabel'
require 'rubabel/molecule/fragmentable'
require 'pry'

$VERBOSE = nil

PRINT_IMAGES = false
PRINT_MASSES = false

describe Rubabel::Molecule::Fragmentable do
  describe "Amide to Cyano dehydration Rule :atcd" do 
    DEBUGGING_FOLDER = "cyano"
    before(:each) do 
      @mol = Rubabel["NC(=O)CCCCCCCCCCCCCCCCCCCCCCC"]
      @fragments = @mol.fragment(rules: [:atcd])
    end
    it "generates fragments" do 
      @fragments.size.>(0).should be_true
    end
    it "makes a cyano" do 
      @fragments.flatten.should_smarts("C#N")
    end
    it "gives expected products" do 
      #@mol.write_file("#{DEBUGGING_FOLDER}/ceramide18-0_24-0.svg", add_atom_index: true)
      resp = resp.flatten.uniq{|a| a.csmiles}
      # output image files for products
      if PRINT_IMAGES
        resp.each_with_index{|mol,i| mol.write_file("#{DEBUGGING_FOLDER}/#{mol.to_s.gsub("/", '-').gsub("\\", '_')}.svg", add_atom_index: true) }
      end
      # Test for the proper products
      resp.include?(Rubabel["N#CCCCCCCCCCCCCCCCCCCCCCCC"]).should be_true
      resp.include?(Rubabel["O"]).should be_true
      #TODO write a function to allow for quick product SMARTS searching.should_smarts("r3")
      
      # Print out the mass
      if PRINT_MASSES
        resp.each_with_index{|mol,i| p mol; p mol.mass }
      end
    end
  end
end
