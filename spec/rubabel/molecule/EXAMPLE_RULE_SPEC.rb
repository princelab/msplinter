require 'spec_helper'

require 'rubabel'
require 'rubabel/molecule/fragmentable'

$VERBOSE = nil

describe Rubabel::Molecule::Fragmentable do

  describe 'rule name here' do 
    # Comment here
    specify 'rule: test_case' do
      mol = Rubabel["NCC(O)CC"]
      frags = mol.fragment(rules: [:cod])
      frags.flatten(1).map{|a| a.map(&:csmiles) }.flatten.should == ["C[NH3+]", "CCC=O", "C([NH3+])C=O", "CC"]
    end

    specify 'rule: test_case' do
      mol = Rubabel["NCC(OO)CC"]
      frags = mol.fragment(rules: [:cod])
      frags.flatten(1).map {|a| a.map(&:csmiles)}.flatten.should == ["OC[NH3+]", "CCC=O", "C([NH3+])C=O", "CCO"]
    end

    specify 'rule: test_case' do
      mol = Rubabel["CCC(=O)O"]
      pieces = mol.fragment(rules: [:cod])
      pieces.flatten(1).map {|a| a.map(&:csmiles)}.flatten.should == ["[CH2-]C", "O=C=O"]
    end

    specify 'rule: test_case' do
      mol = Rubabel["CCC(=O)O"]
      mol.add_h!(1.5)
      pieces = mol.fragment(rules: [:cod])
      pieces.flatten(1).map {|a| a.map(&:csmiles)}.flatten.should == ["CC", "O=C=O"]
    end

  end
end
