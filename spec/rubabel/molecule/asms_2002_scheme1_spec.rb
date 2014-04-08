require 'spec_helper'

require 'rubabel'
require 'rubabel/molecule/fragmentable'

$VERBOSE = nil

describe Rubabel::Molecule::Fragmentable do

  describe 'JASMS 2002 scheme 1' do 
    # Comment here
    specify 'rule: a_c1' do
      mol = Rubabel["OCC(NC(=O)C)C(O)C=C"]
      frags = mol.fragment(rules: [:asms_2002_scheme1_a_c1], rearrange: false)
      frags.map{|a| a.map(&:csmiles) }.flatten.should == ["C=CC=O", "CC(=O)NCCO"].reverse
      second_frags = frags.flatten.map {|a| a.rearrange }
      second_frags.first[:nitrogen_cyclization].map(&:csmiles).uniq.should == ["O", "CC(=O)N1CC1"]
    end

    specify 'rule: b_d1' do
      mol = Rubabel["OCC(NC(=O)C)C(O)C=C"]
      frags = mol.fragment(rules: [:asms_2002_scheme1_b_d1], rearrange: false)
      resp = frags.map {|a| a.map(&:csmiles)}
      frags.map {|a| a.map(&:csmiles)}.flatten.should == ["OC(C1CO1)C=C", "CC(=O)N", "OCC1OC1C=C","CC(=O)N"]
    end

    specify 'rule: b_d1prime' do
      mol = Rubabel["OCC(NC(=O)C)C(O)C=C"]
      pieces = mol.fragment(rules: [:asms_2002_scheme1_b_d1prime])
      pieces.flatten(1).map {|a| a.map(&:csmiles)}.flatten.should == []

    end

    specify 'rearrangement: b_d1prime_loss_formaldehyde' do
      pending
      mol = Rubabel["CCC(=O)O"]
      pieces = mol.fragment(rules: [:cod])
      pieces.flatten(1).map {|a| a.map(&:csmiles)}.flatten.should == ["CC", "O=C=O"]
    end

  end
end

# mol.write 'test.svg'
#frags.flatten.first.write "test.svg"
