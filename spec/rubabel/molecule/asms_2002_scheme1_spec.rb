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
      second_frags = frags.flatten.map {|a| a.rearrange_but_return_sets(rules: [:nitrogen_cyclization]) }
      second_frags.first[:nitrogen_cyclization].flatten.map(&:csmiles).uniq.should == ["O", "CC(=O)N1CC1"]
    end

    specify 'rule: b_d1' do
      mol = Rubabel["OCC(NC(=O)C)C(O)C=C"]
      frags = mol.fragment(rules: [:asms_2002_scheme1_b_d1], rearrange: false)
      resp = frags.map {|a| a.map(&:csmiles)}
      frags.map {|a| a.map(&:csmiles)}.flatten.should == ["OC(C1CO1)C=C", "CC(=O)N", "OCC1OC1C=C","CC(=O)N"]
    end

    specify 'rearrangement: b_d1prime_water_loss' do
      mol = Rubabel["OCC(NC(=O)C)C(O)C=C"]
      frags = mol.fragment(rules: [:asms_2002_scheme1_b_d1], rearrange: false)
      second_frags = frags.flatten.map {|a| a.rearrange(rules: [:asms_2002_scheme1_b_d1_water_loss])}
      second_frags.flatten.map(&:csmiles).uniq.should == [ "CC#N","O", "C=CC=C1CO1"].reverse
    end

    specify 'rearrangement: b_d1prime_formaldehyde_loss' do
      mol = Rubabel["OCC(NC(=O)C)C(O)C=C"]
      frags = mol.fragment(rules: [:asms_2002_scheme1_b_d1], rearrange: false)
      second_frags = frags.flatten.map {|a| a.rearrange(rules: [:asms_2002_scheme1_b_d1_formaldehyde_loss] ) }
      second_frags.flatten.map(&:csmiles).uniq.should == ["C=CC1CO1", "C=O"].reverse
    end
    specify "rule: c_e1" do 
      mol = Rubabel["OCC(NC(=O)CCCl)C(O)C=C"]
      frags = mol.fragment(rules: [:asms_2002_scheme1_c_e1], rearrange: false)
      frags.map {|a| a.map(&:csmiles) }.flatten.should == ["ClCC=C=O", "OCC(C(C=C)O)N"].reverse
    end
    specify "rule: c_e1aprime" do 
      mol = Rubabel["OCC(C(C=C)O)N"]
      frags = mol.fragment(rules: [:asms_2002_scheme1_c_e1aprime], rearrange: false)
      frags.map {|a| a.map(&:csmiles) }.flatten.should == ["O", "NC1COC1C=C"].reverse
    end
    specify "rule: c_e1bprime" do 
      mol = Rubabel["OCC(C(C=C)O)N"]
      frags = mol.fragment(rules: [:asms_2002_scheme1_c_e1bprime], rearrange: false)
      frags.map {|a| a.map(&:csmiles) }.flatten.should == ["O", "OC(C1CN1)C=C"]
    end
    specify "rule: c_e1b_to_d1b" do 
      mol = Rubabel["OCC(C(C=C)O)N"]
      frags = mol.fragment(rules: [:asms_2002_scheme1_c_e1b_to_d1], rearrange: false)
      frags.map {|a| a.map(&:csmiles) }.flatten.should == ["O", "OC(C1CN1)C=C"]
    end
    specify "rule: c_e1b_to_d1bprime" do 
      mol = Rubabel["OCC(C(C=C)O)N"]
      frags = mol.fragment(rules: [:asms_2002_scheme1_c_e1b_to_d1bprime], rearrange: false)
      frags.map {|a| a.map(&:csmiles) }.flatten.should == ["O", "OC(C1CN1)C=C"]
    end
  end
end

# mol.write 'test.svg'
#frags.flatten.map.with_index {|m,i| m.write "test_#{i}.svg" }
