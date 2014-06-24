require 'spec_helper'

require 'rubabel'
require 'rubabel/molecule/fragmentable'

$VERBOSE = nil
MASS_DELTA = 10/1e6

describe Rubabel::Molecule::Fragmentable do

  describe 'JASMS 2002 scheme 1' do 
    # Comment here
    specify 'rule: a_c1' do
      mol = Rubabel["OCC(NC(=O)C)C(O)C=C"]
      mol.adducts << "[Li+]"
      frags = mol.fragment(rules: [:jasms_2002_scheme1_a_c1], rearrange: false)
      frags.map{|a| a.map(&:csmiles) }.flatten.should == ["C=CC=O", "CC(=O)NCCO"].reverse
      second_frags = frags.flatten.map {|a| a.rearrange_but_return_sets(rules: [:nitrogen_cyclization]) }
      second_frags.first[:nitrogen_cyclization].flatten.map(&:csmiles).uniq.should == ["O", "CC(=O)N1CC1"]
    end

    specify 'rule: b_d1' do
      mol = Rubabel["OCC(NC(=O)C)C(O)C=C"]
      mol.adducts << "[Li+]"
      frags = mol.fragment(rules: [:jasms_2002_scheme1_b_d1], rearrange: false)
      resp = frags.map {|a| a.map(&:csmiles)}
      frags.map {|a| a.map(&:csmiles)}.flatten.should == ["OC(C1CO1)C=C", "CC(=O)N", "OCC1OC1C=C","CC(=O)N"]
   #   frags.map {|a| a.map(&:mass)}.flatten.inject(:+).should be_close(mol.mass, MASS_DELTA)
    end

    specify 'rearrangement: b_d1prime_water_loss' do
      mol = Rubabel["OCC(NC(=O)C)C(O)C=C"]
      mol.adducts << "[Li+]"
      frags = mol.fragment(rules: [:jasms_2002_scheme1_b_d1], rearrange: false)
      second_frags = frags.flatten.map {|a| a.rearrange(rules: [:jasms_2002_scheme1_b_d1_water_loss])}
      second_frags.flatten.map(&:csmiles).uniq.should == [ "CC#N","O", "C=CC=C1CO1"].reverse
    end

    specify 'rearrangement: b_d1prime_formaldehyde_loss' do
      mol = Rubabel["OCC(NC(=O)C)C(O)C=C"]
      mol.adducts << "[Li+]"
      frags = mol.fragment(rules: [:jasms_2002_scheme1_b_d1], rearrange: false)
      second_frags = frags.flatten.map {|a| a.rearrange(rules: [:jasms_2002_scheme1_b_d1_formaldehyde_loss] ) }
      second_frags.flatten.map(&:csmiles).uniq.should == ["C=CC1CO1", "C=O"].reverse
    end
    specify "rule: c_e1" do 
      mol = Rubabel["OCC(NC(=O)CCCl)C(O)C=C"]
      mol.adducts << "[Li+]"
      frags = mol.fragment(rules: [:jasms_2002_scheme1_c_e1], rearrange: false)
      frags.map {|a| a.map(&:csmiles) }.flatten.should == ["ClCC=C=O", "OCC(C(C=C)O)N"].reverse
      frags.map {|a| a.map(&:mass)}.flatten.inject(:+).should be_close(mol.mass, MASS_DELTA)
    end
    specify "rule: c_e1aprime" do 
      mol = Rubabel["OCC(C(C=C)O)N"]
      mol.adducts << "[Li+]"
      frags = mol.fragment(rules: [:jasms_2002_scheme1_c_e1aprime], rearrange: false)
      frags.map {|a| a.map(&:csmiles) }.flatten.should == ["O", "NC1COC1C=C"].reverse
      frags.map {|a| a.map(&:mass)}.flatten.inject(:+).should be_close(mol.mass, MASS_DELTA)
    end
    specify "rule: c_e1bprime" do 
      mol = Rubabel["OCC(C(C=C)O)N"]
      mol.adducts << "[Li+]"
      frags = mol.fragment(rules: [:jasms_2002_scheme1_c_e1bprime], rearrange: false)
      frags.map {|a| a.map(&:csmiles) }.flatten.should == ["O", "OC(C1CN1)C=C"]
      frags.map {|a| a.map(&:mass)}.flatten.inject(:+).should be_close(mol.mass, MASS_DELTA)
    end
    specify "rule: c_e1b_to_d1b" do 
      mol = Rubabel["OCC(C(C=C)O)N"]
      mol.adducts << "[Li+]"
      frags = mol.fragment(rules: [:jasms_2002_scheme1_c_e1b_to_d1b], rearrange: false)
      frags.map {|a| a.map(&:csmiles) }.flatten.uniq.should == ["OC(C1CO1)C=C","N", "OCC1OC1C=C"]
  #   frags.map {|a| a.map(&:mass)}.flatten.inject(:+).should be_close(mol.mass, MASS_DELTA)
    end
    specify "rule: c_e1b_to_d1bprime" do 
      mol = Rubabel["OCC(C(C=C)O)N"]
      mol.adducts << "[Li+]"
      frags = mol.fragment(rules: [:jasms_2002_scheme1_c_e1b_to_d1bprime], rearrange: false)
      frags.map {|a| a.map(&:csmiles) }.flatten.uniq.should == ["OC(C1CO1)C=C","N", "OCC1OC1C=C"]
  #   frags.map {|a| a.map(&:mass)}.flatten.inject(:+).should be_close(mol.mass, MASS_DELTA)
    end
    specify "rule: c_e2aprime" do 
      mol = Rubabel["OCC(C(C=C)O)N"]
      mol.adducts << "[Li+]"
      frags = mol.fragment(rules: [:jasms_2002_scheme1_c_e2aprime], rearrange: false)
      # check that the adduct is also gone
      frags.flatten.map(&:adducts).map(&:size).inject(:+).should == 0
      frags.map {|a| a.map(&:csmiles) }.flatten.uniq.should == ["O", "C=CC1OCC1[NH3+]"]
      frags.map {|a| a.map(&:mass)}.flatten.inject(:+).should be_close(mol.mass, MASS_DELTA)
      
    end
    specify "rule: c_e2bprime" do 
      mol = Rubabel["OCC(C(C=C)O)N"]
      mol.adducts << "[Li+]"
      frags = mol.fragment(rules: [:jasms_2002_scheme1_c_e2bprime], rearrange: false)
      frags.flatten[1].write "e2bprime.svg"
      frags.map {|a| a.map(&:csmiles) }.flatten.uniq.should == ["OC(C1C[NH2+]1)C=C", "O"].reverse
      frags.map {|a| a.map(&:mass)}.flatten.inject(:+).should be_close(mol.mass, MASS_DELTA)
      #TODO check that the adduct is also gone
      frags.flatten.map(&:adducts).map(&:size).inject(:+).should == 0
    end
    specify "rule: c_e2aprime_formaldehyde_loss" do 
      mol = Rubabel["C=CC1OCC1[NH3+]"]
      mol.adducts << "[Li+]"
      frags = mol.fragment(rules: [:jasms_2002_scheme1_c_e2aprime_formaldehyde_loss], rearrange: false)
      # TODO, get this to give the right ions... pesky hydrogens!
      frags.map {|a| a.map(&:csmiles) }.flatten.uniq.should == ["[NH3+]C=CC=C", "CO"].reverse
      frags.map {|a| a.map(&:mass)}.flatten.inject(:+).should be_close(mol.mass, MASS_DELTA)
    end
    specify "rule: c_e2bprime_water_loss" do 
      mol = Rubabel["OC(C1C[NH2+]1)C=C"]
      mol.adducts << "[Li+]"
      frags = mol.fragment(rules: [:jasms_2002_scheme1_c_e2bprime_water_loss], rearrange: false)
      frags.map {|a| a.map(&:csmiles) }.flatten.uniq.should == ["O", "[NH2+]1CC1=CC=C"].reverse
      frags.map {|a| a.map(&:mass)}.flatten.inject(:+).should be_close(mol.mass, MASS_DELTA)
    end
    specify "rule: c_e1bbprimeprime" do 
      pending "TODO but not tonight"
      mol = Rubabel["NC1COC1C=C"]
      mol.adducts << "[Li+]"
      frags = mol.fragment(rules: [:jasms_2002_scheme1_c_e1bbprimeprime], rearrange: false)

    end
    specify "rearrangement: c_e1_aprimeprime_formaldehyde_loss" do 
      mol = Rubabel["NC1COC1C=C"]
      mol.adducts << "[Li+]"
      second_frags = mol.rearrange(rules: [:jasms_2002_scheme1_c_e1_aprimeprime_formaldehyde_loss] )
      second_frags.flatten.map(&:csmiles).uniq.should == ["NC=CC=C", "C=O"]
    end
    specify "rearrangement: c_e1_abprimeprime_water_loss" do 
      mol = Rubabel["OC(C1CN1)C=C"]
      second_frags = mol.rearrange(rules: [:jasms_2002_scheme1_c_e1_abprimeprime_water_loss] )
      second_frags.flatten.map(&:csmiles).uniq.should == ["O", "C=CC=C1CN1"]
      second_frags.map {|a| a.map(&:mass)}.flatten.inject(:+).should be_close(mol.mass, MASS_DELTA)
    end
  end
end

# mol.write 'test.svg'
#frags.flatten.map.with_index {|m,i| m.write "test_#{i}.svg" }
