require 'spec_helper'

require 'rubabel'
require 'rubabel/molecule/fragmentable'

$VERBOSE = nil

describe Rubabel::Molecule::Fragmentable do

  describe 'JASMS 2002 scheme 3' do 
    specify 'rule: jasms_2002_3_b1' do
      mol = Rubabel["OCC(NC(=O)CC)C(O)C(O)CC"]
      wt = mol.mol_wt
      frags = mol.fragment(rules: [:jasms_2002_3_b1], rearrange: false)
      frags.map{|a| a.map(&:csmiles) }.flatten.should == ["O", "CCC(C(C1CN1C(=O)CC)O)O"]
      frags.map {|a| a.map(&:mol_wt)}.flatten.inject(:+).should be_close(wt,MASS_DELTA)
    end
    specify 'rule: jasms_2002_3_b2' do
      mol = Rubabel["OCC(NC(=O)CC)C(O)C(O)CC"]
      wt = mol.mol_wt
      frags = mol.fragment(rules: [:jasms_2002_3_b2], rearrange: false)
      frags.map{|a| a.map(&:csmiles) }.flatten.should == ["CCC=O", "OCC(NC(=O)CC)C=O"].reverse # SHOULD ALSO HAVE H2
      frags.map {|a| a.map(&:mol_wt)}.flatten.inject(:+).should be_close(wt-(Rubabel["[H]"].mol_wt*2),MASS_DELTA)
    end
    specify 'rule: jasms_2002_3_b2_water_loss' do
      mol = Rubabel["OCC(NC(=O)CC)C=O"]
      wt = mol.mol_wt
      frags = mol.fragment(rules: [:jasms_2002_3_b2_water_loss], rearrange: false)
      frags.map{|a| a.map(&:csmiles) }.flatten.should == ["O", "CCC(=O)N1CC1C=O"]
      frags.map {|a| a.map(&:mol_wt)}.flatten.inject(:+).should be_close(wt,MASS_DELTA)
    end
    specify 'rule: jasms_2002_3_b2_formaldehyde_loss' do
      mol = Rubabel["OCC(NC(=O)CC)C=O"]
      wt = mol.mol_wt
      frags = mol.fragment(rules: [:jasms_2002_3_b2_formaldehyde_loss], rearrange: false)
      frags.map{|a| a.map(&:csmiles) }.flatten.should == ["C=O", "CCC(=O)NCC=O"]
      frags.map {|a| a.map(&:mol_wt)}.flatten.inject(:+).should be_close(wt,MASS_DELTA)
    end
  end

  describe 'JASMS 2002 scheme 4' do 
    specify 'rule: jasms_2002_4_g1' do
      mol = Rubabel["OCC(NC(=O)C)C(O)C=CC(O)C"]
      wt = mol.mol_wt
      frags = mol.fragment(rules: [:jasms_2002_4_g1], rearrange: false)
      frags.map{|a| a.map(&:csmiles) }.flatten.should == ["CC=O", "OCC(C(C=C)O)NC(=O)C"].reverse
      frags.map {|a| a.map(&:mol_wt)}.flatten.inject(:+).should be_close(wt,MASS_DELTA)
    end
    specify 'rule: jasms_2002_4_g2a' do
      mol = Rubabel["OCC(C(C=C)O)NC(=O)C"]
      wt = mol.mol_wt
      frags = mol.fragment(rules: [:jasms_2002_4_g2a], rearrange: false)
      frags.map{|a| a.map(&:csmiles) }.flatten.should == ["O", "CC(=O)NC1COC1C=C"]
      frags.map {|a| a.map(&:mol_wt)}.flatten.inject(:+).should be_close(wt,MASS_DELTA)
    end
    specify 'rule: jasms_2002_4_g2a_prime' do
      mol = Rubabel["OCC(C(C=C)O)NC(=O)C"]
      wt = mol.mol_wt
      frags = mol.fragment(rules: [:jasms_2002_4_g2a_prime], rearrange: false)
      frags.map{|a| a.map(&:csmiles) }.flatten.should == ["O", "CC(=O)N1CC1C(C=C)O"]
      frags.map {|a| a.map(&:mol_wt)}.flatten.inject(:+).should be_close(wt,MASS_DELTA)
    end
    specify 'rule: jasms_2002_4_g2a_formaldehyde_loss' do
      mol = Rubabel["CC(=O)NC1COC1C=C"]
      wt = mol.mol_wt
      frags = mol.fragment(rules: [:jasms_2002_4_g2a_formaldehyde_loss], rearrange: false)
      frags.map{|a| a.map(&:csmiles) }.flatten.should == ["C=O", "CC(=O)NC=CC=C"].reverse
      frags.map {|a| a.map(&:mol_wt)}.flatten.inject(:+).should be_close(wt,MASS_DELTA)
    end
    specify 'rule: jasms_2002_4_g2a_prime' do
      mol = Rubabel["CC(=O)N1CC1C(C=C)O"]
      wt = mol.mol_wt
      frags = mol.fragment(rules: [:jasms_2002_4_g2a_prime_water_loss], rearrange: false)
      frags.map{|a| a.map(&:csmiles) }.flatten.should == ["O", "CC(=O)N1CC1=CC=C"].reverse
      frags.map {|a| a.map(&:mol_wt)}.flatten.inject(:+).should be_close(wt,MASS_DELTA)
    end
  end
  describe 'JASMS 2002 scheme 5' do 
    specify 'rule: jasms_2002_5_e1_water_loss' do
      mol = Rubabel["CC(=O)NC(CO)C(O)C=CC(O)C"]
      wt = mol.mol_wt
      frags = mol.fragment(rules: [:jasms_2002_5_e1_water_loss], rearrange: false)
      frags.map{|a| a.map(&:csmiles) }.flatten.should == ["O", "CC(=O)NC1COC1C=CC(O)C"].reverse
      frags.map {|a| a.map(&:mol_wt)}.flatten.inject(:+).should be_close(wt,MASS_DELTA)
    end
    specify 'rule: jasms_2002_5_e2' do
      mol = Rubabel["CC(=O)NC1COC1C=CC(O)C"]
      wt = mol.mol_wt
      frags = mol.fragment(rules: [:jasms_2002_5_e2], rearrange: false)
      frags.map{|a| a.map(&:csmiles) }.flatten.should == ["C=C=O", "CC(C=CC1OCC1N)O"]
      frags.map {|a| a.map(&:mol_wt)}.flatten.inject(:+).should be_close(wt,MASS_DELTA)
    end
    specify 'rule: jasms_2002_5_e2_water_loss' do
      mol = Rubabel["CC(C=CC1OCC1N)O"]
      wt = mol.mol_wt
      frags = mol.fragment(rules: [:jasms_2002_5_e2_water_loss], rearrange: false)
      frags.map{|a| a.map(&:csmiles) }.flatten.should == ["O", "NC1COC1C=CC=C"].reverse
      frags.map {|a| a.map(&:mol_wt)}.flatten.inject(:+).should be_close(wt,MASS_DELTA)
    end
    specify 'rule: jasms_2002_5_e5_formaldehyde_loss' do
      mol = Rubabel["NC1COC1C=CC=C"]
      wt = mol.mol_wt
      frags = mol.fragment(rules: [:jasms_2002_5_e5_formaldehyde_loss], rearrange: false)
      frags.map{|a| a.map(&:csmiles) }.flatten.should == ["C=O", "NC=CC=CC=C"].reverse
      frags.map {|a| a.map(&:mol_wt)}.flatten.inject(:+).should be_close(wt,MASS_DELTA)
    end
    specify 'rule: jasms_2002_5_e2_adduct_loss' do
      pending "requires consideration of adduct loss, but actually isn't very complex" # NOTE that this step is actually after the water loss #rearrangement
      mol = Rubabel["CC(C=CC1OCC1N)O"]
      wt = mol.mol_wt
      frags = mol.fragment(rules: [:jasms_2002_5_e2_adduct_loss], rearrange: false)
      frags.map{|a| a.map(&:csmiles) }.flatten.should == ["LiO", "CC(C=CC1OCC1N)O"]
      frags.map {|a| a.map(&:mol_wt)}.flatten.inject(:+).should be_close(wt,MASS_DELTA)
    end
  end
end
