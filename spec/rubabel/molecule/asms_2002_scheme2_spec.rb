require 'spec_helper'

require 'rubabel'
require 'rubabel/molecule/fragmentable'

$VERBOSE = nil

describe Rubabel::Molecule::Fragmentable do

  describe 'JASMS 2002 scheme 2' do 
    # Comment here
    specify 'rule: jasms_2002_2_f1' do
      mol = Rubabel["OCC(NC(=O)C(O)C)C(O)C=C"]
      wt = mol.mol_wt
      frags = mol.fragment(rules: [:jasms_2002_2_f1], rearrange: false)
      frags.map{|a| a.map(&:csmiles) }.flatten.should == ["CC=O", "O=CNC(C(C=C)O)CO"].reverse
      frags.map {|a| a.map(&:mol_wt)}.flatten.inject(:+).should be_close(wt,MASS_DELTA)
    end
    specify 'rule: jasms_2002_2_f1_aprime' do
      mol = Rubabel["O=CNC(C(C=C)O)CO"]
      wt = mol.mol_wt
      frags = mol.fragment(rules: [:jasms_2002_2_f1_aprime], rearrange: false)
      frags.map {|a| a.map(&:csmiles)}.flatten.should == ["O=CNC1COC1C=C", "O"]
      frags.map {|a| a.map(&:mol_wt)}.flatten.inject(:+).should be_close(wt,MASS_DELTA)
    end
    specify 'rule: jasms_2002_2_f1_bprime' do
      mol = Rubabel["O=CNC(C(C=C)O)CO"]
      wt = mol.mol_wt
      frags = mol.fragment(rules: [:jasms_2002_2_f1_bprime], rearrange: false)
      frags.map {|a| a.map(&:csmiles)}.flatten.should == ["C=CC(C1CN1C=O)O", "O"]
      frags.map {|a| a.map(&:mol_wt)}.flatten.inject(:+).should be_close(wt,MASS_DELTA)
    end
    specify 'rearrangement: jasms_2002_2_f2_b_methoxy_loss' do
      mol = Rubabel["O=CNC1COC1C=C"]
      wt = mol.mol_wt
      frags = mol.rearrange(rules: [:jasms_2002_2_f2_b_methoxy_loss])
      frags.map {|a| a.map(&:csmiles)}.flatten.should == ["[CH2-][OH]", "O=CNC=CC=C"].reverse
      frags.map {|a| a.map(&:mol_wt)}.flatten.inject(:+).should be_close(wt,MASS_DELTA)
    end
    specify 'rearrangement: jasms_2002_2_f2_bprime_water_loss' do
      mol = Rubabel["O=CN1CC1C(O)C=C"]
      wt = mol.mol_wt
      frags = mol.rearrange(rules: [:jasms_2002_2_f2_bprime_water_loss])
      frags.map {|a| a.map(&:csmiles)}.flatten.should == ["O", "C=CC=C1CN1C=O"].reverse
      frags.map {|a| a.map(&:mol_wt)}.flatten.inject(:+).should be_close(wt,MASS_DELTA)
    end
    specify 'rule: jasms_2002_2_e1' do
      mol = Rubabel["OCC(NC(=O)C(O)C)C(O)C=C"]
      wt = mol.mol_wt
      frags = mol.fragment(rules: [:jasms_2002_2_e1], rearrange: false)
      frags.map {|a| a.map(&:csmiles)}.flatten.should == ["CC=O", "[O-]#C", "OCC(C(C=C)O)N"].reverse # Should be a C#O without any charge, although I don't find spectral evidence for this ion
      frags.map {|a| a.map(&:mol_wt)}.flatten.inject(:+).should be_close(wt,MASS_DELTA)
    end
    specify 'rule: jasms_2002_2_e1_b' do 
      mol = Rubabel["OCC(C(C=C)O)N"]
      wt = mol.mol_wt
      frags = mol.fragment(rules: [:jasms_2002_2_e1_b], rearrange: false)
      frags.map {|a| a.map(&:csmiles)}.flatten.should == ["O", "NC1COC1C=C"].reverse
      frags.map {|a| a.map(&:mol_wt)}.flatten.inject(:+).should be_close(wt,MASS_DELTA)
    end
    specify 'rule: jasms_2002_2_e1_bprime' do 
      mol = Rubabel["OCC(C(C=C)O)N"]
      wt = mol.mol_wt
      frags = mol.fragment(rules: [:jasms_2002_2_e1_bprime], rearrange: false)
      frags.map {|a| a.map(&:csmiles)}.flatten.should == ["O", "OC(C1CN1)C=C"]
      frags.map {|a| a.map(&:mol_wt)}.flatten.inject(:+).should be_close(wt,MASS_DELTA)
    end
    specify 'rearrangement: jasms_2002_2_e2_b_formaldehyde_loss' do 
      mol = Rubabel["NC1COC1C=C"]
      wt = mol.mol_wt
      frags = mol.rearrange(rules: [:jasms_2002_2_e2_b_formaldehyde_loss])
      frags.map {|a| a.map(&:csmiles)}.flatten.should == ["C=O", "NC=CC=C"].reverse
      frags.map {|a| a.map(&:mol_wt)}.flatten.inject(:+).should be_close(wt,MASS_DELTA)
    end
    specify 'rearrangement: jasms_2002_2_e2_bprime_water_loss' do 
      mol = Rubabel["N1CC1C(O)C=C"]
      wt = mol.mol_wt
      frags = mol.rearrange(rules: [:jasms_2002_2_e2_bprime_water_loss])
      frags.map {|a| a.map(&:csmiles)}.flatten.should == ["O", "C=CC=C1CN1"].reverse
      frags.map {|a| a.map(&:mol_wt)}.flatten.inject(:+).should be_close(wt,MASS_DELTA)
    end
    specify 'rule: jasms_2002_2_e2_bprime_heterocyclic_loss' do 
      mol = Rubabel["N1CC1C(O)C=C"]
      wt = mol.mol_wt
      frags = mol.fragment(rules: [:jasms_2002_2_e2_bprime_heterocyclic_loss], rearrange: false)
      frags.map {|a| a.map(&:csmiles)}.flatten.should == ["N1CC1", "C=CC=O"]
      frags.map {|a| a.map(&:mol_wt)}.flatten.inject(:+).should be_close(wt,MASS_DELTA)
    end
    specify 'rule: jasms_2002_2_c2a' do 
      mol = Rubabel["C1CN1C(=O)C(O)C"]
      wt = mol.mol_wt
      frags = mol.fragment(rules: [:jasms_2002_2_c2a], rearrange: false)
      frags.map {|a| a.map(&:csmiles)}.flatten.should == ["CC=O","C=O", "N1CC1"].reverse
      frags.map {|a| a.map(&:mol_wt)}.flatten.inject(:+).should be_close(wt,MASS_DELTA)
    end
    specify 'rule: jasms_2002_2_a1' do 
      mol = Rubabel["CC(O)C(=O)NC1COC1C=C"]
      wt = mol.mol_wt
      frags = mol.fragment(rules: [:jasms_2002_2_a1], rearrange: false)
      frags.map {|a| a.map(&:csmiles)}.flatten.should == ["CC=O", "N=C=O","C=CC1CCO1"]
      frags.map {|a| a.map(&:mol_wt)}.flatten.inject(:+).should be_close(wt,MASS_DELTA)
    end
  end
end
