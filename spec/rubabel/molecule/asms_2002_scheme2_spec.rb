require 'spec_helper'

require 'rubabel'
require 'rubabel/molecule/fragmentable'

$VERBOSE = nil

describe Rubabel::Molecule::Fragmentable do

  describe 'JASMS 2002 scheme 2' do 
    # Comment here
    specify 'rule: jasms_2002_2_f1' do
      mol = Rubabel["OCC(NC(=O)C(O)C)C(O)C=C"]
      frags = mol.fragment(rules: [:jasms_2002_2_f1], rearrange: false)
      frags.map{|a| a.map(&:csmiles) }.flatten.should == ["CC=O", "O=CNC(C(C=C)O)CO"].reverse
      frags.map {|a| a.map(&:mol_wt)}.flatten.inject(:+).should be_close(mol.mol_wt,MASS_DELTA)
    end
    specify 'rule: jasms_2002_2_f1_aprime' do
      mol = Rubabel["O=CNC(C(C=C)O)CO"]
      frags = mol.fragment(rules: [:jasms_2002_2_f1_aprime], rearrange: false)
      frags.map {|a| a.map(&:csmiles)}.flatten.should == ["O=CNC1COC1C=C", "O"]
      frags.map {|a| a.map(&:mol_wt)}.flatten.inject(:+).should be_close(mol.mol_wt,MASS_DELTA)
    end
    specify 'rule: jasms_2002_2_f1_bprime' do
      mol = Rubabel["O=CNC(C(C=C)O)CO"]
      frags = mol.fragment(rules: [:jasms_2002_2_f1_bprime], rearrange: false)
      frags.map {|a| a.map(&:csmiles)}.flatten.should == ["C=CC(C1CN1C=O)O", "O"]
      frags.map {|a| a.map(&:mol_wt)}.flatten.inject(:+).should be_close(mol.mol_wt,MASS_DELTA)
    end
    specify 'rearrangement: jasms_2002_2_f2_b_methanol_loss' do
      mol = Rubabel["O=CNC1COC1C=C"]
      frags = mol.rearrange(rules: [:jasms_2002_2_f2_b_methanol_loss])
      frags.map {|a| a.map(&:csmiles)}.flatten.should == ["CO", "O=CNC=CC=C"].reverse
      frags.map {|a| a.map(&:mol_wt)}.flatten.inject(:+).should be_close(mol.mol_wt,MASS_DELTA)
    end
    specify 'rearrangement: jasms_2002_2_f2_bprime_water_loss' do
      mol = Rubabel["O=CN1CC1C(O)C=C"]
      frags = mol.rearrange(rules: [:jasms_2002_2_f2_bprime_water_loss])
      frags.map {|a| a.map(&:csmiles)}.flatten.should == ["O", "C=CC=C1CN1C=O"].reverse
      frags.map {|a| a.map(&:mol_wt)}.flatten.inject(:+).should be_close(mol.mol_wt,MASS_DELTA)
    end
    specify 'rule: jasms_2002_2_e1' do
      mol = Rubabel["OCC(NC(=O)C(O)C)C(O)C=C"]
      frags = mol.fragment(rules: [:jasms_2002_2_e1], rearrange: false)
      frags.map {|a| a.map(&:csmiles)}.flatten.should == ["CC=O", "[O-]#C", "OCC(C(C=C)O)N"].reverse # Should be a C#O without any charge
      frags.map {|a| a.map(&:mol_wt)}.flatten.inject(:+).should be_close(mol.mol_wt,MASS_DELTA)
    end
  end
end
