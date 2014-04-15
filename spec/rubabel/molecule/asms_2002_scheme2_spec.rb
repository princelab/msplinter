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

    specify 'rule: test_case' do
      mol = Rubabel["OC(NC(=O)C(O)C)C(O)C=C"]
      frags = mol.fragment(rules: [:word], rearrange: false)
      frags.map {|a| a.map(&:csmiles)}.flatten.should == ["[CH2-]C", "O=C=O"]
      frags.map {|a| a.map(&:mol_wt)}.flatten.inject(:+).should be_close(mol.mol_wt,MASS_DELTA)

    end

    specify 'rule: test_case' do
      mol = Rubabel["OC(NC(=O)C(O)C)C(O)C=C"]
      frags = mol.fragment(rules: [:word], rearrange: false)
      frags.map {|a| a.map(&:csmiles)}.flatten.should == ["CC", "O=C=O"]
    end

  end
end
