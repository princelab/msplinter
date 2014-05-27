
require 'spec_helper'

require 'rubabel'
require 'rubabel/molecule/fragmentable'
require 'pry'

$VERBOSE = nil

describe Rubabel do 
  describe 'Fragmentation of a adduct bond molecule' do 
    describe 'single adduct' do 
      it "carries the adduct into the fragments" do 
        m = Rubabel["LMGP02010003",:lmid]
        m.adducts = ["[Li+]"]
        frags = m.fragment
        frags.flatten.map(&:adducts).uniq.map {|a| a.first.csmiles}.should == ["[Li+]"]
      end
      it 'generates a proper mass for a fragment ion with an adduct' do 
        m = Rubabel["LMGP02010003",:lmid]
        m.adducts = ["[Li+]"]
        frags = m.fragment.flatten
        frags.first.csmiles.should == "CCCCCCCCCCCCCCCCC(=O)OC[CH+]COP(=O)(OCC[NH3+])[O-]"
        frags.first.mass.should be_close(450.29790119499995, 1e-6)
        frags.first.mass_with_adduct.first.should be_close(457.31335714499994, 1e-6)
      end
      it 'generates both the adduct-less mass and the adduct mass when requested' do 
        m = Rubabel["LMGP02010003",:lmid]
        m.adducts = ["[Li+]"]
        frags = m.fragment.flatten
        frags.first.csmiles.should == "CCCCCCCCCCCCCCCCC(=O)OC[CH+]COP(=O)(OCC[NH3+])[O-]"
        frags.first.mass.should be_close(450.29790119499995, 1e-6)
        frags.first.mass_with_adduct.first.should be_close(457.31335714499994, 1e-6)
        frags.first.masses_with_possible_adducts.should == [450.29790119499995, 457.31335714499994]
      end
    end
    describe 'multiple adducts' do 
      it "carries the adduct into the fragments" do 
        m = Rubabel["LMGP02010003",:lmid]
        m.adducts = ["[Li+]", "[Na+]"]
        frags = m.fragment
        frags.flatten.map(&:adducts).uniq.map {|a| a.first.csmiles}.should == ["[Li+]"]
      end
      it 'generates a proper mass for a fragment ion with an adduct' do 
        m = Rubabel["LMGP02010003",:lmid]
        m.adducts = ["[Li+]", "[Na+]"]
        frags = m.fragment.flatten
        frags.first.csmiles.should == "CCCCCCCCCCCCCCCCC(=O)OC[CH+]COP(=O)(OCC[NH3+])[O-]"
        frags.first.mass.should be_close(450.29790119499995, 1e-6)
        frags.first.mass_with_adduct.should == [457.31335714499994, 473.28712187599996]
      end
      it 'generates both the adduct-less mass and the adduct mass when requested' do 
        m = Rubabel["LMGP02010003",:lmid]
        m.adducts = ["[Li+]", "[Na+]"]
        frags = m.fragment.flatten
        frags.first.csmiles.should == "CCCCCCCCCCCCCCCCC(=O)OC[CH+]COP(=O)(OCC[NH3+])[O-]"
        frags.first.mass.should be_close(450.29790119499995, 1e-6)
        frags.first.mass_with_adduct.first.should be_close(457.31335714499994, 1e-6)
        frags.first.mass_with_adduct.last.should be_close(473.28712187599996, 1e-6)
        frags.first.masses_with_possible_adducts.should == [450.29790119499995, 457.31335714499994, 473.28712187599996]
      end
  end
  end
end

  
