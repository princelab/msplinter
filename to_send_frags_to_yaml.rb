require 'msplinter'
require 'yaml'
lmid = ARGV.first
lmid ||= "LMSP02010009"
mol = Rubabel[lmid, :lmid]
mol.adduct << "[Li+]"
frags = mol.fragment add_rule_names: true
# MAKE THE OUTDIR
outdir = lmid
require 'fileutils'
FileUtils.mkdir outdir unless Dir.exists? outdir
# TO YAML?
WRITE_YAML = true
if WRITE_YAML
  frags2 = {}
# CREATE FRAGMENTS
  frags.map {|k,v| next if v.empty? ;frags2[k] = v.flatten.compact.map(&:to_s).zip(v.flatten.compact.map(&:mass_with_adduct)) }
  outfile_name = File.join(outdir, lmid+'.yml')
  File.open(outfile_name, 'w') do |io|
    YAML.dump(frags2, io)  
  end  
end
WRITE_SVGS = false
if WRITE_SVGS
  frags.map do |k,v| 
    v.each_with_index do |mol, i|
      outfile_name = File.join(outdir,k.to_s + "_#{mol.mass_with_adduct.join("_")}.svg")
      p outfile_name
      mol.write_file(outfile_name)
    end
  end
end


puts "HERE ARE THE PREDICTED MASSES:"
puts mol.fragment.flatten.map(&:mass_with_adduct).flatten.uniq(&:to_i).sort
