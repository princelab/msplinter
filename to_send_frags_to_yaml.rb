require 'msplinter'
require 'yaml'
WRITE_YAML = false
WRITE_SVGS = false



lmids = ARGV
lmids << "LMSP02010009" if lmids.size < 1
lmids.each do |lmid|
  mol = Rubabel[lmid, :lmid]
  mol.adduct << "[Li+]"
  frags = mol.fragment add_rule_names: true
  # MAKE THE OUTDIR
  outdir = lmid
  require 'fileutils'
  FileUtils.mkdir outdir unless Dir.exists? outdir
  # TO YAML?
  if WRITE_YAML
    frags2 = {}
    # CREATE FRAGMENTS
    frags.map {|k,v| next if v.empty? ;frags2[k] = v.flatten.compact.map(&:to_s).zip(v.flatten.compact.map(&:mass_with_adduct)) }
    outfile_name = File.join(outdir, lmid+'.yml')
    File.open(outfile_name, 'w') do |io|
      YAML.dump(frags2, io)  
    end  
  end
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
  outfile_name = File.join(outdir, lmid+'.txt')
  File.open(outfile_name, 'w') {|io| io.puts mol.fragment.flatten.map(&:mass_with_adduct).flatten.uniq(&:to_i).sort }
end # lmids
