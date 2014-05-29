require 'mspire/mzml'
require 'mspire/peaklist'
require 'mspire/spectrum'
require 'mspire/bin'
require 'optparse'
require 'pry'

PPM_BIN_DEFAULT = 5

options = {}
opt_parse = OptionParser.new do |opts|
  opts.on_tail("-h", "--help", "Show this message") do 
    puts opts
    exit
  end
  opts.banner = "Usage: #{__FILE__} truth_file.mzML file2.mzML"
  opts.separator "Takes a set of spectral files and kicks out comparison stats for them"
  opts.separator "Output: comparison.csv"

  opts.on("-l", "--list x,y,z", Array, "List of scan numbers to extract") do |scans|
    options[:selected_scans] = scans.map(&:to_i)
  end


end

opt_parse.parse!(ARGV)

if ARGV.size != 2
  puts opt_parse
  exit
end

files = ARGV

def ppm(m1,m2)
  (m2-m1)/m1*1e6
end
def ppm_range(mass, ppm = PPM_BIN_DEFAULT)
  range = (ppm*mass)/1e6
  (mass-range..mass+range)
end

def create_bins(true_values, comparison_values, bin_width: PPM_BIN_DEFAULT)
  min, max = [true_values.minmax, comparison_values.minmax].flatten.minmax
  divisions = []

  puts "using bin width: #{bin_width}" if $VERBOSE
  puts "using ppm for bins: #{use_ppm}" if $VERBOSE

  current_x = min
  loop do
    if current_x >= max
      divisions << max
      break
    else
      divisions << current_x
      current_x += current_x./(1e6).*(bin_width)
    end
  end
  # make each bin exclusive so there is no overlap
  bins = divisions.each_cons(2).map {|pair| Mspire::Bin.new(*pair, true) }
  # make the last bin *inclusive* of the terminating value
  bins[-1] = Mspire::Bin.new(bins.last.begin, bins.last.end)
  bins
end

def simple_matching(true_values, comparison_values, bin_width: PPM_BIN_DEFAULT)
  bins = create_bins(true_values, comparison_values, bin_width: bin_width)
  p bins.size
  puts "MIGHT WANT TO ADJUST YOUR BIN SIZE... this might take a while!  " if bins.size > 200_000
  bins2 = create_bins(true_values, comparison_values, bin_width: bin_width)
  p bins2.size
  Mspire::Bin.bin(bins, true_values)
  bins.reject! {|a| a.data.empty?}
  puts 'cloning...'
  match_bins = bins.clone
  puts 'cloning done'
  Mspire::Bin.bin(match_bins, comparison_values[0..1])
  puts "finishes with only 2 values"
  puts "won't finish with >= 3"
  Mspire::Bin.bin(match_bins, comparison_values)
  puts 'never gets here'
  Mspire::Bin.bin(bins2, comparison_values)
  match_bins.select! {|a| a.data.size > 1}
  bins2.reject! {|a| a.data.empty?}
  [bins, match_bins, bins2]
end

#testing
#  p simple_matching([350.0, 400.00, 450.0],[375.0], bin_width: 5)
#  p simple_matching([350.0, 400.00, 450.0],[350.0], bin_width: 5)

# Takes arrays of mz values
def f1_score(trues, matches, total_searched)
  # sensitivity/precision
  t_count = trues.size
  true_positives = matches.size
  total_searched_count = total_searched.size
  false_positives = total_searched_count - true_positives
  false_negatives = t_count - true_positives
  precision = true_positives/total_searched_count.to_f
  sensitivity = true_positives/t_count.to_f
  p "TP: #{true_positives}"
  p "FP: #{false_positives}"
  p "FN: #{false_negatives}"
  p "Precision: #{precision}"
  p "Sensitivity: #{sensitivity}"
  
end

p f1_score(Array.new(40),Array.new(3), Array.new(5))

bins = simple_matching([350.0, 400.00, 450.0],[350.0,399.0,355.0,424.0], bin_width: 500)
p bins.first.size
p f1_score(*bins)

def top_n_peak_optimization(truth, spectra, method)
  # Here, optimize on a particular method for the best score at N depth
  
end

def andromeda_algorithm(truth_values, spectra_values)
  # 
end
def spectral_contrast_angle(s1, s2)
  # This is from http://www.sciencedirect.com/science/article/pii/S1044030501003270
  # similarity index, and spectral angle...
  
end

