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

  opts.on("-v", "--verbose", "Verbose mode") do |v|
    VERBOSE = v
  end
end

opt_parse.parse!(ARGV)

if ARGV.size != 2
  puts opt_parse
  exit
end

ground_truth_file, comparison_file = ARGV


def ppm(m1,m2)
  (m2-m1)/m1*1e6
end
def ppm_range(mass, ppm = PPM_BIN_DEFAULT)
  range = (ppm*mass)/1e6
  (mass-range..mass+range)
end
def ppm_between(m1, m2)
  ((m2-m1)/m1*1e6).abs
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
  puts "MIGHT WANT TO ADJUST YOUR BIN SIZE... this might take a while!  " if bins.size > 2_000_000
  matches_bins = create_bins(true_values, comparison_values, bin_width: bin_width)
  # Truth and Experiment bins
  t_e_bins = create_bins(true_values, comparison_values, bin_width: bin_width)
  Mspire::Bin.bin(bins, true_values)
  bins.reject! {|a| a.data.empty?}
  Mspire::Bin.bin(t_e_bins, true_values)
  Mspire::Bin.bin(t_e_bins, comparison_values)
  t_e_bins.select! {|a| a.data.size > 1}
  Mspire::Bin.bin(matches_bins, comparison_values)
  matches_bins.reject! {|a| a.data.empty?}
  [bins, t_e_bins, comparison_values]
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
  2*(precision*sensitivity)/(precision+sensitivity)
end


## ANDROMEDA
require_relative 'andromeda'
  # This gives the class AndromedaPscore
 # call the function #optimize to grab the score and q for the best match possible
def spectral_contrast_angle(s1, s2)
  # This is from http://www.sciencedirect.com/science/article/pii/S1044030501003270
  # similarity index, and spectral angle...
  
end

# RUN ANALYSIS
outfile = ARGV.map {|f| File.basename(f).gsub(File.extname(f),'')}.join('_') + '.csv'
File.open(outfile, 'w') do |io|
  output = []
  output << %w(truth_file comparison_file F1-score andromedaPscore @q_depth).join(",")
  gss = Mspire::Mzml.open(ground_truth_file) {|mzml| mzml[0]}
  cs = Mspire::Mzml.open(comparison_file) {|m| m[0]}
  f1 = f1_score(*simple_matching(gss.mzs, cs.mzs))
  score, q = AndromedaPscore.optimize(gss, cs)
  output << [ground_truth_file, comparison_file, f1, score, q].join(",")
  io.puts output
end




# Testing stuff
if $0 == __FILE__
  p f1_score(Array.new(40),Array.new(3), Array.new(5))
  bins = simple_matching([350.0, 400.00, 450.0],[350.0,400.0,355.0,424.0], bin_width: 500)
  p bins
  p bins.first.size
  p f1_score(*bins)

  esp = Mspire::Spectrum.new [[200.0199,431.86,205.35,700],[123,21023,2134,22351]]
  tsp = Mspire::Spectrum.new [[200.02,205.35,232.34,268.245,270.24],Array.new(5,1000)]
  p AndromedaPscore.optimize(tsp,esp)
  p f1_score(*simple_matching(esp.mzs, tsp.mzs))
end



