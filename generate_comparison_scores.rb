require 'mspire/mzml'
require 'mspire/peaklist'
require 'mspire/spectrum'
require 'mspire/bin'
require 'optparse'
require 'pry'

PPM_BIN_DEFAULT = 10
MAXQ = 100
VERBOSE = false

options = {ppm: PPM_BIN_DEFAULT, qmax: MAXQ}
opt_parse = OptionParser.new do |opts|
  opts.on_tail("-h", "--help", "Show this message") do 
    puts opts
    exit
  end
  opts.banner = "Usage: #{__FILE__} truth_file.mzML comparison_file.mzML ... truth_file[n].mzML comparison_file[n].mzML"
  opts.separator "Takes sets of spectral files and kicks out comparison stats for them"
  opts.separator "Output: truth_file_comparison_file.csv ... truth_file[n]_comparison_file[n].csv"

  opts.on("-v", "--verbose", "Verbose mode") do |v|
    VERBOSE = v
  end
  opts.on("-p", "--ppm N", Float, "Set the PPM tolerance for searches") do |p|
    options[:ppm] = p
  end
  opts.on("-q", "--qmax N", Integer, "Set the maximum Q value for searches") do |q|
    options[:qmax] = q 
  end
end

opt_parse.parse!(ARGV)

if ARGV.size % 2 != 0 or ARGV.size == 0
  puts opt_parse
  exit
end


def putsv(thing)
  puts thing if VERBOSE
end
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
  [true_values, t_e_bins, comparison_values]
end
def find_top_peaks_in_window(spectrum, window_bottom, window_top, q_max)
  putsv "window_bottom: #{window_bottom}"
  putsv "window_top: #{window_top}"
  bottom = spectrum.find_nearest_index(window_bottom)
  top = spectrum.find_nearest_index(window_top)
  putsv "window: #{spectrum[bottom..top]}"
  output = spectrum[bottom..top].transpose.sort_by(&:last).reverse[0..q_max]
  range = window_bottom..window_top
  output.select {|a| range.include?(a.first) }
end
def divide_masses_into_windows(masses)
  min, max = masses.minmax
  # return window start values to divide the spectrum into ~evenly sized windows < 100 Th
  range = max-min
  bins = (range / 100.0).ceil
  size = range / bins.to_f
  resp = (bins).times.map {|i| min+i*size}
  resp << max
end

def matching(gold_standard_spectrum, tested_spectrum, ppm: PPM_BIN_DEFAULT, q_max: MAXQ)
  # Only need masses
  tested_masses = tested_spectrum.mzs
  gold_standard_masses = gold_standard_spectrum.mzs
  # Filter each and every 100 m/z region into top Q peaks (optimized by score)
  windows_for_q = divide_masses_into_windows(tested_masses)
  top_peak_lists = windows_for_q.each_cons(2).map {|b, t| find_top_peaks_in_window(tested_spectrum, b, t, q_max)}
  q_lists = (1..q_max).to_a.map do |q|
    # return a list of the top q peaks per region
    top_peak_lists.map{|arr| arr[0...q]}.flatten(1)
  end
  # Matched by ppm or AMU tolerance 
  matches = q_lists.map.with_index do |list, i|
    scores = list.map do |peak|
      # finds closest match in the GS spectrum
      closest_match = gold_standard_spectrum[gold_standard_spectrum.find_nearest_index(peak.first)]

      #putsv "Potential problem with #{closest_match}: bigger peak nearby" if larger_peaks_nearby(gold_standard_spectrum, closest_match)
      ppm_between(peak.first,closest_match.first) < ppm ? peak.first : nil
    end.compact
  end
  matches
end

#testing
#  p simple_matching([350.0, 400.00, 450.0],[375.0], bin_width: 5)
#  p simple_matching([350.0, 400.00, 450.0],[350.0], bin_width: 5)
#  esp = Mspire::Spectrum.new [[200.0198,200.01998, 200.0199998,205.35],[123,2123,2134,22351]]
#  tsp = Mspire::Spectrum.new [[200.02,205.35,232.34,268.245,270.24],Array.new(5,1000)]
#  p matching(tsp, esp, q_max: 5)

# Takes arrays of mz values
def calculate_f1_score(trues, matches, tests)
  # sensitivity/precision
  false_positives = tests - trues
  false_negatives = trues - matches
  precision = matches/(matches + false_positives).to_f || 0
  sensitivity = matches/trues.to_f
  ans = 2*(precision*sensitivity)/(precision+sensitivity) || 0
  ans = 0 if ans.nan?
  ans
end

def stats_from_matches(trues, matches, attempts)
  [trues,matches,attempts].map(&:size)
end

def optimize_f1_score(gss, cs, ppm: PPM_BIN_DEFAULT, q_max: MAXQ)
  matches = matching(gss, cs, ppm: ppm, q_max: q_max)
  scores = (1..q_max).to_a.map.with_index do |q,i|
    [calculate_f1_score(*stats_from_matches(gss.mzs, matches[i], cs.mzs))*q_max/q.to_f, q]
  end
  scores.sort_by {|a| a.first}.reverse
end

## ANDROMEDA
require_relative 'andromeda'
  # This gives the class AndromedaPscore
 # call the function #optimize to grab the score and q for the best match possible
def spectral_contrast_angle(s1, s2)
  # This is from http://www.sciencedirect.com/science/article/pii/S1044030501003270
  # similarity index, and spectral angle...
end

## HIT COUNT = X*hits-Y*misses+X*unpredicted  where X is a scaling factor 10>=X>=1 and Y is a scaling factor 2>X>=1
# hits 
def hit_count_from_numbers(hits, misses, unpredicted, x: 10, y: 1, z: 1)
  x*hits - misses*y - z*unpredicted
end
def calculate_hit_count(matches, predicted, q)
  match_count = matches.size
  total_checked = predicted.size
  hit_count_from_numbers(match_count, total_checked-match_count, q-total_checked)
end
def optimize_hit_count(gss,cs, q_max: MAXQ, ppm: PPM_BIN_DEFAULT)
  matches = matching(gss,cs, ppm: ppm, q_max: q_max)
  scores = (1..q_max).to_a.map.with_index do |q,i|
    [calculate_hit_count(matches[i], cs.mzs, q), q]
  end
  scores.sort_by(&:first).reverse
end

# RUN ANALYSIS
ARGV.each_slice(2) do |ground_truth_file, comparison_file|
  gss = Mspire::Mzml.open(ground_truth_file) {|mzml| mzml[0]}
  cs = Mspire::Mzml.open(comparison_file) {|m| m[0]}

  outfile = [ground_truth_file, comparison_file].map {|f| File.basename(f).gsub(File.extname(f),'')}.join('_') + '.csv'
  p outfile
  File.open(outfile, 'w') do |io|
    output = []
    output << %w(truth_file comparison_file optimizedF1-score @qdepth andromedaPscore @q_depth hitcountscore @q_h_depth).join(",")
    gss = Mspire::Mzml.open(ground_truth_file) {|mzml| mzml[0]}
    cs = Mspire::Mzml.open(comparison_file) {|m| m[0]}
    f1 = calculate_f1_score(*stats_from_matches(*simple_matching(gss.mzs, cs.mzs, bin_width: options[:ppm])))
    scores = AndromedaPscore.optimize(gss, cs, ppm: options[:ppm], q_max: options[:qmax])
    optimized_f1, q_f1 = optimize_f1_score(gss, cs, ppm: options[:ppm], q_max: options[:qmax]).first
    a_score, a_q = scores.first
    hitcount, h_q = optimize_hit_count(gss, cs, ppm: options[:ppm], q_max: options[:qmax]).first
    puts "F1: #{f1}"
    puts "optimized_F1: #{optimized_f1} @ q=#{q_f1}"
    puts "Andromeda: #{a_score} @ q=#{a_q}"
    puts "HitCount: #{hitcount} @ q=#{h_q}"
    output << [ground_truth_file, comparison_file, optimized_f1, q_f1, a_score, h_q, hitcount, h_q].join(",")
    io.puts output
  end
end



# Testing stuff
if false #$0 == __FILE__
  p calculate_f1_score(Array.new(40),Array.new(3), Array.new(5))
  bins = simple_matching([350.0, 400.00, 450.0],[350.0,400.0,355.0,424.0], bin_width: 500)
  p bins
  p bins.first.size
  p calculate_f1_score(*bins)

  esp = Mspire::Spectrum.new [[200.0199,431.86,205.35,700],[123,21023,2134,22351]]
  tsp = Mspire::Spectrum.new [[200.02,205.35,232.34,268.245,270.24],Array.new(5,1000)]
  p AndromedaPscore.optimize(tsp,esp)
  p calculate_f1_score(*simple_matching(esp.mzs, tsp.mzs))
end



