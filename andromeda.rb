require 'mspire/mzml'
require 'mspire/spectrum'

VERBOSE = false
class Integer
  def factorial
    (1..self).inject(:*) || 1
  end
end


def putsv(thing)
  puts thing if VERBOSE
end

PPM_TOLERANCE = 10
MAXQ = 10
class AndromedaPscore
  def self.ppm_between(m1, m2)
    ((m2-m1)/m1*1e6).abs
  end
  def self.larger_peaks_nearby(spectrum, target_peak, nearby_window_amu: 50)
    target_mass, target_intensity = target_peak
    bottom = spectrum.find_nearest_index(target_mass - nearby_window_amu)
    top = spectrum.find_nearest_index(target_mass + nearby_window_amu)
    spectrum[bottom..top].transpose.select {|a| a.last > target_intensity }
  end
  def self.divide_masses_into_windows(masses)
    min, max = masses.minmax
    # return window start values to divide the spectrum into ~evenly sized windows < 100 Th
    range = max-min
    bins = (range / 100.0).ceil
    size = range / bins.to_f
    resp = (bins).times.map {|i| min+i*size}
    resp << max
  end
  def self.find_top_peaks_in_window(spectrum, window_bottom, window_top, q_max)
    putsv "window_bottom: #{window_bottom}"
    putsv "window_top: #{window_top}"
    bottom = spectrum.find_nearest_index(window_bottom)
    top = spectrum.find_nearest_index(window_top)
    putsv "window: #{spectrum[bottom..top]}"
    spectrum[bottom..top].transpose.sort_by(&:last)[0..q_max]
  end
  def self.calculate_ks(tested_spectrum, gold_standard_spectrum, q_max: MAXQ )
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
    q_scores = q_lists.map.with_index do |list, i|
      scores = list.map do |peak|
        closest_match = gold_standard_spectrum[gold_standard_spectrum.find_nearest_index(peak.first)]
        
        #putsv "Potential problem with #{closest_match}: bigger peak nearby" if larger_peaks_nearby(gold_standard_spectrum, closest_match)
        ppm_between(peak.first,closest_match.first) < PPM_TOLERANCE ? closest_match.first : nil
      end.compact
    end
    q_scores
  end
  def self.calculate_score(tested_spectrum, gold_standard_spectrum, k, q)
    n = gold_standard_spectrum.mzs.size
    permutations = (n.factorial / (k.factorial * (n-k).factorial).to_f)
    sum = (k..n).to_a.map {|j| permutations*(q/100.0)**j * (1-q/100.0)**(n-j) }.inject(:+)
    s = -10*Math::log10(sum)  # sum from j=k to n of (n j)(q/100)^j(1-q/100)^(n-j)
  end
  def self.optimize(tested_spectrum, gold_standard_spectrum, q_max: MAXQ)
    ks = calculate_ks(tested_spectrum, gold_standard_spectrum).map(&:size)
    scores = (1..q_max).to_a.map.with_index do |q,i|
      [calculate_score(tested_spectrum, gold_standard_spectrum, ks[i], q), q].flatten
    end
    scores.sort_by {|a| a.last}.first
  end
end

if __FILE__ == $0
  require 'pry'
  a = AndromedaPscore.new
  r = a.divide_masses_into_windows [100, 699]
  #p r
  #p r.each_cons(2).map {|a,b| b-a}
  esp = Mspire::Spectrum.new [[200.0199,431.86,205.35,700],[123,21023,2134,22351]]
  tsp = Mspire::Spectrum.new [[200.02,205.35,232.34,268.245,270.24],Array.new(5,1000)]
  #p a.larger_peaks_nearby(esp, esp.first)
  #p a.calculate_ks(tsp, esp)
  #p a.calculate_score(tsp, esp, 2, 1)
  score, q = a.optimize(tsp,esp)
  p score
  p q
  
end
