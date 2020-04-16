import itertools
import time
import math
import xalglib
# import the .NET System namespace so we can create .NET arrays from Python
import System

# Defines some helper functions for the algorithm

def get_trace(i, j):
    """The one-dimensional array of values for the cube's trace"""
    return cube.column(i, j).get_rawvalues()

def manhattan_coords(ci, cj, d):
    """A generator of i,j tuples which have a manhattan distance less than 
    or equal to d from (ci, cj), not including (ci, cj).  E.g. 

            3
           323
          32123
         321X123
          32123
           323
            3
    """
    for i in range(ci-d, ci+d+1):
        for j in range(cj-d, cj+d+1):
            if i == ci and j == cj:
                continue
            if abs(i - ci) + abs(j - cj) <= d:
                yield (i, j)
        
def hzs_to_samples(freqs_upper, f_sample):
    """Returns the (lower, upper) sample indexes for the frequency bins, from 0Hz, given a particular sample frequency"""
    s_lower = 0
    for f_u in freqs_upper:
        s_upper = math.floor(f_u/f_sample)
        yield (int(s_lower), int(s_upper))
        s_lower = s_upper + 1


def reorder(xs, l):
    """Used only by correlate"""
    if l % 2 == 0:
        raise ArgumentError("Length must be odd")
    return list(reversed(xs[0:l/2 + 1])) + list(reversed(xs[-l/2 + 1:]))

def correlate(xs, ys):
    """The cross-correlation of ys with xs, with output length the same as the input.
    Pass the same values twice to give auto-correlation"""
    if len(xs) != len(ys):
        raise ArgumentError("Lengths must be the same")
    l = len(xs)
    vs = xalglib.corrr1d(list(ys), l, list(xs), l)
    rs = reorder(vs, l)
    return rs

def smooth(xs, smoothlen):
    """Box filter values"""
    ones = [1.0/smoothlen] * smoothlen
    l = len(xs)
    if l % 2 == 0:
        raise ArgumentError("Length must be odd")
    return xalglib.convr1d(xs, l, ones, smoothlen)[smoothlen/2:-smoothlen/2+1]

# Helper functions which perform arithmetic on sequences

def add(a, b):
    """Adds two sequences, returns an array"""
    if len(a) != len(b):
        raise ArgumentError("lengths do not match")
    s = System.Array.CreateInstance(System.Double, len(a))
    for i in range(0, len(a)):
        s[i] = a[i] + b[i]
    return s


def div(a, b):
    """Divides a sequence by a scalar, returns an array"""
    s = System.Array.CreateInstance(System.Double, len(a))
    for i in range(0, len(a)):
        s[i] = a[i] / b
    return s

def diva(a, b):
    """Divides two sequences element-wise.  If b has zeros, the result element is 0 """
    s = System.Array.CreateInstance(System.Double, len(a))
    for i in range(0, len(a)):
        if b[i] != 0:
            s[i] = a[i] / b[i]
        else:
            s[i] = 0.0
    return s

def sub(a, b):
    """Subtracts two sequences, returns an array"""
    s = System.Array.CreateInstance(System.Double, len(a))
    for i in range(0, len(a)):
        s[i] = a[i] - b[i]
    return s


# Perform the SNR calculation in a particular region. Avoid edges as 
# the algorithm assumes all surrounding traces are available
min_i = 5
max_i = 345
min_j = 5
max_j = 395

sample_rate_hz = 250  # Taken by inspection of the input cube
num_samples = len(get_trace(min_i, min_j))
sample_width = num_samples/sample_rate_hz

# Output three histogram maps, [0-10Hz], [10-40Hz, 40-120Hz]
bins_hz = [10, 40, 120]

bins_samples = list(hzs_to_samples(bins_hz, sample_width))

# Store the values for the result maps provided in a dict
# so they can be easily set as the calculation proceeds
bin0_values = bin0.all().get_rawvalues()
bin1_values = bin1.all().get_rawvalues()
bin2_values = bin2.all().get_rawvalues()

results_dict = {bins_samples[0] : bin0_values,
                bins_samples[1] : bin1_values,
                bins_samples[2] : bin2_values}



def get_mean(ijs):
    """Returns the mean of the traces given by the sequence of (i, j) coordinate tuples in ijs"""
    t_sum = [0] * num_samples
    for i, j in ijs:
        t_sum = add(t_sum, get_trace(i, j))
    return div(t_sum, len(ijs))

# Book-keeping for tracking progress
num_traces = 0 
total_traces = (max_i - min_i) * (max_j - min_j)
start = time.time()

# Iterate over all the traces in the specified region
for t_i in range(min_i, max_i):
    for t_j in range(min_j, max_j):
        trace = get_trace(t_i, t_j)
        # Auto-correlation
        t_ac = correlate(trace, trace) 
    
        # Cross-correlation with the surrounding traces (manhattan distance of 3)
        t_mean = get_mean(list(manhattan_coords(t_i, t_j, 3)))
        t_xc = correlate(trace, t_mean)

        # The signal+noise power spectrum (from the auto-correlated trace).  Assumes data is already 0-phase, so we throw away the imaginary component
        snp = [abs(c.real) for c in xalglib.fftr1d(t_ac)]
        snp_smooth = smooth(snp, 9)

        # The signal power spectrum (from the cross-correlation with surrounding traces).  Assumes data is 0-phase.
        sp = [abs(c.real) for c in xalglib.fftr1d(t_xc)]
        sp_smooth= smooth(sp, 9)

        # SNR = 1 - sp/(sp + snp)  
        snr = sub([1] * num_samples, diva(sp_smooth, add(sp_smooth, snp_smooth)))

        # Divide the snr samples into the bins specified at the top, and assign the mean snr for the bin to the appropriate result map 
        for bin in bins_samples:
            snr_slice = snr[bin[0]:bin[1]+1]
            snr_for_bin = xalglib.samplemean(list(snr_slice))
            results_dict[bin][t_i, t_j] = snr_for_bin
        
        num_traces = num_traces + 1
        if num_traces % 500 == 0:
            print("%d of %d at %d tps" % (num_traces, total_traces, num_traces/ (time.time() - start)))
        
print("Finished calculating.  Setting result map values...")
bin0.all().set(bin0_values)
bin1.all().set(bin1_values)
bin2.all().set(bin2_values)
print("Done.")


