import gzip
import argparse

import subprocess
from subprocess import PIPE, Popen, STDOUT

from sys import stdin

from collections import defaultdict

import numpy as np

###############################################################################
#### DEFAULT VALUES

#KMER_LENGTH = 6

KNOWN_ADAPTER_RATIO_THRESHOLD = 0.4

SEQ_DICT    = { 0: "A", 1: "C", 2: "G", 3: "T" }

# Order of the known adapters is important
# They are searched in the given order.
KNOWN_ADAPTERS = [\
                     ## The most commonly used adapter sequence
                    "CTGTAGGCACCATCAAT",

                    ## This adapter has been mentioned in experiments and their GEO page :
                    ## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1660950
                    ## It is essentially Tru-Seq adapter ligated to a pre-Adenylated molecule
                    ## Hence we see the extra A on the left.
                    "AAGATCGGAAGAGCACACGTCT",

                    # TruSeq Adapter
                    "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",

                    # Truseq small RNA adapter
                    "TGGAATTCTCGGGTGCCAAGG",

                    # Poly A Tail
                    #"AAAAAAAAAA"
                 ]

#######################################################

###############################################################################
###############################################################################
### Data Structures to Read Fastq Files #######################################

class FastqEntry:
    """
    Helper function for Fastq File
    In each iteration, a FastqEntry is returned for each read
    """
    def __init__(self , header , sequence , plus , quality):
        self.header   = header
        self.sequence = sequence
        self.plus     = plus
        self.quality  = quality

    def reverse_comp_sequence(self):
        complement_table = {"A" : "T" , "C" : "G" , "G" : "C", "T" : "A", "N" : "N" ,
                            "a" : "T" , "c" : "G" , "g" : "C", "t" : "C", "n" : "N"}
        raw_result = self.sequence[::-1]
        result = list()
        for nucleotide in raw_result:
            result.append(complement_table[nucleotide])

        return ''.join(result)

    def reverse_complement(self):
        self.sequence = self.reverse_comp_sequence()
        self.quality  = self.quality[::-1]

    def __str__(self ):
        return("\n".join( ('@' + self.header , self.sequence , self.plus , self.quality ) ) )


class FastqFile:
    """
    Fastq File Reader Class
    """
    def __init__(self , file):
        myopen = open
        if file.endswith(".gz"):
            myopen = gzip.open

        if(file):
            self.f = myopen(file , "rt")
        else:
            #self.f = "a"
            self.f = stdin

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        pass

    def __del__(self):
        self.f.close()

    def __getitem__(self , index ):
        header   = self.f.readline().rstrip()
        if(header == ""):
            self.f.close()
            raise(IndexError)
        if header[0] != '@':
            raise(IOError('The first character of the fastq file must me a @'))
        sequence = self.f.readline().rstrip()
        plus     = self.f.readline().rstrip()
        quality  = self.f.readline().rstrip()
        return(FastqEntry( header[1:] , sequence , plus , quality  ))


###############################################################################
###############################################################################
###############################################################################
###########      F U N C T I O N S      #######################################


def convert_int_to_kmer(i, k):
    """converts integer i to a sequence of length k
    Example: 21 -> 000111 when k = 6  """
    digits = list()

    powers = list( reversed( range(0, k) ) )

    for x in powers:
        this_power = 4 ** x
        divisor    = i // this_power
        i          = i - ( divisor * this_power )

        digits.append( SEQ_DICT[divisor] )

    return "".join(digits)


########################################################

def generate_kmers(k):
    """Generates all DNA sequences of length k
    for k = 6, the result looks as follows:

    AAAAAA
    AAAAAC
    ......
    TTTTTT
    """

    kmers = list()

    for m in range(0, 4**k):
        kmers.append( convert_int_to_kmer(m, k) )

    return kmers


############################################################

############################################################

def get_reads(fastq_file, samplesize = 1000, skippedreads = 1000000):
    """
    Gets the nucleotide sequences from the FASTQ file
    """

    file_path        = fastq_file
    read_buffer      = []
    pre_read_max     = samplesize + skippedreads

    # First we check if there are sufficiently many reads in the file
    # if the number of reads is less than the sample size
    # the file is rejected.

    # If the file has less than skipped_reads + sample size reads
    # then we take the last "sample_size"

    # This is not memory efficient as we keep the reads in memory
    # But it is faster as we are reading the file only once
    # Otherwise we had to do two passes
    # First for counting reads and second for grabbing the reads.
    with FastqFile(file_path) as fastq_input:
        for entry in  fastq_input:

            # If we collected sufficient reads, we can stop reading
            if len(read_buffer) > pre_read_max:
                break

            read_buffer.append(entry.sequence)

    if len(read_buffer) >= samplesize:
         return read_buffer[len(read_buffer) - samplesize:]
    else:
        print("{} has insufficient reads ({})".format(fastq_file, len(read_buffer)))
        print("At least {} reads are required!".format(samplesize))
        exit(1)


    """
    fastq_reads  = list()
    read_counter = 0


    sequence_pointer = 0

    ### First read the sequences into an array
    ### For this, we skip the first SKIPPED_READS
    ### Then fill the array with the next samplesize many reads
    with FastqFile(file_path) as fastq_input:
        for entry in fastq_input:
            if read_counter < skippedreads:
                read_counter += 1
                continue

            fastq_reads.append(entry.sequence)

            sequence_pointer += 1
            if sequence_pointer >= samplesize:
                break

    return fastq_reads
    """

#############################################################

def count_all_kmers(fastq_reads, kmers):
    kmer_counters = { x: 0 for x in kmers  }
    # Count all kmers
    for this_read in fastq_reads:
        for k in kmers:
            search_result = this_read.find(k)
            if search_result >= 0:
                kmer_counters[k] += 1

    return kmer_counters
#############################################################


def determine_anchor_sequence(fastq_reads, kmer_counters,
                              minlength, kmerlength,
                              matchratio, verbose = False):
    """
    Determines the anchor sequence to be extended to the adapter
    """


    # seq length -> { sequence -> frequency  }
    top_hits    = dict()
    top_hits[kmerlength] = \
        dict( (sorted(kmer_counters.items(),
        key = lambda kv: (kv[1], kv[0]) ) )[-1*minlength:] )
    #print(top_hits[kmerlength])

    pre_adapter = sorted(top_hits[kmerlength].items(), key = lambda kv: (kv[1], kv[0] ) )[-1]

    # loop over the given adapter length range
    for i in range(kmerlength + 1, minlength + 1):
        tmp_counter = dict()

        # loop over the previous kmers
        for thismer in top_hits[i-1].keys():
            this_extension = [thismer + c for c in ["A", "C", "G", "T"] ]

            # loop over the extended kmer
            for e in this_extension:
                tmp_counter[e] = 0

                # search the kmer in the fastq reads
                for r in fastq_reads:
                     search_result = r.find(e)
                     if search_result >= 0:
                         tmp_counter[e] += 1

        top_hits[i] = \
            dict( sorted(tmp_counter.items(), key = lambda kv: (kv[1], kv[0]) ) [-1*minlength:])

        sorted_hits = sorted(top_hits[i].items(), key = lambda kv: (kv[1], kv[0] ) )
        this_ratio  = sorted_hits[-1][1] / len(fastq_reads)

        if this_ratio < matchratio:
            break

        pre_adapter = sorted_hits[-1]



    if verbose:
        print("Extending the sequence: ", pre_adapter)

    # Now we find the extensions of the pre_adapter

    pre_sequence = pre_adapter[0]
    return pre_sequence

##############################################################################

def get_trails(fastq_reads, pre_sequence):
    """
    Trails are the subsequences of the reads starting with the searched kmer
    """
    trails = []

    #print("get_trails:", pre_sequence)

    for this_read in fastq_reads:
        search_result = this_read.find(pre_sequence)

        if search_result < 0:
            continue

        trails.append( this_read[search_result:] )

    return trails

def count_extension_nucleotides(pre_sequence, trails, maxlength):
    ###
    # seq length -> { sequence -> frequency  }
    nuc_counter    = defaultdict(int)

    for i in range(len(pre_sequence), maxlength):
        nuc_counter[i] = defaultdict(int)
        for this_t in trails:
            #print(this_t)
            if len(this_t) > i:
                nuc_counter[i][this_t[i]] += 1

    return nuc_counter

##############################################################################

def _guess_adapter(fastq_reads,
                   minlength, maxlength,
                   kmerlength, matchratio,
                   verbose = False):

    kmers         = generate_kmers(kmerlength)

    kmer_counters = count_all_kmers(fastq_reads, kmers)

    pre_sequence  = determine_anchor_sequence(fastq_reads, kmer_counters, minlength, kmerlength, matchratio,
                                              verbose = verbose)

    # Each entry of trails holds a partial read that starts with the pre_sequence
    trails        = get_trails(fastq_reads, pre_sequence)

    ###
    # seq length -> { sequence -> frequency  }
    nuc_counter   = count_extension_nucleotides(pre_sequence, trails, maxlength)

    ## At this point "nuc_counter" holds the frequency of
    ## each nucleotide for a given length.
    ## So at each length we determine the most frequent nucleotide
    ## if its percentage is above a threshold, we keep extending.
    ## we stop if any of the following happen
    ##  a) The percentage is below the THRESHOLD
    ##  b) We reach maxlength

    if len(pre_sequence) < minlength:
        print("Failed to find the anchor sequence with length {}. ".format(minlength) +\
              "Found \"{}\"".format(pre_sequence) +\
              " length = {}".format(len(pre_sequence)))
        exit(1)

    #adapter_guess = pre_sequence[0]
    adapter_guess = pre_sequence

    if verbose:
        print("Pre Sequence is ", pre_sequence)
        print(nuc_counter)

    for i in range(len(pre_sequence), maxlength):
        sorted_nucs = sorted(nuc_counter[i].items(),
                              key = lambda kv: (kv[1], kv[0] ) )

        total_count = sum( list(map(lambda x: x[1], sorted_nucs))  )

        if len(sorted_nucs) == 0:
            break
        
        max_nuc   = sorted_nucs[-1][0]
        max_count = sorted_nucs[-1][1]

        this_ratio = max_count / total_count

        if this_ratio >= matchratio and len(adapter_guess) < maxlength:
            adapter_guess += max_nuc
        else:
            break

    return adapter_guess

################################################################################

def are_read_lens_different(reads):
    """
    We expect the reads to be of the same size.
    If their sizes are different, this may indicate that the adapters are pre-trimmed.
    """

    if len(reads) == 0:
        print(" Error! No reads given to check_read_len_variation")
        exit(1)

    min_len = len(reads[0])
    max_len = len(reads[0])

    for r in reads[1:]:
        if len(r) > max_len:
            max_len = len(r)
        if len(r) < min_len:
            min_len = len(r)

    if min_len == max_len:
        return False
    else:
        return( [min_len, max_len])

################################################################################

def get_parameters():
    parser = argparse.ArgumentParser(
                  description = "This script tries to guess the 3p adapter after sampling reads from a given FASTQ File. "
                                "The default values of the parameters have been optimized and should work for the majority of the cases. "
                                "The input file can be gzipped or plain. The format is guessed from the file extension.")

    parser.add_argument("fastq",
                        type = str)

    parser.add_argument('-v', '--verbose',
                        action = 'store_true')

    parser.add_argument("--samplesize",
                        type     = int,
                        default  = 10000,
                        required = False,
                        help     = "Number of reads to be sampled from the fastq file."
                                   "Default is 10000")

    parser.add_argument("--skippedreads",
                        type     = int,
                        default  = 1000000,
                        required = False,
                        help     = "Number of reads to be skipped in the fastq file."
                                   "Default is 1M, i.e., first 1M reads are skipped.")

    parser.add_argument("--minlength",
                        type     = int,
                        default  = 10,
                        required = False,
                        help     = "Minimum adapter length"
                                   " as long as match threshold is satisfied. "
                                   "Default is 10.")

    parser.add_argument("--maxlength",
                        type     = int,
                        default  = 20,
                        required = False,
                        help     = "Maximum adapter length. "
                                   "Default is 20.")

    parser.add_argument("--kmerlength",
                        type     = int,
                        default  = 6,
                        required = False,
                        help     = "Number of nucleotides to be matched to initiate the search. "
                                   "After kmer matches, they are extended to an anchor sequence of "
                                   "length <= minlength. "
                                   "Default is 6.")

    parser.add_argument("--matchratio",
                        type     = float,
                        default  = 0.5,
                        required = False,
                        help     = "Minimum required match ratio to extend the kmer or anchor sequence. "
                                   "Default is 0.5.")

    parser.add_argument("--skipnucs",
                        type     = int,
                        default  = 25,
                        required = False,
                        help     = "Number of nucleotides to be skipped from the 5' end of the reads. "
                                   "In other words, the adapter is  searched ONLY AFTER the first skipnucs nucleotides in the reads. "
                                   "Default value is 25.")

    parser.add_argument('--cutadapt',
                        action = 'store_true',
                        help   =  "Run cutadapt on  the sample using the guessed adapter. "
                                  "It will display the statistics and warnings of cutadapt.")

    args = parser.parse_args()
    return args


################################################################################

def _look_for_known_adapter(adapter,
                            fastq_reads, kmerlength,
                            matchratio,
                            skipnucs):

    this_kmer = adapter[:kmerlength]
    frequency = 0

    for r in fastq_reads:
        if r[skipnucs:].find(this_kmer) >= 0:
            frequency += 1

    return (adapter, frequency)


################################################################################


def look_for_known_adapters(fastq_reads, kmerlength, minlength,
                            matchratio,
                            skipnucs):
    """
    Searches for the known adapter sequences
    IN ALL KNOWN ADAPTERS

    It tries to anchor the first kmerlength sequences first.
    """

    # Each element is of the form (adapter, number_of_matches)
    lookup_results = []

    # First we determine the positions of the first kmer sequences of the known adapters
    for k_a in KNOWN_ADAPTERS:
        this_result = _look_for_known_adapter(
                                    adapter     = k_a,
                                    fastq_reads = fastq_reads,
                                    kmerlength  = kmerlength,
                                    matchratio  = matchratio,
                                    skipnucs    = skipnucs)

        lookup_results.append(this_result)

    selected_adapter = lookup_results[0][0]
    adapter_freq     = lookup_results[0][1]

    # Select the adapter with highest frequency
    for this_adapter, this_freq in lookup_results:
        if this_freq > adapter_freq:
            selected_adapter = this_adapter
            adapter_freq     = this_freq


    observed_ratio = (adapter_freq / len(fastq_reads) )

    # print("obs ration", observed_ratio)
    # If the frequency of the adapter is low,
    # return nothing
    """
    if observed_ratio > matchratio:
        return (selected_adapter, adapter_freq)
    else:
        return None
    """
    if observed_ratio < matchratio:
        return None


    ### Now extend the adapter to minlength and search again
    adapter_anchor = selected_adapter[:minlength]
    anchor_matches = 0

    for r in fastq_reads:
        if r.find(adapter_anchor) >= 0:
            anchor_matches += 1

    anchor_ratio = ( anchor_matches / len( fastq_reads ) )

    # print("anch ration", anchor_ratio)

    if anchor_ratio >= matchratio:
        return selected_adapter
    else:
        return None

################################################################################



def guess_adapter(fastq, verbose, samplesize, skippedreads,
                  minlength, maxlength,
                  kmerlength, matchratio,
                  skipnucs, cutadapt):
    """Guesses the adapter from the given FASTQ File"""
    fastq_reads = get_reads(fastq, samplesize, skippedreads)

    read_lens_different = are_read_lens_different(fastq_reads)
    if read_lens_different and ( read_lens_different[0] != read_lens_different[1]) :
        print("Reads are variable length!")
        print("Min length is {} and max length is {}".format(
                 read_lens_different[0], read_lens_different[1]) )
        return None

    if verbose:
        print("Read lnegth is {}".format(len(fastq_reads[0])))

    known_adapter_guess = look_for_known_adapters(
                                fastq_reads = fastq_reads,
                                kmerlength  = kmerlength,
                                minlength   = minlength,
                                matchratio  = KNOWN_ADAPTER_RATIO_THRESHOLD ,
                                skipnucs    = skipnucs)

    if known_adapter_guess:
        # print(known_adapter_guess)
        return known_adapter_guess

    if verbose:
        print("No known adapters were found.")

    if verbose:
        read_lens_different = are_read_lens_different(fastq_reads)

        if read_lens_different:
            print("Warning! Different read lengths ( {} and {}) were observed!".\
                   format(read_lens_different[0], read_lens_different[1]))

    return _guess_adapter(fastq_reads,
                          minlength  = minlength,
                          maxlength  = maxlength,
                          kmerlength = kmerlength,
                          matchratio = matchratio,
                          verbose    = verbose)

################################################################################

def run_cutadapt(fastq_file, adapter,
                 skipped_reads, samplesize ,
                 minimum_read_length = 15):

    """
    It runs cutadapt on the sampled reads used to guess the adapter.
    As the 3p adapter, it uses the guessed adapter.

    """

    head_lines = (4 * skipped_reads) + (4 * samplesize)
    tail_lines = 4 * samplesize

    mycommand = "zcat {fastq} | head -n {head_lines} | tail -n {tail_lines}".\
                  format(fastq      = fastq_file,
                         head_lines = head_lines,
                         tail_lines = tail_lines)

    mycommand += " | cutadapt -a {adapter} - 1>/dev/null ".format(adapter = adapter)

    print(mycommand)

    #ps = Popen(mycommand,shell=True,stdout=PIPE,stderr=STDOUT)
    ps = subprocess.run(mycommand,
                        shell  = True,
                        stdout = PIPE,
                        stderr = STDOUT)

    #output = ps.communicate()[0]
    #ps.stdout.close()
    output = ps.stdout.decode("utf-8")
    print(output)


################################################################################


def main():
    params  = get_parameters()

    adapter = guess_adapter( **vars(params) )

    print("\nAdapter Guess is:\n{}".format(adapter) )

    if params.cutadapt:
        run_cutadapt(fastq_file    = params.fastq,
                     skipped_reads = params.skippedreads,
                     samplesize    = params.samplesize,
                     adapter       = adapter)


##############################################################

if __name__ == "__main__":
    main()
