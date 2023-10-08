#main_functions

This file contains 3 main functions used for analyze biological data.

* run_dna_rna_tools - use for analysis DNA data. Details are described below 
* process_seqs - use for analysis Protein data. Details are described below 
* read_fastq - use for analysis DNA data from fastq file. Details are described below 

BioinfoTool contains 3 modules:

* DNATool.py
* DNA_fastq_filter.py
* ProteinTool.py

# DNATool

The DNATool contains 6 functions, which we can use to analyze out DNA data. 

##Functions:

transcribe — return the transcribed sequence
reverse — return the reversed sequence
complement — return the complementary sequence
reverse_complement — return the reverse complementary sequence
is_dna - used as an internal function to determine the type of sequence (DNA or RNA)

## DNATool using examples

Execute script (you should be on directory with script):
```bash
python3
>>> from main_functions import run_dna_rna_tools
>>>print(run_dna_rna_tools(__command__, __sequence or list of sequences__))
```

run_dna_rna_tools('transcribe', 'ATG') # 'AUG'
run_dna_rna_tools('reverse', 'ATG') # 'GTA'
run_dna_rna_tools('complement', 'AtG') # 'TaC'
run_dna_rna_tools('reverse_complement', 'ATg') # 'cAT'
run_dna_rna_tools( 'reverse', 'ATG', 'aT') # ['GTA', 'Ta']

# ProteinTool

## Tool for PROtein SEQuences Operation

This tool can perform some simple operations on amino acid sequences:
* help you calculate protein lengths, molecular weights, isoelectric points and GRAVY values
* find and show you heaviest and lightest proteins
* rewrite 1-letter sequence to 3-letter sequence

## How use process_seqs()
Execute script (you should be on directory with script):
```bash
python3
>>> from ProtSeqO import process_seqs
>>>print(process_seqs(__command__, __sequence or list of sequences__))
```

You can input to `process_seqs()` sequence as string or list with any strings of sequences. __Pay attention__ that your sequence(s) should contain 1-letter symbols (case does not matters) of 20 common amino acids ('U' for selenocysteine and 'O' for pyrrolysine doesn't allowed).

Command must be a string with one of followed options.

## process_seqs options
* 'lengths' - return list with numbers of AA in each sequence(s)
* 'molw' - return list of protein molecular weight (use the average molecular weight of AA, 110 Da)
* 'iso' - return list of approximate isoelectric point of given amino acids sequence
* 'gravy' - return list of GRAVY (grand average of hydropathy) values
* 'rename' - return list of sequences in 3-letter AA code (AA separated by hyphens)
* 'heavy' - return the sequence(s) with maximum molecular weight and weigth value
* 'light' - return the sequence(s) with minimum molecular weight and weigth value

# DNA_fastq_filter

DNA_fastq_filter contains 4 functions that are used to assess the quality of data and filter according to specified criteria.

##Functions:

cal_gc_content — return percent of GC content
gc_bounds — define filter boundatry and return dna if it passed the filtering using cal_gc_content
length_bounds — efine filter boundatry and return dna if it passed the filtering by set range
quality_threshold — Define quality of each nucleotide using ASCII table and Phred-33 scale


## DNATool using examples

Execute script (you should be on directory with script):
```bash
python3
>>> from main_functions import read_fastq
>>>print(read_fastq(__dic__, __gc_bounds__, __length_bounds, __quality__))
```

read_fastq(dict_name, (20, 80), (0, 2**32), 0) 
read_fastq(dict_name, (75), (1000), 0) 
