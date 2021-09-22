# flu_fasta_check
Flu fasta files checking for flu database. 

Inf_seg_check.py accepts 8 filenames as input, each of which contains sequences from one of the eight Influenza A/B virus segments. The script checks the validity of fasta files and match sequences to a corresponding segment.

```
usage: Inf_seg_check.py [-h] [--max_N_ratio MAX_N_RATIO] [--min_chunk_length MIN_CHUNK_LENGTH] [--max_segment_length MAX_SEGMENT_LENGTH] input [input ...] {A,B}

positional arguments:
  input                 Path to 8 fasta files to analyze. Each file corresponds to the Influenza virus segment: segment1.fasta, segment2.fasta...segment8.fasta
  {A,B}                 Type of Influenza virus (A or B).

optional arguments:
  -h, --help            show this help message and exit
  --max_N_ratio MAX_N_RATIO
                        Maximum permissible content of N in each sequence , default: 0.5
  --min_chunk_length MIN_CHUNK_LENGTH
                        minimum length for each chunk of segment, default: 50
  --max_segment_length MAX_SEGMENT_LENGTH
                        maximum length of segment, default: 3000


```

