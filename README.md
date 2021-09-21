# flu_fasta_check
Flu fasta files checking for flu database. 

Inf_seg_check.py checks the confirmity of the sequence to the segment of Influenza A/B virus. Thie console runningg is for demonstration. The presented functions are supposed to be called from other script. 
```
usage: python Inf_seg_check.py [-h] [-seg {1,2,3,4,5,6,7,8} [{1,2,3,4,5,6,7,8} ...]]
                        input [input ...] {A,B} > output.csv

The script contains the set of functions needed to check if a sequence is the
certain segment of the Influenza A/B genome. Use it with parameter -help for
more detailed info.

positional arguments:
  input                 Path to fasta file to analyze. In the example all
                        sequences from the file will be checked.
  {A,B}                 Type of Influenza virus (A or B).

optional arguments:
  -h, --help            show this help message and exit
  -seg {1,2,3,4,5,6,7,8} [{1,2,3,4,5,6,7,8} ...]
                        Number of the segement(s) to check the sequence for.
                        PB2: 1 PB1 and PB1-F2: 2 PA: 3 HA: 4 NP: 5 NA: 6 M2
                        and M1: 7 NS2 and NS1: 8

```

