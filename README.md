[![Go Report Card](https://goreportcard.com/badge/github.com/TGenNorth/vipr)](https://goreportcard.com/report/github.com/TGenNorth/vipr)
[![DOI](https://zenodo.org/badge/140470631.svg)](https://zenodo.org/badge/latestdoi/140470631)

# ViPR

Virtual PCR

ViPR is intented to function using a list of forward and reverse primers, a max allowed mismatch count, and a list of sequences from a fasta file. The script creates a list of the full possible primer mismatches plus any possiblities from partial call symbols such as K being G or T. This full list of possible versions of forward and reverse primers are then searched for using a regex search on each sequence. ViPR mainly looks to count pairs of forward and reverse primers which both appear in the same sequence. The ViPR script can then output either a list of what it found in each file and sequence, or a file which shows a matrix table of which forwards hit with each of reverse primers. This script is also designed with multithreading for handling the large possiblity space of all the primers and should be used inside a job handler like Slurm.

```
# The following is a example of running vipr in a directory of fastas
for assembly_fasta in *.fasta; do
  vipr -min 200 -max 400 -primers <primer_list> [-probe <probe_sequence>] $assembly_fasta >> <output_file>
done
```
