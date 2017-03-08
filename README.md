# isPCR
in silico PCR

[![Build Status](https://circleci.com/gh/TGenNorth/isPCR/tree/master.svg?style=shield&circle-token=897111c87d78438dffb2bb5924437a42e4fc3a11)](https://circleci.com/gh/TGenNorth/isPCR)

```
# The following is a example of running neben in a directory of fastas
for assembly_fasta in *.fasta; do
  neben -min 200 -max 400 -primers <primer_list> [-probe <probe_sequence>] $assembly_fasta >> <output_file>
done
```
