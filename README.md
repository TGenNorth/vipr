[![Go Report Card](https://goreportcard.com/badge/github.com/TGenNorth/vipr)](https://goreportcard.com/report/github.com/TGenNorth/vipr)
[![DOI](https://zenodo.org/badge/140470631.svg)](https://zenodo.org/badge/latestdoi/140470631)

# ViPR

Virtual PCR

```
# The following is a example of running neben in a directory of fastas
for assembly_fasta in *.fasta; do
  neben -min 200 -max 400 -primers <primer_list> [-probe <probe_sequence>] $assembly_fasta >> <output_file>
done
```
