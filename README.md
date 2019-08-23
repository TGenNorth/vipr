[![Go Report Card](https://goreportcard.com/badge/github.com/TGenNorth/vipr)](https://goreportcard.com/report/github.com/TGenNorth/vipr)

# ViPR

Virtual PCR

```
# The following is a example of running neben in a directory of fastas
for assembly_fasta in *.fasta; do
  neben -min 200 -max 400 -primers <primer_list> [-probe <probe_sequence>] $assembly_fasta >> <output_file>
done
```
