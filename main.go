package main

import (
	"bufio"
	"bytes"
	"encoding/json"
	"errors"
	"flag"
	"fmt"
	"index/suffixarray"
	"path"
	"sort"
	"strings"
	"unicode"
	"io"
	"log"
	"os"
	"sync"
	"time"
	"regexp"
	"io/ioutil"
)

/*
A Adenine
C Cytosine
G Guanine
T Thymine
U Uracil
W Weak A/T
S Strong C/G
M aMino A/C
K Keto G/T
R puRine A/G
Y pYrimidine C/T
B not A
D not C
H not G
V not T
N any Nucleotide
NOTE: . or - are valid characters to indicate a gap, but are excluded here as they are not relevant as primers.
*/
const IupacNucleotideCode = "ACGTUMRWSYKVHDBN"

const (
	MaxUint = ^uint(0)
	MinUint = 0
	MaxInt  = int(MaxUint >> 1)
	MinInt  = -MaxInt - 1
)

var (
	//debugFlag             bool
	maxMismatchFlag       int
	contigWorkersFlag     uint
	indexWorkersFlag      uint
	searchWorkersFlag     uint
	maxSequenceLengthFlag int
	minSequenceLengthFlag int
	primersFlag           PrimerList
	primersForwardFastaFlag       string
	primersReverseFastaFlag       string
	probeFlag             ProbeFlag
	profileFlag           string
	typeFlag              string
	regexJsonFileFlag     string
	seqToRegexMap         map[string][]string
)

func init() {
	// Disable concurrency flags to minimize arguments
	//flag.UintVar(&contigWorkersFlag, "concurrent-contig", 5, "max pending contigs")
	//flag.UintVar(&indexWorkersFlag, "concurrent-index", 5, "concurrent contig index builders")
	//flag.UintVar(&searchWorkersFlag, "concurrent-search", 10, "concurrent contig index searches")
	contigWorkersFlag = 5
	indexWorkersFlag = 5
	searchWorkersFlag = 10

	// TODO: do not allow negative
	flag.IntVar(&minSequenceLengthFlag, "min", 0, "minimum sequence length")
	flag.IntVar(&maxSequenceLengthFlag, "max", MaxInt, "maximum sequence length")
	flag.IntVar(&maxMismatchFlag, "max-mismatch", 0, "Maximum mismatches allowed")
	flag.Var(&primersFlag, "primers", "`PrimerList` is a filename or comma delimited list of forward followed by reverse primers to locate in the source contigs")
	flag.StringVar(&primersForwardFastaFlag, "forward-fasta-primers", "", "A filename for forward primers")
	flag.StringVar(&primersReverseFastaFlag, "reverse-fasta-primers", "", "A filename for reverse primers")
	// Disable probeFlag until the code is tested
	//flag.Var(&probeFlag, "probe", "`PROBE` is an optional, comma delimited list of DNA sequences that must be present for an allele to be considered a match")
	// Disable development profileFlag
	//flag.StringVar(&profileFlag, "profile", "", "(dev) enable profiling one of `cpu|mem|block`")
	// Disable typeFlag until the code is tested
	flag.StringVar(&typeFlag, "type", "allele", "one of `allele|stat`")
	// TODO: support multiple log levels
	//flag.BoolVar(&debugFlag, "debug", false, "print log messages to stderr")

	// By default, the flag package will render the usage message default value from the type String() method.
	// The PrimerList.String() is useful for debugging, not as a usage message default value.
	// Override the default value with an empty string.
	flag.Lookup("primers").DefValue = ""

	flag.StringVar(&regexJsonFileFlag, "regex-json", "./regex_vipr.json", "file to pull or write regex to for continous use across files")
}

var filename string

func contains(s []text, e text) bool {
  for _, a := range s {
      if string(a) == string(e) {
          return true
      }
  }
  return false
}

func intContained(s []int, e int) bool {
  for _, a := range s {
      if a == e {
          return true
      }
  }
  return false
}

func main() {
	//var err error

	flag.Usage = func() {
		// TODO: Add ./neben help PrimerList subcommand documentation.
		fmt.Fprintf(os.Stderr, "usage: %s --primers PrimerList.txt Assembly.fasta\n", os.Args[0])
		flag.PrintDefaults()
	}
	flag.Parse()
	args := flag.Args()

	//switch strings.ToLower(profileFlag) {
	//case "":
	//case "cpu":
	//	defer profile.Start(profile.CPUProfile, profile.ProfilePath(".")).Stop()
	//case "mem":
	//	defer profile.Start(profile.MemProfile, profile.ProfilePath(".")).Stop()
	//case "block":
	//	defer profile.Start(profile.BlockProfile, profile.ProfilePath(".")).Stop()
	//default:
	//	fatalf("invalid profile: %s\n", profileFlag)
	//}

	if minSequenceLengthFlag < 0 || maxSequenceLengthFlag < 0 {
		fatalf("min and max must be positive integers\n")
	}

	seqToRegexMap = make(map[string][]string)
	if _, err := os.Stat(regexJsonFileFlag); err == nil {
		// read our opened jsonFile as a byte array.
		byteValue, _ := ioutil.ReadFile(regexJsonFileFlag)
		// we unmarshal our byteArray which contains our
		// jsonFile's content into 'seqToRegexMap' which we defined as a flag
		json.Unmarshal(byteValue, &seqToRegexMap)
	}

	// if forward and reverse rpimer flag were set
	// check that the files are valid
	// create a map of primer to name for each primer in primerFlag
	// if primersFlag is empty use these primers to create it

	// if primersFlag is bad, throw err

	//if !debugFlag {
	//log.SetOutput(ioutil.Discard)
	//}

	if len(primersFlag.forward) == 0 || len(primersFlag.reverse) == 0 {
		// If the primer flag is set, its parser is responsible for validating the input.
		// This guard exists in the event the flag was not set.
		// The "flag needs an argument" error message mirrors the error returned by
		// the flag package when a flag is set that expects an argument.
		fatalf("flag needs an argument: -primers\n")
	}

	sequenceChan := make(chan *Contig, contigWorkersFlag)
	indexChan := make(chan Contig, indexWorkersFlag)
	matchChan := make(chan ContigMatch, searchWorkersFlag)

	// Read from either a fasta file
	var sequenceFile *os.File
	stat, _ := os.Stdin.Stat()
	if (stat.Mode() & os.ModeCharDevice) == 0 {
		// data is being piped to stdin
		log.Printf("Reading from stdin\n")
		if len(args) > 0 {
			log.Fatalf("unexpected arguments: %v\n", args)
		}
		sequenceFile = os.Stdin
		filename = "stdin"
	} else {
		// stdin is from a terminal
		var err error
		if len(args) == 0 {
			log.Fatal("expected a fasta file either as a commandline argument or piped through stdin")
		} else if len(args) > 1 {
			log.Fatalf("unexpected arguments: %v\n", args[1:])
		}

		filename = path.Base(args[0])
		log.Printf("Reading from %s\n", filename)
		sequenceFile, err = os.Open(args[0])
		if err != nil {
			log.Fatal(err)
		}
		defer sequenceFile.Close()
	}

	go readFasta(sequenceChan, sequenceFile)

	// Build suffix array
	go suffixarrayWorkers(indexChan, sequenceChan, int(indexWorkersFlag))

	newPrimerFlag := primersFlag

	// Build seq to regex map into the primer lists
	for i := range primersFlag.forward {
		for _, sequence := range primersFlag.forward[i].Sequences {
			// Regex Match
			var regexList []string
			if _, ok := seqToRegexMap[string(sequence)]; ok {
				for _, sequence := range seqToRegexMap[string(sequence)] {
					if strings.Count(sequence,"{") <= maxMismatchFlag {
						if !contains(newPrimerFlag.forward[i].Sequences, text(sequence)) {
							newPrimerFlag.forward[i].Sequences = append(newPrimerFlag.forward[i].Sequences, text(sequence))
						}
					}
				}
				continue
			}else{
				regexList = getAllRegexStr(string(sequence), maxMismatchFlag)
				seqToRegexMap[string(sequence)] = regexList
				for _, sequence := range seqToRegexMap[string(sequence)] {
					if strings.Count(sequence,"{") <= maxMismatchFlag {
						if !contains(newPrimerFlag.forward[i].Sequences, text(sequence)) {
							newPrimerFlag.forward[i].Sequences = append(newPrimerFlag.forward[i].Sequences, text(sequence))
						}
					}
				}
			}
		}
		for _, sequence := range primersFlag.forward[i].RcSequences {
			var regexList []string
			if _, ok := seqToRegexMap[string(sequence)]; ok {
				for _, sequence := range seqToRegexMap[string(sequence)] {
					if strings.Count(sequence,"{") <= maxMismatchFlag {
						if !contains(newPrimerFlag.forward[i].Sequences, text(sequence)) {
							newPrimerFlag.forward[i].RcSequences = append(newPrimerFlag.forward[i].RcSequences, text(sequence))
						}
					}
				}
				continue
			}else{
				regexList = getAllRegexStr(string(sequence), maxMismatchFlag)
				seqToRegexMap[string(sequence)] = regexList
				for _, sequence := range seqToRegexMap[string(sequence)] {
					if strings.Count(sequence,"{") <= maxMismatchFlag {
						if !contains(newPrimerFlag.forward[i].Sequences, text(sequence)) {
							newPrimerFlag.forward[i].RcSequences = append(newPrimerFlag.forward[i].RcSequences, text(sequence))
						}
					}
				}
			}
		}
	}

	for i := range primersFlag.reverse {
		for _, sequence := range primersFlag.reverse[i].Sequences {
			var regexList []string
			if _, ok := seqToRegexMap[string(sequence)]; ok {
				for _, sequence := range seqToRegexMap[string(sequence)] {
					if strings.Count(sequence,"{") <= maxMismatchFlag {
						if !contains(newPrimerFlag.reverse[i].Sequences, text(sequence)) {
							newPrimerFlag.reverse[i].Sequences = append(newPrimerFlag.reverse[i].Sequences, text(sequence))
						}
					}
				}
				continue
			}else{
				regexList = getAllRegexStr(string(sequence), maxMismatchFlag)
				seqToRegexMap[string(sequence)] = regexList
				for _, sequence := range seqToRegexMap[string(sequence)] {
					if strings.Count(sequence,"{") <= maxMismatchFlag {
						if !contains(newPrimerFlag.reverse[i].Sequences, text(sequence)) {
							newPrimerFlag.reverse[i].Sequences = append(newPrimerFlag.reverse[i].Sequences, text(sequence))
						}
					}
				}
			}
		}
		for _, sequence := range primersFlag.reverse[i].RcSequences {
			var regexList []string
			if _, ok := seqToRegexMap[string(sequence)]; ok {
				for _, sequence := range seqToRegexMap[string(sequence)] {
					if strings.Count(sequence,"{") <= maxMismatchFlag {
						if !contains(newPrimerFlag.reverse[i].RcSequences, text(sequence)) {
							newPrimerFlag.reverse[i].RcSequences = append(newPrimerFlag.reverse[i].RcSequences, text(sequence))
						}
					}
				}
				continue
			}else{
				regexList = getAllRegexStr(string(sequence), maxMismatchFlag)
				seqToRegexMap[string(sequence)] = regexList
				for _, sequence := range seqToRegexMap[string(sequence)] {
					if strings.Count(sequence,"{") <= maxMismatchFlag {
						if !contains(newPrimerFlag.reverse[i].RcSequences, text(sequence)) {
							newPrimerFlag.reverse[i].RcSequences = append(newPrimerFlag.reverse[i].RcSequences, text(sequence))
						}
					}
				}
			}
		}
	}

	primersFlag = newPrimerFlag

	jsonStr, jsonErr := json.Marshal(seqToRegexMap)
	if jsonErr != nil {
			fmt.Printf("Error: %s", jsonErr.Error())
	} else {
			_ = ioutil.WriteFile(regexJsonFileFlag, []byte(jsonStr), 0644)
	}

	go matchWorkers(matchChan, indexChan, probeFlag, primersFlag, int(searchWorkersFlag))

	bw := bufio.NewWriter(os.Stdout)
	defer bw.Flush()

	switch typeFlag {
	case "allele":
		writeMatches(bw, matchChan, minSequenceLengthFlag, maxSequenceLengthFlag)
	case "stat":
		writeMatchMatrix(bw, matchChan, primersFlag, minSequenceLengthFlag, maxSequenceLengthFlag)
	default:
		fatalf("type flag value must be either allele or stat\n")
	}
	//fmt.Printf("ViPR Finished!")
}

func writeMatchMatrix(w io.Writer, matchChan chan ContigMatch, primers PrimerList, minLen, maxLen int) {
	// {fwd,rev}PrimerMatchCounter count how many times each primer was found.
	// It does not account for whether the primer actually paired with anything.
	fwdPrimerMatchCounter := make([]int, len(primers.forward))
	revPrimerMatchCounter := make([]int, len(primers.reverse))

	// Initialize a 2D matrix with a counter for each forward/reverse primer pair
	matrix := make([][]int, len(primers.forward))
	for i := range matrix {
		matrix[i] = make([]int, len(primers.reverse))
	}

	// Aggregate matches
	for match := range matchChan {
		log.Printf("aggregate contig results: %s\n", match.contig.descriptor)

		// Count all potential forward primer hits
		for _, fwd := range match.forward {
			fwdPrimerMatchCounter[fwd.primer.Idx] += len(fwd.indices) + len(fwd.rcIndices)
		}
		// Count all potential reverse primer hits
		for _, rev := range match.reverse {
			revPrimerMatchCounter[rev.primer.Idx] += len(rev.indices) + len(rev.rcIndices)
		}

		// Pair every forward/reverse primer combination.
		for _, fwd := range match.forward {
			fLen := len(fwd.primer.Sequence)
			for _, rev := range match.reverse {
				rLen := len(rev.primer.Sequence)
				// Locate all the potential sense strand alleles pairing all the permutations of each primer.
				for _, fIdx := range fwd.indices {
					for _, rIdx := range rev.rcIndices {
						if !match.isAllele(fIdx, fLen, rIdx, rLen, minLen, maxLen, true) {
							continue
						}
						matrix[fwd.primer.Idx][rev.primer.Idx]++

					}
				}
				// Locate all the potential antisense strand alleles pairing all the permutations of each primer.
				for _, fIdx := range fwd.rcIndices {
					for _, rIdx := range rev.indices {
						if !match.isAllele(fIdx, fLen, rIdx, rLen, minLen, maxLen, false) {
							continue
						}
						matrix[fwd.primer.Idx][rev.primer.Idx]++
					}
				}

			}
		}
	}

	// Print summary
	fmt.Fprintf(w, "Source: %s\n", filename)
	fmt.Fprintf(w, "Min: %d\n", minSequenceLengthFlag)
	fmt.Fprintf(w, "Max: %d\n", maxSequenceLengthFlag)

	fmt.Fprintf(w, "\n\nForward primer hits:\n")
	for i := range primers.forward {
		fmt.Fprintf(w, "%s\t%d\n", primers.forward[i].Label, fwdPrimerMatchCounter[i])
	}

	fmt.Fprintf(w, "\n\nReverse primer hits:\n")
	for i := range primers.reverse {
		fmt.Fprintf(w, "%s\t%d\n", primers.reverse[i].Label, revPrimerMatchCounter[i])
	}

	// Print matrix

	// Print header
	fmt.Fprintf(w, "\t")
	for column := range primers.reverse {
		fmt.Fprintf(w, "%s\t", primers.reverse[column].Label)
	}
	fmt.Fprintf(w, "\n")

	// Print rows/columns
	for row := range matrix {
		fmt.Fprintf(w, "%s\t", primers.forward[row].Label)
		for column := range matrix[row] {
			fmt.Fprintf(w, "%d\t", matrix[row][column])
		}
		fmt.Fprintf(w, "\n")
	}
}

func writeMatches(w io.Writer, matchChan chan ContigMatch, minLen, maxLen int) {
	//matches := make(map[string]map[int]struct{})
	for match := range matchChan {
		for _, fwd := range match.forward {
			log.Printf("MATCH FORWARD: %s %s %s %d %d\n", match.contig.descriptor, fwd.primer.Sequence, reverseComplement(fwd.primer.Sequence), len(fwd.indices), len(fwd.rcIndices))
		}
		for _, rev := range match.reverse {
			log.Printf("MATCH REVERSE: %s %s %s %d %d\n", match.contig.descriptor, rev.primer.Sequence, reverseComplement(rev.primer.Sequence), len(rev.indices), len(rev.rcIndices))
		}

		contigIdentifier := []byte(match.contig.descriptor)
		if idx := bytes.IndexByte(contigIdentifier, ' '); idx != -1 {
			contigIdentifier = contigIdentifier[:idx]
		}
		for _, fwd := range match.forward {
			fLen := len(fwd.primer.Sequence)
			for _, rev := range match.reverse {
				rLen := len(rev.primer.Sequence)
				//for fIdx := range fwd.indices {
				for _, fIdx := range fwd.indices {
					//for rIdx := range rev.rcIndices {
					for _, rIdx := range rev.rcIndices {
						if !match.isAllele(fIdx, fLen, rIdx, rLen, minLen, maxLen, true) {
							continue
						}
						start := fIdx
						end := rIdx + len(rev.primer.Sequence)
						fwdStartSeqMatches := fwd.matchRegex[fIdx]
						// fwdEndSeqMatches := fwd.matchRegex[rIdx]
						// fwdRStartSeqMatches := fwd.rMatchRegex[fIdx]
						// fwdREndSeqMatches := fwd.rMatchRegex[rIdx]
						// revStartSeqMatches := rev.matchRegex[fIdx]
						// revEndSeqMatches := rev.matchRegex[rIdx]
						// revRStartSeqMatches := rev.rMatchRegex[fIdx]
						revREndSeqMatches := rev.rMatchRegex[rIdx]
						if (end-start) < 1001 {
							fmt.Fprintf(w, "MATCHES\t%d\t%d\n", fwdStartSeqMatches, revREndSeqMatches)
							fmt.Fprintf(w, "+\t%s\t%s\t%s\t%d\t%d\t%d\t%s\t%s\n\n", fwd.primer.Label, rev.primer.Label, match.contig.sequence[start:end], start, end, end-start, filename, contigIdentifier)
						}
					}
				}
				//for fIdx := range fwd.rcIndices {
				for _, fIdx := range fwd.rcIndices {
					//for rIdx := range rev.indices {
					for _, rIdx := range rev.indices {
						if !match.isAllele(fIdx, fLen, rIdx, rLen, minLen, maxLen, false) {
							continue
						}
						start := rIdx
						end := fIdx + len(fwd.primer.Sequence)
						// fwdStartSeqMatches := fwd.matchRegex[fIdx]
						// fwdEndSeqMatches := fwd.matchRegex[rIdx]
						fwdRStartSeqMatches := fwd.rMatchRegex[fIdx]
						// fwdREndSeqMatches := fwd.rMatchRegex[rIdx]
						// revStartSeqMatches := rev.matchRegex[fIdx]
						revEndSeqMatches := rev.matchRegex[rIdx]
						// revRStartSeqMatches := rev.rMatchRegex[fIdx]
						// revREndSeqMatches := rev.rMatchRegex[rIdx]
						if (end-start) < 1001 {
							fmt.Fprintf(w, "MATCHES\t%d\t%d\n", fwdRStartSeqMatches, revEndSeqMatches)
							rcSequence := reverseComplement(match.contig.sequence[start:end])
							fmt.Fprintf(w, "-\t%s\t%s\t%s\t%d\t%d\t%d\t%s\t%s\n\n", fwd.primer.Label, rev.primer.Label, rcSequence, start, end, end-start, filename, contigIdentifier)
						}
					}
				}
			}
		}
	}

}

type Contig struct {
	descriptor string
	// NOTE: sequence MUST NOT be altered once assigned
	sequence []byte
	index    *suffixarray.Index
}

func NewContig(descriptor string) *Contig {
	return &Contig{
		descriptor: descriptor,
		// The buffer passed to the Write method can change.
		// Explicitly allocate a new buffer so we do not accidentally
		// keep the buffer that was passed in.
		sequence: make([]byte, 0, 4096),
	}
}

func (c *Contig) Free() {
	c.sequence = nil
	c.index = nil
}

func (c *Contig) Write(b []byte) {
	c.sequence = append(c.sequence, b...)
}

func readFasta(sequenceChan chan<- *Contig, r io.Reader) {
	defer close(sequenceChan)

	var line []byte
	var err error
	var contig *Contig

	br := bufio.NewReader(r)

	for line, _, err = br.ReadLine(); err == nil || err != io.EOF; line, _, err = br.ReadLine() {
		if len(line) > 0 && line[0] == '>' {
			if contig != nil {
				// A contig descriptor indicates the end of the prior contig for all except the first contig.
				sequenceChan <- contig
			}
			contig = NewContig(string(line[1:]))
		} else if contig != nil {
			line = bytes.ToUpper(line)
			contig.Write(line)
		}
	}

	if err != io.EOF {
		fatalf("%s\n", err)
	}

	if contig != nil {
		sequenceChan <- contig
	}
}

// PrimerMatch pairs a primer with all the indices where it was found on either for sense or antisense strand of a contig.
// See ContigMatch
type PrimerMatch struct {
	primer Primer
	//indices   map[int]struct{}
	//rcIndices map[int]struct{}
	indices    []int
	matchRegex map[int][]string
	rcIndices  []int
	rMatchRegex map[int][]string
}

func (c PrimerMatch) String() string {
	//return fmt.Sprintf("%s %d %d\n%s\n%s\n", primer.sequence, primerMatch.indices, primerMatch.rcIndices, primer.sequences, primer.rcSequences)
	str, err := json.Marshal(c)
	if err != nil {
		return err.Error()
	}
	return string(str)
}

type ContigMatch struct {
	contig              Contig
	probeForwardIndices []int
	probeReverseIndices []int
	forward             []PrimerMatch
	reverse             []PrimerMatch
}

func NewContigMatch(contig Contig, probes ProbeFlag, primers PrimerList) *ContigMatch {
	//log.Printf("New Contig Match Started\n")
	//log.Printf("Max Mismatches:%d\n",maxMismatchFlag)
	var probeForward, probeReverse []int
	for _, probe := range probes {
		for i := range probe.Sequences {
			regexpFinal, _ := regexp.Compile(string(probe.Sequences[i]))
			//log.Printf("Checking Regex:\n%v\n", regexpFinal)
			mismatchInds := contig.index.FindAllIndex(regexpFinal, -1)
			//log.Printf("Mismatch Ind List:\n%v\n", mismatchInds)
			for _, indxPair := range mismatchInds {
				probeForward = append(probeForward, indxPair[0])
			}
			// //Exact match check
			// probeForward = append(probeForward, contig.index.Lookup(probe.Sequences[i], -1)...)
			// //Regex check for mismatches
			// var regexList []string
			// if _, ok := seqToRegexMap[string(probe.Sequences[i])]; ok {
			// 	regexList = seqToRegexMap[string(probe.Sequences[i])]
			// }else{
			// 	regexList = getAllRegexStr(string(probe.Sequences[i]), maxMismatchFlag)
			// 	seqToRegexMap[string(probe.Sequences[i])] = regexList
			// }
			// //log.Printf("Got rgex of:\n%v\n", regexList)
			// for _, regexStr := range regexList {
			// 	regexpFinal, _ := regexp.Compile(regexStr)
			// 	//log.Printf("Checking Regex:\n%v\n", regexpFinal)
			// 	mismatchInds := contig.index.FindAllIndex(regexpFinal, -1)
			// 	//log.Printf("Mismatch Ind List:\n%v\n", mismatchInds)
			// 	for _, indxPair := range mismatchInds {
			// 		// found := false
			// 		for _, existingLoc := range probeForward {
			// 			if existingLoc == indxPair[0] {
			// 				// found = true
			// 				break
			// 			}
			// 		}
			// 		if true {
			// 			probeForward = append(probeForward, indxPair[0])
			// 		}
			// 	}
			// }

		}

		for i := range probe.RcSequences {
			regexpFinal, _ := regexp.Compile(string(probe.RcSequences[i]))
			//log.Printf("Checking Regex:\n%v\n", regexpFinal)
			mismatchInds := contig.index.FindAllIndex(regexpFinal, -1)
			//log.Printf("Mismatch Ind List:\n%v\n", mismatchInds)
			for _, indxPair := range mismatchInds {
				probeReverse = append(probeReverse, indxPair[0])
			}
			// //Exact Match check
			// probeReverse = append(probeReverse, contig.index.Lookup(probe.RcSequences[i], -1)...)
			//
			// //Regex check for mismatches
			// var regexList []string
			// if _, ok := seqToRegexMap[string(probe.RcSequences[i])]; ok {
			// 	regexList = seqToRegexMap[string(probe.RcSequences[i])]
			// }else{
			// 	regexList = getAllRegexStr(string(probe.RcSequences[i]), maxMismatchFlag)
			// 	seqToRegexMap[string(probe.RcSequences[i])] = regexList
			// }
			// //log.Printf("Got rgex of:\n%v\n", regexList)
			// for _, regexStr := range regexList {
			// 	regexpFinal, _ := regexp.Compile(regexStr)
			// 	//log.Printf("Checking Regex:\n%v\n", regexpFinal)
			// 	mismatchInds := contig.index.FindAllIndex(regexpFinal, -1)
			// 	//log.Printf("Mismatch Ind List:\n%v\n", mismatchInds)
			// 	for _, indxPair := range mismatchInds {
			// 		// found := false
			// 		for _, existingLoc := range probeReverse {
			// 			if existingLoc == indxPair[0] {
			// 				// found = true
			// 				break
			// 			}
			// 		}
			// 		if true {
			// 			probeReverse = append(probeReverse, indxPair[0])
			// 		}
			// 	}
			// }

		}
	}

	// The probe indices MUST be sorted for efficient searching when checking
	// for their presence in a sequence.
	sort.Ints(probeForward)
	sort.Ints(probeReverse)

	m := ContigMatch{
		contig:              contig,
		probeForwardIndices: probeForward,
		probeReverseIndices: probeReverse,
	}

	//log.Printf("Primer search started\n")
	// TODO: is primer search required is probe not found?
	for i := range primers.forward {
		m.addPrimer(primers.forward[i], true)
	}

	for i := range primers.reverse {
		m.addPrimer(primers.reverse[i], false)
	}

	return &m
}

func (c ContigMatch) String() string {
	//return fmt.Sprintf("%s %d %d\n%s\n%s\n", primer.sequence, primerMatch.indices, primerMatch.rcIndices, primer.sequences, primer.rcSequences)
	str, err := json.Marshal(c)
	if err != nil {
		return err.Error()
	}
	return string(str)
}

//This function generates all possible mismatch regex for a given seq deterined by the max-mismatch flag
func getAllRegexStr(sequence string, max_mismatch int) []string {
	var regList []string
	var indexList []int
	startCount := 1
	indexList = append(indexList,startCount)
	//Deletion/Replacment mutation regex list
	for len(indexList) <= max_mismatch && len(indexList) <= len(sequence) {
		regexByteArrRep := sequence
		for indexInd := 0; indexInd < len(indexList); indexInd++ {
			//deletion/replacement
			deleteByteArr := ".{0,1}" //nothing or 1 anything
			indexVal := indexList[indexInd]-1 + indexInd*(len(deleteByteArr)-1)
			if indexVal >= len(regexByteArrRep) {
				regexByteArrRep = regexByteArrRep+deleteByteArr
				break
			}
			regexByteArrRep = regexByteArrRep[0:indexVal]+deleteByteArr+regexByteArrRep[indexVal+1:len(regexByteArrRep)]
		}
		//log.Printf("New Regex: %v",regexByteArrRep)
		if regexByteArrRep != sequence {
			regList = append(regList, regexByteArrRep)
		}


		regexByteArrIn := sequence
		for indexInd := 0; indexInd < len(indexList); indexInd++ {
			//Insertion
			insertByteArr := ".{1,1}"
			indexVal := indexList[indexInd] + indexInd*len(insertByteArr)
			if indexVal < len(regexByteArrIn) {
				regexByteArrIn = regexByteArrIn[0:indexVal]+ insertByteArr+ regexByteArrIn[indexVal:len(regexByteArrIn)]
			}
		}
		//log.Printf("New Insert Regex: %v",regexByteArrIn)
		if regexByteArrIn != sequence {
			regList = append(regList, regexByteArrIn)
		}

		increase := true
		for indexInc := len(indexList)-1; indexInc >= 0; indexInc-- {
			if increase {
				indexList[indexInc] = indexList[indexInc]+1
			}
			if indexList[indexInc] <= len(sequence) {
				increase = false
				break;
			}else{
				startCount = startCount + 1
				indexList[indexInc] = startCount + indexInc
			}
		}
		if increase {
			startCount = 1
			indexList = append(indexList,0)
			for indexInd,_ := range indexList {
				indexList[indexInd] = startCount + indexInd
			}
		}
		//log.Printf("Index List: %v",indexList)
	}
	return regList
}

func (m *ContigMatch) addPrimer(primer Primer, isForward bool) {
	// log.Printf("addPrimer func started\n")
	primerMatch := PrimerMatch{
		primer: primer,
		//indices:   make(map[int]struct{}),
		//rcIndices: make(map[int]struct{}),
		indices:   []int{},
		matchRegex:   make(map[int][]string),
		rcIndices: []int{},
		rMatchRegex:   make(map[int][]string),
	}

	// primerMatch.matchRegex = make(map[int][]string)
	// primerMatch.rMatchRegex = make(map[int][]string)

	for _, sequence := range primer.Sequences {
		//for _, idx := range m.contig.index.Lookup(sequence, -1) {
		//	primerMatch.indices[idx] = struct{}{}
		//}
		// log.Printf("Checking Seq:\n%s\n", sequence)

		regexpFinal, _ := regexp.Compile(string(sequence))
		//log.Printf("Checking Regex:\n%v\n", regexpFinal)
		mismatchInds := m.contig.index.FindAllIndex(regexpFinal, -1)
		//log.Printf("Mismatch Ind List:\n%v\n", mismatchInds)
		for _, indxPair := range mismatchInds {
			if !intContained(primerMatch.indices, indxPair[0]) {
				primerMatch.indices = append(primerMatch.indices, indxPair[0])
			}
			primerMatch.matchRegex[indxPair[0]] = append(primerMatch.matchRegex[indxPair[0]], string(sequence))
		}

		// // Exact match
		// primerMatch.indices = append(primerMatch.indices, m.contig.index.Lookup(sequence, -1)...)
		// // Regex Match
		// var regexList []string
		// if _, ok := seqToRegexMap[string(sequence)]; ok {
		// 	regexList = seqToRegexMap[string(sequence)]
		// }else{
		// 	regexList = getAllRegexStr(string(sequence), maxMismatchFlag)
		// 	seqToRegexMap[string(sequence)] = regexList
		// }
		//
		// // log.Printf("Got rgex of:\n%v\n", regexList)
		// for _, regexStr := range regexList {
		// 	regexpFinal, _ := regexp.Compile(regexStr)
		// 	//log.Printf("Checking Regex:\n%v\n", regexpFinal)
		// 	mismatchInds := m.contig.index.FindAllIndex(regexpFinal, -1)
		// 	//log.Printf("Mismatch Ind List:\n%v\n", mismatchInds)
		// 	// Only include new indexes to avoid repeats
		// 	for _, indxPair := range mismatchInds {
		// 		// found := false
		// 		for _, existingLoc := range primerMatch.indices {
		// 			if existingLoc == indxPair[0] {
		// 				// found = true
		// 				break
		// 			}
		// 		}
		// 		if true {
		// 			// log.Printf("New Mismatch found:\n%v\n", regexpFinal)
		// 			primerMatch.indices = append(primerMatch.indices, indxPair[0])
		// 		}
		// 	}
		// }

	}

	for _, sequence := range primer.RcSequences {
		//for _, idx := range m.contig.index.Lookup(sequence, -1) {
		//	primerMatch.rcIndices[idx] = struct{}{}
		//}

		regexpFinal, _ := regexp.Compile(string(sequence))
		//log.Printf("Checking Regex:\n%v\n", regexpFinal)
		mismatchInds := m.contig.index.FindAllIndex(regexpFinal, -1)
		//log.Printf("Mismatch Ind List:\n%v\n", mismatchInds)
		for _, indxPair := range mismatchInds {
			if !intContained(primerMatch.rcIndices, indxPair[0]) {
				primerMatch.rcIndices = append(primerMatch.rcIndices, indxPair[0])
			}
			primerMatch.rMatchRegex[indxPair[0]] = append(primerMatch.rMatchRegex[indxPair[0]], string(sequence))
		}

		// // Exact Match
		// primerMatch.rcIndices = append(primerMatch.rcIndices, m.contig.index.Lookup(sequence, -1)...)
		// // Regex Match
		// var regexList []string
		// if _, ok := seqToRegexMap[string(sequence)]; ok {
		// 	regexList = seqToRegexMap[string(sequence)]
		// }else{
		// 	regexList = getAllRegexStr(string(sequence), maxMismatchFlag)
		// 	seqToRegexMap[string(sequence)] = regexList
		// }
		// // log.Printf("Got rgex of:\n%v\n", regexList)
		// for _, regexStr := range regexList {
		// 	regexpFinal, _ := regexp.Compile(regexStr)
		// 	//log.Printf("Checking Regex:\n%v\n", regexpFinal)
		// 	mismatchInds := m.contig.index.FindAllIndex(regexpFinal, -1)
		// 	//log.Printf("Mismatch Ind List:\n%v\n", mismatchInds)
		// 	// Only include new indexes to avoid repeats
		// 	for _, indxPair := range mismatchInds {
		// 		// found := false
		// 		for _, existingLoc := range primerMatch.rcIndices {
		// 			if existingLoc == indxPair[0] {
		// 				// found = true
		// 				break
		// 			}
		// 		}
		// 		if true {
		// 			// log.Printf("New Mismatch found:\n%v\n", regexpFinal)
		// 			primerMatch.rcIndices = append(primerMatch.rcIndices, indxPair[0])
		// 		}
		// 	}
		// }
	}

	if isForward {
		m.forward = append(m.forward, primerMatch)
	} else {
		m.reverse = append(m.reverse, primerMatch)
	}
}

func (m *ContigMatch) isAllele(fIdx, fLen, rIdx, rLen, min, max int, isForward bool) bool {
	// sequenceLength is the length of the entire allele sequence including the forward/reverse primer
	var startIdx, startLen, endIdx, endLen int

	if isForward {
		startIdx, startLen, endIdx, endLen = fIdx, fLen, rIdx, rLen
	} else {
		startIdx, startLen, endIdx, endLen = rIdx, rLen, fIdx, fLen
	}

	sequenceLength := endIdx + endLen - startIdx

	// An allele is assumed to be any sequence:
	// - bounded on either side by any forward/reverse primer pair
	// - contains any probe sequence
	// - does not exceed the users min/max lengths
	return startIdx < endIdx &&
		sequenceLength > min &&
		sequenceLength < max &&
		//(len(probeFlag) == 0 || m.sequenceContainsProbe(startIdx+startLen, endIdx-1, isForward))
		m.sequenceContainsProbe(startIdx+startLen, endIdx-1, isForward)
}

// sequenceContainsProbe returns true if a probe was found within the (inclusive) range of position.
// the start/end index boundaries of a sequence in the contig.
func (m *ContigMatch) sequenceContainsProbe(startIdx, endIdx int, isForward bool) bool {
	var probeIndices []int

	// FIXME: panic if the start/end indices are not valid indices of the contig sequence.
	// if startIdx > endIdx || startIdx < 0 || endIdx < 0 {
	//	panic(fmt.Sprintf("invalid sequence range %d - %d in contig %s length %d\n",
	//	startIdx, endIdx, m.contig.descriptor, len(contig.sequence)))
	//}

	if isForward {
		probeIndices = m.probeForwardIndices
	} else {
		probeIndices = m.probeReverseIndices
	}

	// If the user did not specify a probe, it does not matter if the sequence contains a probe.
	// This guards against an index out of bounds error.
	if len(probeIndices) == 0 {
		return true
	}

	// Assuming the probe indices are sorted, find the first probe index that at least matches the startIdx position.
	// Return the length of the array if no match is found.
	i := sort.Search(len(probeIndices), func(i int) bool {
		return startIdx <= probeIndices[i]
	})

	// FIXME: The function finds a probe that starts within the sequence, but no where does it check the end position.
	// The function does not verify the probe is contained within the sequence.
	return i < len(probeIndices) && probeIndices[i] <= endIdx
}

func suffixarrayWorkers(indexChan chan<- Contig, sequenceChan <-chan *Contig, threads int) {
	defer close(indexChan)
	var wg sync.WaitGroup
	wg.Add(threads / 2)
	for i := 0; i < threads/2; i++ {
		go func(indexChan chan<- Contig, sequenceChan <-chan *Contig, i int) {
			defer wg.Done()
			for contig := range sequenceChan {
				log.Printf("Start index %s %d/%d\n", contig.descriptor, len(sequenceChan), cap(sequenceChan))
				start := time.Now()
				contig.index = suffixarray.New(contig.sequence)
				log.Printf("End index %s %d %fs\n", contig.descriptor, len(contig.sequence), time.Since(start).Seconds())
				indexChan <- *contig
			}
			log.Printf("Shutdown index worker %d\n", i)
		}(indexChan, sequenceChan, i)
	}
	wg.Wait()
	log.Println("Shutdown index WaitGroup")
}

func matchWorkers(matchChan chan ContigMatch, indexChan <-chan Contig, probes ProbeFlag, primers PrimerList, threads int) {
	defer close(matchChan)
	var wg sync.WaitGroup
	wg.Add(threads)
	for i := 0; i < threads; i++ {
		go func(matchChan chan<- ContigMatch, indexChan <-chan Contig, primers PrimerList, i int) {
			defer wg.Done()
			for contig := range indexChan {
				log.Printf("Start match %s %d/%d\n", contig.descriptor, len(indexChan), cap(indexChan))
				start := time.Now()

				matchChan <- *NewContigMatch(contig, probes, primers)

				log.Printf("End match %s %d %fs\n", contig.descriptor, len(contig.sequence), time.Since(start).Seconds())
			}
			log.Printf("Shutdown match worker %d\n", i)
		}(matchChan, indexChan, primers, i)

	}

	wg.Wait()
	log.Println("Shutdown match WaitGroup")
}

// text is an alias for byte because the json package encodes []byte as a base64-encoded string by default.
type text []byte

func (t text) String() string {
	return string(t)
}

// MarshalJSON returns text as a quoted string.
func (t text) MarshalJSON() ([]byte, error) {
	// Duplicate the text with padding for the quotes
	x := make([]byte, len(t)+2)
	copy(x[1:], t)
	x[0] = '"'
	x[len(x)-1] = '"'
	return x, nil
}

type Primer struct {
	// TODO: a unique sequential number assigned to each primer for the purpose of sorting matches
	Idx int `json:"idx"`
	// label is an optional description
	Label text `json:"label"`
	// sequence is the original sequence with degenerarcies
	Sequence text `json:"sequence"`
	// sequences are the sequence permutations with degeneracies expanded
	Sequences []text `json:"sequences"`
	// rcSequences are the reverse complement of sequences
	RcSequences []text `json:"rcSequences"`
}

func expandDegeneratePosition(primers []text, position int, l ...byte) []text {
	for j := range primers {
		primers[j][position] = l[0]
	}
	primers_len := len(primers)
	for _, m := range l[1:] {
		for j := 0; j < primers_len; j++ {
			p := make([]byte, len(primers[j]))
			copy(p, primers[j])
			p[position] = m
			primers = append(primers, p)
		}
	}
	return primers
}

func expandDegenerateSequence(sequence []byte) []text {
	var primers []text

	if len(sequence) == 0 {
		return nil
	}

	primer := make(text, len(sequence))
	copy(primer, sequence)

	primers = append(primers, primer)
	for i, nt := range sequence {
		switch nt {
		default:
			// TODO: should unrecognized characters panic?
			//primers = expandDegeneratePosition(primers, i, nt)
			// This could return an error.
			panic(ErrInvalidNucleotide{
				position: i + 1,
				sequence: string(sequence),
			})
		case 'A', 'C', 'G', 'T', 'U':
			primers = expandDegeneratePosition(primers, i, nt)
		case 'W':
			primers = expandDegeneratePosition(primers, i, 'A', 'T')
		case 'S':
			primers = expandDegeneratePosition(primers, i, 'G', 'C')
		case 'M':
			primers = expandDegeneratePosition(primers, i, 'A', 'C')
		case 'K':
			primers = expandDegeneratePosition(primers, i, 'G', 'T')
		case 'R':
			primers = expandDegeneratePosition(primers, i, 'A', 'G')
		case 'Y':
			primers = expandDegeneratePosition(primers, i, 'C', 'T')
		case 'B':
			primers = expandDegeneratePosition(primers, i, 'C', 'G', 'T')
		case 'D':
			primers = expandDegeneratePosition(primers, i, 'A', 'G', 'T')
		case 'H':
			primers = expandDegeneratePosition(primers, i, 'A', 'C', 'T')
		case 'V':
			primers = expandDegeneratePosition(primers, i, 'A', 'C', 'G')
		case 'N':
			primers = expandDegeneratePosition(primers, i, 'G', 'A', 'T', 'C')
		case '.', '-':
			panic(fmt.Errorf("%q is a valid IUPAC Nucleotide Code, but not a valid nucleotide at position %d in sequence %s", nt, i+1, sequence))
		}
	}
	return primers
}

type PrimerList struct {
	forward []Primer
	reverse []Primer
}

func (p *PrimerList) Append(sequence, label []byte, isForwardPrimer bool) error {
	// TODO: Should the sequence be appended if it is a duplicate?
	// NOTE: There should be no side-effects modifying the sequence/label parameters.

	// Normalize the sequence as capital letters.
	sequence = bytes.ToUpper(sequence)

	// If the user does not specify a label, the sequence is the default label.
	// It is possible to keep the original sequence capitalization by explicitly
	// copying the sequence before normalizing, but I decided to enforce all caps.
	if label == nil {
		label = sequence
	}

	// Validate sequence contains only characters in the IUPAC nucleotide code character set.
	for i, nt := range sequence {
		// TODO: Should '-' be replaced with 'N'?
		if idx := bytes.IndexByte([]byte(IupacNucleotideCode), nt); idx == -1 {
			return ErrInvalidNucleotide{
				position: i + 1,
				sequence: string(sequence),
			}
		}
	}

	primer := Primer{
		Label:       label,
		Sequence:    sequence,
		Sequences:   expandDegenerateSequence(sequence),
		RcSequences: expandDegenerateSequence(reverseComplement(sequence)),
	}

	if isForwardPrimer {
		primer.Idx = len(p.forward)
		p.forward = append(p.forward, primer)
		if s, err := json.Marshal(primer); err == nil {
			log.Printf("Add forward primer: %s\n", s)
		} else {
			return fmt.Errorf("TODO: error context %s", err)
		}
	} else {
		primer.Idx = len(p.reverse)
		p.reverse = append(p.reverse, primer)
		if s, err := json.Marshal(primer); err == nil {
			log.Printf("Add reverse primer: %s\n", s)
		} else {
			return fmt.Errorf("TODO: error context %s", err)
		}
	}

	return nil
}

func (p *PrimerList) Read(r io.Reader) error {
	var isForwardPrimer = true
	var line, sequence, label []byte

	// TODO: validate primer character-set
	scanner := bufio.NewScanner(r)
	for scanner.Scan() {
		if len(scanner.Bytes()) == 0 {
			// It is assumed the primer list is a file with
			// the forward reads listed on separate lines followed by
			// a blank line followed by
			// the reverse reads listed on separate lines
			// FIXME: Because the flag is toggled with every empty line,
			// it may toggle more than intended if the user uses more than
			// one line to delimit the forward/reverse primers.
			isForwardPrimer = !isForwardPrimer
			continue
		}

		// Split line into primer sequence and label
		line = make([]byte, len(scanner.Bytes()))
		copy(line, scanner.Bytes())

		if idx := bytes.IndexFunc(line, unicode.IsSpace); idx != -1 {
			sequence = bytes.ToUpper(line[:idx])
			label = bytes.TrimSpace(line[idx:])
		} else {
			sequence = bytes.ToUpper(line)
			label = line
		}

		if err := p.Append(sequence, label, isForwardPrimer); err != nil {
			return err
		}
	}

	if err := scanner.Err(); err != nil {
		return err
	}

	// TODO: unit test entire file is read (especially the last line), all combinations are built, if raises an error if the character-set is invalid or both forward/reverse primers are not present (format error or empty file)

	if len(p.forward) == 0 || len(p.reverse) == 0 {
		return ErrInvalidFormat("at least one forward and one reverse primer sequence is required")
	}

	return nil
}

func (p *PrimerList) Set(s string) error {
	// TODO: return error if both forward and reverse primers are not present

	if len(s) == 0 {
		return errors.New("expected a filename or list of forward/reverse primers")
	}

	// First try reading as a file containing a list of primers.
	if primerListFile, err := os.Open(s); err == nil {
		defer primerListFile.Close()
		if err := p.Read(primerListFile); err != nil {
			return err
		}
		return nil
	}

	// Fallback to reading as a list of primers
	// It is assumed the string will be:
	// - a comma-delimited list of forward primers
	// - a colon to delimit the forward/reverse primer lists
	// - a comma-delimited list of reverse primers
	//
	// example: GATC,ATCG:GGGG,ATGC,CCTA

	// Split the string into forward/reverse primers
	primers := bytes.Split([]byte(s), []byte(":"))
	if len(primers) != 2 {
		return errors.New("expected a colon delimiting the list of forward primers from the reverse primers")
	}

	forward := bytes.Split(primers[0], []byte(","))
	if len(forward) == 1 && len(forward[0]) == 0 {
		return errors.New("the list of forward primers is empty")
	}
	for i := range forward {
		if err := p.Append(forward[i], forward[i], true); err != nil {
			return err
		}
	}

	reverse := bytes.Split(primers[1], []byte(","))
	if len(reverse) == 1 && len(reverse[0]) == 0 {
		return errors.New("the list of reverse primers is empty")
	}
	for i := range reverse {
		if err := p.Append(reverse[i], reverse[i], false); err != nil {
			return err
		}
	}

	return nil
}

func (p *PrimerList) String() string {
	//return ""
	forwardLabels := make([][]byte, len(p.forward))
	for i := range p.forward {
		forwardLabels[i] = p.forward[i].Label
	}
	reverseLabels := make([][]byte, len(p.reverse))
	for i := range p.reverse {
		reverseLabels[i] = p.reverse[i].Label
	}
	return fmt.Sprintf("{forward: %q, reverse: %q}", forwardLabels, reverseLabels)
}

type Probe struct {
	Sequence    text   `json:"sequence"`
	Sequences   []text `json:"sequences"`
	RcSequences []text `json:"rcSequences"`
}

type ProbeFlag []Probe

func (p *ProbeFlag) Set(s string) error {
	probes := bytes.Split([]byte(s), []byte(","))

	for i := range probes {
		sequence, err := iupacNucleotideSequence(probes[i])
		if err != nil {
			return err
		}

		probe := Probe{
			Sequence:    sequence,
			Sequences:   expandDegenerateSequence(sequence),
			RcSequences: expandDegenerateSequence(reverseComplement(sequence)),
		}

		if s, err := json.Marshal(probe); err == nil {
			log.Printf("Add probe: %s\n", s)
		} else {
			return fmt.Errorf("TODO: error context %s", err)
		}

		*p = append(*p, probe)
	}

	return nil
}

func (p *ProbeFlag) String() string {
	this := *p
	b := make([][]byte, len(this))
	for i := range this {
		b[i] = this[i].Sequence
	}
	return string(bytes.Join(b, []byte(",")))
}

var complementLookupTable = [256]uint8{
	'A': 'T', 'a': 'T',
	'C': 'G', 'c': 'G',
	'G': 'C', 'g': 'C',
	'T': 'A', 't': 'A',
	'U': 'A', 'u': 'A',
	'M': 'K', 'm': 'K',
	'R': 'Y', 'r': 'Y',
	'W': 'W', 'w': 'W',
	'S': 'S', 's': 'S',
	'Y': 'R', 'y': 'R',
	'K': 'M', 'k': 'M',
	'V': 'B', 'v': 'B',
	'H': 'D', 'h': 'D',
	'D': 'H', 'd': 'H',
	'B': 'V', 'b': 'V',
	'N': 'N', 'n': 'N',
	'-': '-', '.': '.',
}

// reverseComplement returns the reverse complement of a nucleotide sequence.
func reverseComplement(s []byte) []byte {
	sLen := len(s)
	rc := make([]byte, sLen)
	for i := 0; i < sLen; i++ {
		rc[i] = complementLookupTable[s[sLen-i-1]]
		// FIXME: all invalid nucleotide character codes will silently be replaced with
		// a zero value. Should this trigger a panic or throw an error?
		//if rc[i] == 0 {
		//	panic(fmt.Errorf(""))
		//}
	}
	return rc
}

type ErrInvalidFormat string

func (e ErrInvalidFormat) Error() string {
	return string(e)
}

func iupacNucleotideSequence(s []byte) (text, error) {
	s = bytes.ToUpper(s)

	for i, nt := range s {
		if strings.IndexByte(IupacNucleotideCode, nt) == -1 {
			return nil, ErrInvalidNucleotide{
				position: i + 1,
				sequence: string(s),
			}
		}
	}

	return text(s), nil
}

type ErrInvalidNucleotide struct {
	position int
	sequence string
}

func (e ErrInvalidNucleotide) Error() string {
	if e.position == 0 || len(e.sequence) == 0 {
		panic(fmt.Errorf("TODO: ErrInvalidSequence improperly instantiated: %#v", e))
	}
	return fmt.Sprintf("invalid nucleotide %q at position %d in sequence %s", e.sequence[e.position-1], e.position, e.sequence)
}

func fatalf(fmtMsg string, a ...interface{}) {
	fmt.Fprintf(os.Stderr, fmtMsg, a...)
	flag.Usage()
	os.Exit(1)
}
