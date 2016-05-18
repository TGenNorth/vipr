package main

import (
	"bufio"
	"bytes"
	"encoding/json"
	"flag"
	"fmt"
	"io/ioutil"
	"path"
	"unicode"

	"github.com/corburn/neben/fmi"
	//"github.com/pkg/profile"

	"io"
	"log"
	"os"
	//"strings"
	"sync"
	"time"
)

var (
	threadsFlag     uint
	maxMismatchFlag int
	maxSequenceFlag int
	primersFlag     string
	debugFlag       bool
)

func init() {
	flag.UintVar(&threadsFlag, "threads", 10, "")
	// TODO: do not allow negative
	flag.IntVar(&maxMismatchFlag, "max-mismatch", 0, "")
	flag.IntVar(&maxSequenceFlag, "max-sequence", 200, "")
	flag.StringVar(&primersFlag, "primers", "", "path to a file listing forward followed by reverse primers")
	flag.BoolVar(&debugFlag, "debug", false, "print log messages to stderr")
}

var filename string

func main() {
	flag.Usage = func() {
		fmt.Fprintf(os.Stderr, "Usage of %s:\n", os.Args[0])
		flag.PrintDefaults()
	}
	flag.Parse()
	//defer profile.Start(profile.CPUProfile, profile.ProfilePath(".")).Stop()
	args := flag.Args()

	if !debugFlag {
		log.SetOutput(ioutil.Discard)
	}

	if len(args) != 1 {
		flag.Usage()
		os.Exit(1)
	}

	threads := int(threadsFlag)
	// TODO: validate arguments
	sequenceFilename := args[0]
	filename = path.Base(sequenceFilename)
	primerListFilename := primersFlag
	sequenceChan := make(chan *Contig, threads/2)
	indexChan := make(chan Contig, threads)
	matchChan := make(chan ContigMatch, threads)

	primerListFile, err := os.Open(primerListFilename)
	if err != nil {
		log.Fatal(err)
	}
	defer primerListFile.Close()
	var primers PrimerList
	if err := primers.Read(primerListFile); err != nil {
		log.Fatal(err)
	}

	// Parse fasta into contigs
	sequenceFile, err := os.Open(sequenceFilename)
	if err != nil {
		log.Fatal(err)
	}
	defer sequenceFile.Close()

	go readFasta(sequenceChan, sequenceFile)

	// Build suffix array
	go suffixarrayWorkers(indexChan, sequenceChan, threads)

	go matchWorkers(matchChan, indexChan, primers, threads)

	writeMatches(matchChan)
}

func writeMatches(matchChan chan ContigMatch) {
	//matches := make(map[string]map[int]struct{})
	bw := bufio.NewWriter(os.Stdout)
	defer bw.Flush()
	for match := range matchChan {
		/*for _, fwd := range match.forward {
			log.Printf("MATCH FORWARD: %s %s %s %d %d\n", match.contig.descriptor, fwd.primer.sequence, reverseComplement(fwd.primer.sequence), len(fwd.indices), len(fwd.rcIndices))
		}
		for _, rev := range match.reverse {
			log.Printf("MATCH REVERSE: %s %s %s %d %d\n", match.contig.descriptor, rev.primer.sequence, reverseComplement(rev.primer.sequence), len(rev.indices), len(rev.rcIndices))
		}*/

		contigIdentifier := []byte(match.contig.descriptor)
		if idx := bytes.IndexByte(contigIdentifier, ' '); idx != -1 {
			contigIdentifier = contigIdentifier[:idx]
		}
		for _, fwd := range match.forward {
			for _, rev := range match.reverse {
				for fIdx := range fwd.indices {
					for rIdx := range rev.rcIndices {
						sequenceLength := rIdx + len(rev.primer.Sequence) - fIdx
						if fIdx > rIdx || (maxSequenceFlag > 0 && sequenceLength > maxSequenceFlag) {
							continue
						}
						start := fIdx
						end := rIdx + len(rev.primer.Sequence)
						fmt.Fprintf(bw, "+\t%s\t%s\t%s\t%d\t%d\t%d\t%s\t%s\n", fwd.primer.Label, rev.primer.Label, match.contig.sequence[start:end], start, end, end-start, filename, contigIdentifier)

					}
				}
				for fIdx := range fwd.rcIndices {
					for rIdx := range rev.indices {
						sequenceLength := fIdx + len(fwd.primer.Sequence) - rIdx
						if rIdx > fIdx || (maxSequenceFlag > 0 && sequenceLength > maxSequenceFlag) {
							continue
						}
						start := rIdx
						end := fIdx + len(fwd.primer.Sequence)
						rcSequence, err := reverseComplement(match.contig.sequence[start:end])
						if err != nil {
							// TODO: add context to error
							log.Fatal(err)
						}
						fmt.Fprintf(bw, "-\t%s\t%s\t%s\t%d\t%d\t%d\t%s\t%s\n", fwd.primer.Label, rev.primer.Label, rcSequence, start, end, end-start, filename, contigIdentifier)
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
	//index    *suffixarray.Index
	index *fmi.Index
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

/*type Fasta struct {
	r  io.Reader
	br bufio.Reader
}

func (f *Fasta) ReadContig() (Contig, error) {
	var line []byte
	var err error
	var contig *Contig

	var id uint
	for line, _, err = f.br.ReadLine(); err == nil || err != io.EOF; line, _, err = f.br.ReadLine() {
		line = bytes.ToUpper(line)
		if line[0] == '>' {
			if contig != nil {
				// A contig descriptor indicates the end of the prior contig for all except the first contig.
				sequenceChan <- contig
			}
			contig = NewContig(string(line[1:]))
		} else {
			contig.Write(line)
		}
	}

	if err != io.EOF {
		return err
	}

	return *contig, nil
}*/

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
		log.Fatal(err)
	}

	if contig != nil {
		sequenceChan <- contig
	}
}

type PrimerMatch struct {
	primer    Primer
	indices   map[int]struct{}
	rcIndices map[int]struct{}
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
	contig  Contig
	forward []PrimerMatch
	reverse []PrimerMatch
}

func (c ContigMatch) String() string {
	//return fmt.Sprintf("%s %d %d\n%s\n%s\n", primer.sequence, primerMatch.indices, primerMatch.rcIndices, primer.sequences, primer.rcSequences)
	str, err := json.Marshal(c)
	if err != nil {
		return err.Error()
	}
	return string(str)
}

func (m *ContigMatch) addPrimer(primer Primer, isForward bool) {
	primerMatch := PrimerMatch{
		primer:    primer,
		indices:   make(map[int]struct{}),
		rcIndices: make(map[int]struct{}),
	}

	for _, sequence := range primer.Sequences {
		for _, idx := range m.contig.index.Lookup(sequence, maxMismatchFlag) {
			primerMatch.indices[idx] = struct{}{}
		}
	}

	for _, sequence := range primer.RcSequences {
		for _, idx := range m.contig.index.Lookup(sequence, maxMismatchFlag) {
			primerMatch.rcIndices[idx] = struct{}{}
		}
	}

	if isForward {
		m.forward = append(m.forward, primerMatch)
	} else {
		m.reverse = append(m.reverse, primerMatch)
	}
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
				//contig.index = suffixarray.New(contig.sequence)
				contig.index = fmi.New(contig.sequence)
				indexChan <- *contig
				log.Printf("End index %s %d %fs\n", contig.descriptor, len(contig.sequence), time.Since(start).Seconds())
			}
			log.Printf("Shutdown index worker %d\n", i)
		}(indexChan, sequenceChan, i)
	}
	wg.Wait()
	log.Println("Shutdown index WaitGroup")
}

func matchWorkers(matchChan chan ContigMatch, indexChan <-chan Contig, primers PrimerList, threads int) {
	defer close(matchChan)
	var wg sync.WaitGroup
	wg.Add(threads)
	for i := 0; i < threads; i++ {
		go func(matchChan chan<- ContigMatch, indexChan <-chan Contig, primers PrimerList, i int) {
			defer wg.Done()
			for contig := range indexChan {
				log.Printf("Start match %s %d/%d\n", contig.descriptor, len(indexChan), cap(indexChan))
				start := time.Now()

				contigMatch := ContigMatch{
					contig: contig,
				}

				log.Printf("Scan FORWARD primers in %s\n", contig.descriptor)
				for _, primer := range primers.forward {
					contigMatch.addPrimer(primer, true)
				}
				log.Printf("Scan REVERSE primers in %s\n", contig.descriptor)
				for _, primer := range primers.reverse {
					contigMatch.addPrimer(primer, false)
				}

				matchChan <- contigMatch

				log.Printf("End match %s %d %fs\n", contig.descriptor, len(contig.sequence), time.Since(start).Seconds())
			}
			log.Printf("Shutdown match worker %d\n", i)
		}(matchChan, indexChan, primers, i)

	}

	wg.Wait()
	log.Println("Shutdown match WaitGroup")
}

// text is an alias for byte because []byte encodes as a base64-encoded string
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
	// label is an optional description
	Label text `json:"label"`
	// sequence is the original sequence with degenerarcies
	Sequence text `json:"sequence"`
	// sequences are the sequence permutations with degeneracies expanded
	Sequences []text `json:"sequences"`
	// rcSequences are the reverse complement of sequences
	RcSequences []text `json:"rcSequences"`
}

// String is for unit test error messages
func (p Primer) String() string {
	if b, err := json.Marshal(p); err != nil {
		return err.Error()
	} else {
		return string(b)
	}
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
			panic(fmt.Errorf("error expanding primer sequence: %s ", sequence))
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
		case 'N', '-':
			primers = expandDegeneratePosition(primers, i, 'G', 'A', 'T', 'C')
		}
	}
	return primers
}

type PrimerList struct {
	forward []Primer
	reverse []Primer
}

func (p PrimerList) String() string {
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

func (p *PrimerList) Append(sequence, label []byte, isForwardPrimer bool) error {
	if label == nil {
		label = sequence
	}

	rcSequence, err := reverseComplement(sequence)
	if err != nil {
		return err
	}

	primer := Primer{
		Label:       label,
		Sequence:    sequence,
		Sequences:   expandDegenerateSequence(sequence),
		RcSequences: expandDegenerateSequence(rcSequence),
	}

	if isForwardPrimer {
		p.forward = append(p.forward, primer)
		log.Printf("Add forward primer:\n%s\n", primer)
	} else {
		p.reverse = append(p.reverse, primer)
		log.Printf("Add reverse primer:\n%s\n", primer)
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

var reverseComplementTable = map[byte]byte{
	'A': 'T',
	'C': 'G',
	'G': 'C',
	'T': 'A',
	'U': 'A',
	'M': 'K',
	'R': 'Y',
	'W': 'W',
	'S': 'S',
	'Y': 'R',
	'K': 'M',
	'V': 'B',
	'H': 'D',
	'D': 'H',
	'B': 'V',
	'N': 'N',
}

type ErrInvalidFormat string

func (e ErrInvalidFormat) Error() string {
	return string(e)
}

type ErrInvalidSequence string

func (e ErrInvalidSequence) Error() string {
	return string(e)
}

func reverseComplement(s []byte) ([]byte, error) {
	sLen := len(s)
	rc := make([]byte, sLen)
	for i := 0; i < sLen; i++ {
		idx := sLen - i - 1
		if c, ok := reverseComplementTable[s[idx]]; ok {
			rc[i] = c
		} else {
			return nil, ErrInvalidSequence(fmt.Sprintf("unrecognized nucleotide %q at index %d in sequence %q", s[idx], idx, s))
		}
	}
	return rc, nil
}
