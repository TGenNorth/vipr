package main

import (
	"bufio"
	"bytes"
	"flag"
	"fmt"
	"path"
	//"github.com/pkg/profile"
	"index/suffixarray"
	"io"
	"log"
	"os"
	//"strings"
	"sync"
	"time"
)

var (
	threadsFlag     uint
	maxSequenceFlag uint
	primersFlag     string
)

func init() {
	flag.UintVar(&threadsFlag, "threads", 10, "")
	flag.UintVar(&maxSequenceFlag, "max-sequence", 200, "")
	flag.StringVar(&primersFlag, "primers", "", "")
}

var filename string

func main() {
	flag.Parse()
	//defer profile.Start(profile.CPUProfile, profile.ProfilePath(".")).Stop()
	args := flag.Args()

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
	primers.Read(primerListFile)

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

	//matches := make(map[string]map[int]struct{})
	bw := bufio.NewWriter(os.Stdout)
	defer bw.Flush()
	for match := range matchChan {
		for _, fwd := range match.forward {
			log.Printf("MATCH FORWARD: %s %s %s %d %d\n", match.contig.identifier, fwd.primer.sequence, reverseComplement(fwd.primer.sequence), len(fwd.indices), len(fwd.rcIndices))
		}
		for _, rev := range match.reverse {
			log.Printf("MATCH REVERSE: %s %s %s %d %d\n", match.contig.identifier, rev.primer.sequence, reverseComplement(rev.primer.sequence), len(rev.indices), len(rev.rcIndices))
		}

		contigIdentifier := []byte(match.contig.identifier)
		if idx := bytes.IndexByte(contigIdentifier, ' '); idx != -1 {
			contigIdentifier = contigIdentifier[:idx]
		}
		for _, fwd := range match.forward {
			for _, rev := range match.reverse {
				for _, fIdx := range fwd.indices {
					for _, rIdx := range rev.rcIndices {
						sequenceLength := rIdx + len(rev.primer.sequence) - fIdx
						if fIdx > rIdx || sequenceLength > 200 {
							continue
						}
						start := fIdx
						end := rIdx + len(rev.primer.sequence)
						fmt.Fprintf(bw, "F\t%s\t%s\t%s\t%d\t%d\t%d\t%s\t%s\n", fwd.primer.sequence, rev.primer.sequence, match.contig.sequence[start:end], start, end, end-start, filename, contigIdentifier)

					}
				}
				for _, fIdx := range fwd.rcIndices {
					for _, rIdx := range rev.indices {
						sequenceLength := fIdx + len(fwd.primer.sequence) - rIdx
						if rIdx > fIdx || sequenceLength > 200 {
							continue
						}
						start := rIdx
						end := fIdx + len(fwd.primer.sequence)
						fmt.Fprintf(bw, "R\t%s\t%s\t%s\t%d\t%d\t%d\t%s\t%s\n", fwd.primer.sequence, rev.primer.sequence, reverseComplement(match.contig.sequence[start:end]), start, end, end-start, filename, contigIdentifier)
					}
				}
			}
		}
	}
}

type Contig struct {
	identifier string
	// NOTE: sequence MUST NOT be altered once assigned
	sequence []byte
	index    *suffixarray.Index
}

func NewContig(identifier string) *Contig {
	return &Contig{
		identifier: identifier,
		// The buffer passed to the Write method can change.
		// Explicitly allocate a new buffer so we do not accidentally
		// keep the buffer that was passed in.
		sequence: make([]byte, 0, 4096),
		//sequence:   bufferPool.Get().([]byte)[:0],
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
				// A contig identifier indicates the end of the prior contig for all except the first contig.
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
		line = bytes.ToUpper(line)
		if line[0] == '>' {
			if contig != nil {
				// A contig identifier indicates the end of the prior contig for all except the first contig.
				sequenceChan <- contig
			}
			contig = NewContig(string(line[1:]))
		} else {
			contig.Write(line)
		}
	}

	if err != io.EOF {
		log.Fatal(err)
	}

	sequenceChan <- contig
}

type PrimerMatch struct {
	primer    Primer
	indices   []int
	rcIndices []int
}

type ContigMatch struct {
	contig  Contig
	forward []PrimerMatch
	reverse []PrimerMatch
}

func (m *ContigMatch) addPrimer(primer Primer, isForward bool) {
	primerMatch := PrimerMatch{
		primer: primer,
	}

	for _, sequence := range primer.sequences {
		primerMatch.indices = append(primerMatch.indices, m.contig.index.Lookup(sequence, -1)...)
	}

	for _, sequence := range primer.rcSequences {
		primerMatch.rcIndices = append(primerMatch.rcIndices, m.contig.index.Lookup(sequence, -1)...)
	}

	log.Printf("ContigMatch.addPrimer %s %d %d\n%s\n%s\n", primer.sequence, primerMatch.indices, primerMatch.rcIndices, primer.sequences, primer.rcSequences)

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
				log.Printf("Start index %s %d/%d\n", contig.identifier, len(sequenceChan), cap(sequenceChan))
				start := time.Now()
				contig.index = suffixarray.New(contig.sequence)
				indexChan <- *contig
				log.Printf("End index %s %d %fs\n", contig.identifier, len(contig.sequence), time.Since(start).Seconds())
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
				log.Printf("Start match %s %d/%d\n", contig.identifier, len(indexChan), cap(indexChan))
				start := time.Now()

				contigMatch := ContigMatch{
					contig: contig,
				}

				log.Printf("Scan FORWARD primers in %s\n", contig.identifier)
				for _, primer := range primers.forward {
					contigMatch.addPrimer(primer, true)
				}
				log.Printf("Scan REVERSE primers in %s\n", contig.identifier)
				for _, primer := range primers.reverse {
					contigMatch.addPrimer(primer, false)
				}

				matchChan <- contigMatch

				log.Printf("End match %s %d %fs\n", contig.identifier, len(contig.sequence), time.Since(start).Seconds())
			}
			log.Printf("Shutdown match worker %d\n", i)
		}(matchChan, indexChan, primers, i)

	}

	wg.Wait()
	log.Println("Shutdown match WaitGroup")
}

type Primer struct {
	sequence    []byte
	sequences   [][]byte
	rcSequences [][]byte
}

func expandDegeneratePosition(primers [][]byte, position int, l ...byte) [][]byte {
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

func expandDegenerateSequence(sequence []byte) [][]byte {
	var primers [][]byte

	if len(sequence) == 0 {
		return nil
	}

	primer := make([]byte, len(sequence))
	copy(primer, sequence)

	primers = append(primers, primer)
	for i, nt := range sequence {
		switch nt {
		default:
			// TODO: should unrecognized characters panic?
			primers = expandDegeneratePosition(primers, i, nt)
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
	forward   []Primer
	rcReverse []Primer
	reverse   []Primer
	rcForward []Primer
}

func (p *PrimerList) Read(r io.Reader) error {
	var isForwardPrimer bool = true

	// TODO: validate primer character-set
	scanner := bufio.NewScanner(r)
	for scanner.Scan() {
		if len(scanner.Bytes()) == 0 {
			// It is assumed the primer list is a file with
			// the forward reads listed on separate lines followed by
			// a blank line followed by
			// the reverse reads listed on separate lines
			isForwardPrimer = !isForwardPrimer
			continue
		}
		sequence := make([]byte, len(scanner.Bytes()))
		copy(sequence, bytes.ToUpper(scanner.Bytes()))
		reverseComplement := reverseComplement(sequence)
		log.Printf("%s\n", reverseComplement)
		if isForwardPrimer {
			p.forward = append(p.forward, Primer{
				sequence:    sequence,
				sequences:   expandDegenerateSequence(sequence),
				rcSequences: expandDegenerateSequence(reverseComplement),
			})
			/*p.rcForward = append(p.rcForward, Primer{
				sequence:  sequence,
				sequences: expandDegenerateSequence(reverseComplement),
			})*/
			log.Printf("Add forward primer %s\n", sequence)
		} else {
			p.reverse = append(p.reverse, Primer{
				sequence:    sequence,
				sequences:   expandDegenerateSequence(sequence),
				rcSequences: expandDegenerateSequence(reverseComplement),
			})
			/*p.rcReverse = append(p.rcReverse, Primer{
				sequence:  sequence,
				sequences: expandDegenerateSequence(reverseComplement),
			})*/
			log.Printf("Add reverse primer %s\n", sequence)
		}
	}

	if err := scanner.Err(); err != nil {
		return err
	}

	// TODO: unit test entire file is read (especially the last line), all combinations are built, if raises an error if the character-set is invalid or both forward/reverse primers are not present (format error or empty file)
	/*if len(*p)%2 != 0 {
		return fmt.Errorf("The primer list file must contain alternating both forward and reverse primers")
	}*/

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

func reverseComplement(s []byte) []byte {
	sLen := len(s)
	rc := make([]byte, sLen)
	for i := 0; i < sLen; i++ {
		rc[i] = reverseComplementTable[s[sLen-i-1]]
	}
	return rc
}
