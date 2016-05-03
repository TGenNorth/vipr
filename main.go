package main

import (
	"bufio"
	"bytes"
	"flag"
	"fmt"
	//"github.com/pkg/profile"
	"index/suffixarray"
	"io"
	"log"
	"os"
	//"strings"
	"sync"
	"time"
)

type Contig struct {
	identifier string
	// NOTE: sequence MUST NOT be altered once assigned
	sequence []byte
	index    *suffixarray.Index
}

func NewContig(identifier string) *Contig {
	return &Contig{
		identifier: identifier,
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

type Match struct {
	contig         *Contig
	primer         *Primer
	forwardIndices []int
	reverseIndices []int
}

func (m *Match) Write(w io.Writer) {
	fmt.Fprintf(w, "%s %d %d\n", m.contig.identifier, m.forwardIndices, m.reverseIndices)
	for _, fIdx := range m.forwardIndices {
		for _, rIdx := range m.reverseIndices {
			if rIdx < fIdx || rIdx-fIdx > 500 {
				//if rIdx < fIdx {
				continue
			}
			//fmt.Printf("%s sequence %d fwd %s %d rev %s %d allele %d\n", m.contig.identifier, len(m.contig.sequence), m.primer.forward, fIdx, m.primer.reverse, rIdx, rIdx-fIdx)
			fmt.Fprintf(w, "%s\t%s\t%d\t%d\t%s\n", m.primer.forward, m.primer.reverse, rIdx, fIdx, m.contig.sequence[fIdx:rIdx+len(m.primer.reverse)])
			fmt.Fprintf(w, "%s\t%s\n", m.contig.sequence[fIdx:fIdx+len(m.primer.forward)], m.contig.sequence[rIdx:rIdx+len(m.primer.reverse)])
		}
	}
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

func init() {

}

func main() {
	flag.Parse()
	//defer profile.Start(profile.CPUProfile, profile.ProfilePath(".")).Stop()
	args := flag.Args()

	var threads int = 10
	// TODO: validate arguments
	sequenceFilename := args[0]
	primerListFilename := args[1]
	sequenceChan := make(chan *Contig, threads/2)
	indexChan := make(chan *Contig, threads)
	matchChan := make(chan []Match, threads)

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
	go func(indexChan chan<- *Contig, sequenceChan <-chan *Contig) {
		defer close(indexChan)
		var wg sync.WaitGroup
		wg.Add(threads / 2)
		for i := 0; i < threads/2; i++ {
			go func(indexChan chan<- *Contig, sequenceChan <-chan *Contig, i int) {
				defer wg.Done()
				for contig := range sequenceChan {
					log.Printf("Start index %s %d/%d\n", contig.identifier, len(sequenceChan), cap(sequenceChan))
					start := time.Now()
					contig.index = suffixarray.New(contig.sequence)
					indexChan <- contig
					log.Printf("End index %s %d %fs\n", contig.identifier, len(contig.sequence), time.Since(start).Seconds())
				}
				log.Printf("Shutdown index worker %d\n", i)
			}(indexChan, sequenceChan, i)
		}
		wg.Wait()
		log.Println("Shutdown index WaitGroup")
	}(indexChan, sequenceChan)

	go matchWorker(matchChan, indexChan, primers, threads)

	//matches := make(map[string]map[int]struct{})
	f := bufio.NewWriter(os.Stdout)
	defer f.Flush()
	for matches := range matchChan {
		for _, match := range matches {
			match.Write(f)
		}
	}
}

func matchWorker(matchChan chan []Match, indexChan chan *Contig, primers PrimerList, threads int) {
	defer close(matchChan)
	var wg sync.WaitGroup
	wg.Add(threads)
	for i := 0; i < threads; i++ {
		go func(matchChan chan<- []Match, indexChan <-chan *Contig, primers PrimerList, i int) {
			defer wg.Done()
			for contig := range indexChan {
				log.Printf("Start match %s %d/%d\n", contig.identifier, len(indexChan), cap(indexChan))
				start := time.Now()

				var forwardIndices, reverseIndices []int

				var matches []Match

				for i, primer := range primers {
					log.Printf("Lookup %s %s %s\n", contig.identifier, primer.forward, primer.reverse)
					for _, forwardPrimer := range primer.forwardExpanded {
						forwardIndices = append(forwardIndices, contig.index.Lookup(forwardPrimer, -1)...)
					}
					// No point in searching for reverse primers if the forward primer didn't match
					if len(forwardIndices) == 0 {
						forwardIndices = nil
						continue
					}
					for _, reversePrimer := range primer.reverseExpanded {
						reverseIndices = append(reverseIndices, contig.index.Lookup(reversePrimer, -1)...)
					}
					// No point in reporting match if the reverse primer didn't match
					if len(reverseIndices) == 0 {
						reverseIndices = nil
						continue
					}

					log.Printf("fwd: %s %s\n", primer.forward, bytes.Join(primer.forwardExpanded, []byte(",")))
					matches = append(matches, Match{
						contig:         contig,
						primer:         &primers[i],
						forwardIndices: forwardIndices,
						reverseIndices: reverseIndices,
					})
				}

				matchChan <- matches

				log.Printf("End match %s %d %fs\n", contig.identifier, len(contig.sequence), time.Since(start).Seconds())
			}
			log.Printf("Shutdown match worker %d\n", i)
		}(matchChan, indexChan, primers, i)

	}

	wg.Wait()
	log.Println("Shutdown match WaitGroup")
}

type Primer struct {
	forward         []byte
	reverse         []byte
	forwardExpanded [][]byte
	reverseExpanded [][]byte
}

func NewPrimer(forward, reverse []byte) Primer {
	p := Primer{
		forward: forward,
		reverse: reverse,
	}
	p.forwardExpanded = p.expandDegeneratePrimer(forward)
	p.reverseExpanded = p.expandDegeneratePrimer(reverse)
	return p
}

func (p Primer) expandDegeneratePosition(primers [][]byte, position int, l ...byte) [][]byte {
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

func (p Primer) expandDegeneratePrimer(sequence []byte) [][]byte {
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
			primers = p.expandDegeneratePosition(primers, i, nt)
		case 'A', 'C', 'G', 'T', 'U':
			primers = p.expandDegeneratePosition(primers, i, nt)
		case 'W':
			primers = p.expandDegeneratePosition(primers, i, 'A', 'T')
		case 'S':
			primers = p.expandDegeneratePosition(primers, i, 'G', 'C')
		case 'M':
			primers = p.expandDegeneratePosition(primers, i, 'A', 'C')
		case 'K':
			primers = p.expandDegeneratePosition(primers, i, 'G', 'T')
		case 'R':
			primers = p.expandDegeneratePosition(primers, i, 'A', 'G')
		case 'Y':
			primers = p.expandDegeneratePosition(primers, i, 'C', 'T')
		case 'B':
			primers = p.expandDegeneratePosition(primers, i, 'C', 'G', 'T')
		case 'D':
			primers = p.expandDegeneratePosition(primers, i, 'A', 'G', 'T')
		case 'H':
			primers = p.expandDegeneratePosition(primers, i, 'A', 'C', 'T')
		case 'V':
			primers = p.expandDegeneratePosition(primers, i, 'A', 'C', 'G')
		case 'N', '-':
			primers = p.expandDegeneratePosition(primers, i, 'G', 'A', 'T', 'C')
		}
	}
	return primers
}

type PrimerList []Primer

func (p *PrimerList) Read(r io.Reader) error {
	var isForwardPrimer bool = true
	var forwardPrimers, reversePrimers [][]byte

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
		line := make([]byte, len(scanner.Bytes()))
		copy(line, bytes.ToUpper(scanner.Bytes()))
		if isForwardPrimer {
			forwardPrimers = append(forwardPrimers, line)
			log.Printf("Add forward primer %s\n", line)
		} else {
			reversePrimers = append(reversePrimers, line)
			log.Printf("Add reverse primer %s\n", line)
		}
	}

	if err := scanner.Err(); err != nil {
		return err
	}

	log.Println("Building all forward/reverse primer combinations")
	for _, forwardPrimer := range forwardPrimers {
		for _, reversePrimer := range reversePrimers {
			log.Printf("%s\t%s\n", forwardPrimer, reversePrimer)
			*p = append(*p, NewPrimer(forwardPrimer, reversePrimer))
		}
	}

	// TODO: unit test entire file is read (especially the last line), all combinations are built, if raises an error if the character-set is invalid or both forward/reverse primers are not present (format error or empty file)
	/*if len(*p)%2 != 0 {
		return fmt.Errorf("The primer list file must contain alternating both forward and reverse primers")
	}*/

	return nil
}

/*
func reverseComplement(b []byte) []byte {
	rc := make([]byte, b)
	for _, c := range b {
		rc[i] = reverseComplementTable[c]
	}
}
*/
