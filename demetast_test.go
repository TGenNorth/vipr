package main

import (
	"bytes"
	"testing"
)

/*func TestExpandDegeneratePrimer(t *testing.T) {
	var testtable = []struct {
		in  []byte
		out [][]byte
	}{
		{[]byte(""), nil},
		{[]byte("GATC"), [][]byte{[]byte("GATC")}},
		{[]byte("WATC"), [][]byte{[]byte("AATC"), []byte("TATC")}},
		{[]byte("SATC"), [][]byte{[]byte("GATC"), []byte("CATC")}},
		{[]byte("DATC"), [][]byte{[]byte("AATC"), []byte("GATC"), []byte("TATC")}},
		{[]byte("N"), [][]byte{[]byte("G"), []byte("A"), []byte("T"), []byte("C")}},
		// It should expand multiple degenerarcies
		{[]byte("WWTC"), [][]byte{[]byte("AATC"), []byte("TATC"), []byte("ATTC"), []byte("TTTC")}},
		{[]byte("SATS"), [][]byte{[]byte("GATG"), []byte("CATG"), []byte("GATC"), []byte("CATC")}},
		{[]byte("TGTTAGGTAATCCAACTMGCACYT"), [][]byte{[]byte("TGTTAGGTAATCCAACTAGCACCT"), []byte("TGTTAGGTAATCCAACTCGCACCT"), []byte("TGTTAGGTAATCCAACTAGCACTT"), []byte("TGTTAGGTAATCCAACTCGCACTT")}},
	}
	p := Primer{}
	for _, tt := range testtable {
		primers := p.expandDegeneratePrimer(tt.in)
		if len(primers) != len(tt.out) {
			t.Fatalf("expandPrimerDegeneracies(%s) => %d primers %s, expected %d %s\n", tt.in, len(primers), primers, len(tt.out), tt.out)
		}
		for i := range tt.out {
			if !bytes.Equal(primers[i], tt.out[i]) {
				t.Errorf("expandPrimerDegeneracies(%s) => primer #%d %s, expected %s\n", tt.in, i, primers[i], tt.out[i])
			}
		}
	}
}

func TestNewPrimer(t *testing.T) {
	var testtable = []struct {
		forward         string
		reverse         string
		forwardExpanded []string
		reverseExpanded []string
	}{
		{
			forward:         "TGTTAGGTAATCCAACTMGCACYT",
			reverse:         "GATGGAGGACTCGTWYGCTTGT",
			forwardExpanded: []string{"TGTTAGGTAATCCAACTAGCACCT", "TGTTAGGTAATCCAACTCGCACCT", "TGTTAGGTAATCCAACTAGCACTT", "TGTTAGGTAATCCAACTCGCACTT"},
			reverseExpanded: []string{"GATGGAGGACTCGTACGCTTGT", "GATGGAGGACTCGTTCGCTTGT", "GATGGAGGACTCGTATGCTTGT", "GATGGAGGACTCGTTTGCTTGT"},
		},
	}
	for _, tt := range testtable {
		primer := NewPrimer([]byte(tt.forward), []byte(tt.reverse))
		if !bytes.Equal([]byte(tt.forward), primer.forward) {
			t.Errorf("NewPrimer(%s, %s) => primer.forward %s, expected %s\n", tt.forward, tt.reverse, primer.forward, tt.forward)
		}
		if !bytes.Equal([]byte(tt.reverse), primer.reverse) {
			t.Errorf("NewPrimer(%s, %s) => primer.reverse %s, expected %s\n", tt.forward, tt.reverse, primer.reverse, tt.reverse)
		}
		if len(primer.forwardExpanded) != len(tt.forwardExpanded) {
			t.Errorf("NewPrimer(%s, %s) => %d primer.forwardExpanded %s, expected %d %s\n", tt.forward, tt.reverse, len(primer.forwardExpanded), primer.forwardExpanded, len(tt.forwardExpanded), tt.forwardExpanded)
		}
		for i := range tt.forwardExpanded {
			if !bytes.Equal([]byte(tt.forwardExpanded[i]), primer.forwardExpanded[i]) {
				t.Errorf("NewPrimer(%s, %s) => #%d primer.forwardExpanded %s, expected %s\n", tt.forward, tt.reverse, i, primer.forwardExpanded[i], tt.forwardExpanded[i])
			}
		}
		if len(primer.reverseExpanded) != len(tt.reverseExpanded) {
			t.Errorf("NewPrimer(%s, %s) => %d primer.reverseExpanded %s, expected %d %s\n", tt.reverse, tt.reverse, len(primer.reverseExpanded), primer.reverseExpanded, len(tt.reverseExpanded), tt.reverseExpanded)
		}
		for i := range tt.reverseExpanded {
			if !bytes.Equal([]byte(tt.reverseExpanded[i]), primer.reverseExpanded[i]) {
				t.Errorf("NewPrimer(%s, %s) => #%d primer.reverseExpanded %s, expected %s\n", tt.reverse, tt.reverse, i, primer.reverseExpanded[i], tt.reverseExpanded[i])
			}
		}
	}
}*/

func TestMatchWorker(t *testing.T) {

}

func TestReverseComplement(t *testing.T) {
	var testtable = []struct {
		in  []byte
		out []byte
	}{
		{[]byte(""), []byte("")},
		{[]byte("GATC"), []byte("GATC")},
		{[]byte("ACGTUMRWSYKVHDBN"), []byte("NVHDBMRSWYKAACGT")},
		{[]byte("GKTARGTAATCCAACTAGCACCT"), []byte("AGGTGCTAGTTGGATTACYTAMC")},
	}

	for _, tt := range testtable {
		rc := reverseComplement(tt.in)
		if !bytes.Equal(tt.out, rc) {
			t.Errorf("reverseComplement(%s) => %s, expected %s\n", tt.in, rc, tt.out)
		}
	}
}

func TestReadFasta(t *testing.T) {
	identifier := "Hello World"
	// 5000 was chosen as a value larger than the default buffer size
	s := make([]byte, 5000)
	for i := range s {
		s[i] = 'A'
	}
	copy(s[:], ">"+identifier+"\n")
	ch := make(chan *Contig)
	go readFasta(ch, bytes.NewReader(s))
	select {
	case contig, ok := <-ch:
		if !ok {
			t.Fatalf("readFasta channel was closed reading the first contig")
		}
		if contig.identifier != identifier {
			t.Errorf("readFasta() => identifier %s, expected %s", contig.identifier)
		}
		if len(contig.sequence) != 5000-len(identifier)+2 {
			t.Errorf("readFasta() => sequence length %d, expected %d", len(contig.sequence), 5000-len(identifier)+2)
		}
	default:
		t.Fatalf("readFasta channel was blocked reading the first contig")
	}
	select {
	case contig, ok := <-ch:
		t.Errorf("%#v, %#v", contig, ok)
	}
}
