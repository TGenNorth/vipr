package main

import (
	"bytes"
	"reflect"
	"testing"
)

func TestExpandDegenerateSequence(t *testing.T) {
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
	for _, tt := range testtable {
		primers := expandDegenerateSequence(tt.in)
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

/*func TestNewPrimer(t *testing.T) {
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

func TestReadFasta(t *testing.T) {
	descriptor := "ContigA"
	// 5000 was chosen as a value larger than the default buffer size
	s := make([]byte, 5000)
	for i := range s {
		s[i] = 'A'
	}
	copy(s[:], ">"+descriptor+"\n")
	ch := make(chan *Contig)
	go readFasta(ch, bytes.NewReader(s))

	contig, ok := <-ch

	if !ok {
		t.Fatalf("readFasta() channel was closed reading the first contig")
	}
	if contig.descriptor != descriptor {
		t.Errorf("readFasta() => descriptor %s, expected %s", contig.descriptor, descriptor)
	}
	if len(contig.sequence) != 5000-len(descriptor)-2 {
		t.Errorf("readFasta() => sequence length %d, expected %d", len(contig.sequence), 5000-len(descriptor)-2)
	}

	// FIXME: this test could lock
	contig, ok = <-ch
	if ok {
		t.Fatalf("readFasta() channel returned an unexpected contig: %s", contig)
	}
	/*select {
	case contig, ok := <-ch:
		if !ok {
			t.Fatalf("readFasta channel was closed reading the first contig")
		}
		if contig.descriptor != descriptor {
			t.Errorf("readFasta() => descriptor %s, expected %s", contig.descriptor, descriptor)
		}
		if len(contig.sequence) != 5000-len(descriptor)+2 {
			t.Errorf("readFasta() => sequence length %d, expected %d", len(contig.sequence), 5000-len(descriptor)+2)
		}
	default:
		//t.Fatalf("readFasta channel was blocked reading the first contig")
	}
	select {
	case contig, ok := <-ch:
		t.Errorf("%#v, %#v", contig, ok)
	}*/
}

func TestPrimerListRead(t *testing.T) {
	t.SkipNow()

	twoPrimerList := PrimerList{
		forward: []Primer{
			Primer{
				Label:    []byte("gatc"),
				Sequence: []byte("GATC"),
			},
		},
		reverse: []Primer{
			Primer{
				Label:    []byte("ctag"),
				Sequence: []byte("CTAG"),
			},
		},
	}

	labeledTwoPrimerList := PrimerList{
		forward: []Primer{
			Primer{
				Label:    []byte("ForwardPrimerLabel"),
				Sequence: []byte("GATC"),
			},
		},
		reverse: []Primer{
			Primer{
				Label:    []byte("ReversePrimerLabel"),
				Sequence: []byte("CTAG"),
			},
		},
	}

	var testtable = []struct {
		in  string
		out PrimerList
		err error
	}{
		// It should require at least one forward and one reverse primer sequence.
		// It should accept any sequence iff it is composed of the nucleic acid notation character set.
		// It should accept an optional primer descriptor in the UTF-8 character set following the primer sequence.
		// It should use the primer sequence as an descriptor if a primer descriptor is absent.
		// It should not modify the descriptor.
		{
			// empty file base case
			in:  "",
			out: PrimerList{},
			err: ErrInvalidFormat("at least one forward and one reverse primer sequence is required"),
		}, {
			in:  "gatc",
			out: PrimerList{},
			err: ErrInvalidFormat("at least one forward and one reverse primer sequence is required"),
		}, {
			// *nix line ending
			in:  "gatc\ngatc",
			out: PrimerList{},
			err: ErrInvalidFormat("at least one forward and one reverse primer sequence is required"),
		}, {
			// windows line ending
			in:  "gatc\r\ngatc",
			out: PrimerList{},
			err: ErrInvalidFormat("at least one forward and one reverse primer sequence is required"),
		},
		// It should require at least one blank line delimiter between the forward and reverse primer lists
		// It should t the number of delimiters
		{
			in:  "gatc\n\nctag",
			out: twoPrimerList,
			err: nil,
		}, {
			in:  "gatc\n\n\nctag",
			out: twoPrimerList,
			err: nil,
		}, {
			in:  "gAtC\n\ngatc\ngatc",
			out: PrimerList{},
			err: nil,
		}, {
			in:  "gatc ForwardPrimerLabel\n\nctag ReversePrimerLabel",
			out: labeledTwoPrimerList,
			err: nil,
		}, {
			in:  "gatc ForwardPrimerLabel\n\n\nctag ReversePrimerLabel",
			out: labeledTwoPrimerList,
			err: nil,
		}, {
			in:  "gatc ForwardPrimerLabel\n\n\n\nctag ReversePrimerLabel",
			out: labeledTwoPrimerList,
			err: nil,
		}, {
			in:  "gatc ForwardPrimerLabel\r\n\r\nctag ReversePrimerLabel",
			out: labeledTwoPrimerList,
			err: nil,
		}, {
			in:  "gatc Primer Label",
			out: PrimerList{},
			err: nil,
		}, {
			in:  "gatc Primer Label",
			out: PrimerList{},
			err: nil,
		}, {
			in:  "Primer Label gatc",
			out: PrimerList{},
			err: nil,
		},
	}

	for _, tt := range testtable {
		pl := PrimerList{}

		if err := pl.Read(bytes.NewReader([]byte(tt.in))); err != tt.err {
			t.Errorf("PrimerList.Read(%q) => %q, expected %q", tt.in, err, tt.err)
		}

		if !reflect.DeepEqual(pl, tt.out) {
			t.Errorf("PrimerList.Read(%q) => %q, expected %q\n", tt.in, pl, tt.out)
		}
	}
}

func TestReverseComplement(t *testing.T) {
	var testtable = []struct {
		in  []byte
		out []byte
		err error
	}{
		{[]byte(""), []byte(""), nil},
		{[]byte("GATC"), []byte("GATC"), nil},
		{[]byte("ACGTUMRWSYKVHDBN"), []byte("NVHDBMRSWYKAACGT"), nil},
		{[]byte("GKTARGTAATCCAACTAGCACCT"), []byte("AGGTGCTAGTTGGATTACYTAMC"), nil},
		//{[]byte("GATQ"), nil, ErrInvalidSequence("unrecognized nucleotide 'Q' at index 3 in sequence \"GATQ\"")},
		// FIXME: should it handle mixed-case?
		//{[]byte("gatq"), nil, ErrInvalidSequence("unrecognized nucleotide 'q' at index 3 in sequence GATQ")},
	}

	for _, tt := range testtable {
		//if rc, err := reverseComplement(tt.in); err != tt.err || !bytes.Equal(tt.out, rc) {
		if rc := reverseComplement(tt.in); !bytes.Equal(tt.out, rc) {
			t.Errorf("reverseComplement(%q) => %q expected %q\n", tt.in, rc, tt.out)
		}
	}
}
