package main

import (
	"bytes"
	"errors"
	"index/suffixarray"
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

func TestPrimerListSet(t *testing.T) {
	t.SkipNow()
	var testtable = []struct {
		in  string
		out PrimerList
		err error
	}{
		{"", PrimerList{}, errors.New("expected a filename or list of forward/reverse primers")},
		{",", PrimerList{}, errors.New("expected a colon delimiting the list of forward primers from the reverse primers")},
		{":", PrimerList{}, nil},
		{"#", PrimerList{}, nil},
		{"$", PrimerList{}, nil},
		{"invalidCharacters", PrimerList{}, nil},
		{"invalid:Characters", PrimerList{}, nil},
		{"GATC:invalidCharacters", PrimerList{}, nil},
		{"invalidCharacters:GATC", PrimerList{}, nil},
		{"GATC:", PrimerList{}, nil},
		{":GATC", PrimerList{}, nil},
		{"GATC", PrimerList{}, nil},
		{"GGTT,AACC", PrimerList{}, nil},
		{"GGTT,AACC:", PrimerList{}, nil},
		{":GGTT,AACC", PrimerList{}, nil},
		{":GGTT,AACC:", PrimerList{}, nil},
		{":GGTT:AACC:", PrimerList{}, nil},
		{"GGTT AACC", PrimerList{}, nil},
		{"GGTT\nAACC", PrimerList{}, nil},
	}

	for _, tt := range testtable {
		p := PrimerList{}
		if err := p.Set(tt.in); err != tt.err {
			t.Errorf("PrimerList.Set(%q) => \"%v\" expected \"%v\"\n", tt.in, err, tt.err)
		}
		if !reflect.DeepEqual(p, tt.out) {
			t.Errorf("PrimerList.Set(%q) %q expected %q\n", tt.in, p, tt.out)
		}
	}
}

func TestProbeSet(t *testing.T) {
	t.SkipNow()
	var testtable = []struct {
		in  string
		out ProbeFlag
		err error
	}{
		{"", ProbeFlag{}, nil},
		{",", ProbeFlag{}, nil},
		{":", ProbeFlag{}, nil},
		{"#", ProbeFlag{}, nil},
		{"$", ProbeFlag{}, nil},
		{"invalidCharacters", ProbeFlag{}, nil},
		{"invalid:Characters", ProbeFlag{}, nil},
		{"GATC:invalidCharacters", ProbeFlag{}, nil},
		{"invalidCharacters:GATC", ProbeFlag{}, nil},
		{"GATC:", ProbeFlag{}, nil},
		{":GATC", ProbeFlag{}, nil},
		{"GATC", ProbeFlag{}, nil},
		{"GGTT,AACC", ProbeFlag{}, nil},
		{"GGTT,AACC:", ProbeFlag{}, nil},
		{":GGTT,AACC", ProbeFlag{}, nil},
		{":GGTT,AACC:", ProbeFlag{}, nil},
		{":GGTT:AACC:", ProbeFlag{}, nil},
		{"GGTT AACC", ProbeFlag{}, nil},
		{"GGTT\nAACC", ProbeFlag{}, nil},
	}

	for _, tt := range testtable {
		p := ProbeFlag{}
		if err := p.Set(tt.in); err != tt.err {
			t.Errorf("ProbeFlag.Set(%q) => \"%v\" expected \"%v\"\n", tt.in, err, tt.err)
		}
		if !reflect.DeepEqual(p, tt.out) {
			t.Errorf("ProbeFlag.Set(%q) %q expected %q\n", tt.in, p, tt.out)
		}
	}
}

func TestContigMatchSequenceContainsProbe(t *testing.T) {
	sequence := []byte("ABCDEFGHIJKLMNOPQRSTUVWXYZ")
	contig := Contig{
		descriptor: "TestContig",
		sequence:   sequence,
		index:      suffixarray.New(sequence),
	}

	// index is a helper to allow letters to be used as an index in the alphabet.
	index := func(letter byte) int {
		return int(letter - 'A')
	}

	// seq is a helper to express a range as the first and last letter (inclusive)
	seq := func(a, b byte) []byte {
		return sequence[index(a) : index(b)+1]
	}

	var testtable = []struct {
		desc   string
		start  int
		end    int
		probe  []byte
		expect bool
	}{
		// The following test probes fully inside the sequence:
		{
			desc:   "it should accept the empty set of probes (complete contig)",
			start:  index('A'),
			end:    index('Z'),
			probe:  nil,
			expect: true,
		}, {
			desc:   "it should accept the empty set of probes (partial contig)",
			start:  index('D'),
			end:    index('W'),
			probe:  nil,
			expect: true,
		}, {
			desc: "it should accept a probe fully contained within the sequence (complete contig)",
			// In this variation of the fully contained probe, the start and end index of the sequence
			// include the entire contig. The result should not be affected by contig boundaries.
			start: index('A'),
			end:   index('Z'),
			// Probe start/end is any point within the sequence and not adjacent to a boundary.
			probe:  seq('L', 'P'),
			expect: true,
		}, {
			desc: "it should accept a probe fully contained within the sequence (partial contig)",
			// In this variation of the fully contained probe, the end index is moved
			// an arbitrary distance from the contig boundary.
			start: index('D'),
			end:   index('W'),
			// Probe start/end is any point within the sequence and not adjacent to a boundary.
			probe:  seq('L', 'P'),
			expect: true,
		}, {
			desc:   "it should accept a probe that starts at the sequence boundary and ends within the sequence",
			start:  index('A'),
			end:    index('C'),
			probe:  seq('A', 'B'),
			expect: true,
		}, {
			desc:   "it should accept a probe that starts within the sequence and ends at the sequence boundary",
			start:  'A' - 'A',
			end:    'C' - 'A',
			probe:  seq('B', 'C'),
			expect: true,
		}, {
			desc: "it should accept a probe that starts and ends at the sequence boundaries",
			// start/end is any valid index.
			start: 'A' - 'A',
			end:   'C' - 'A',
			// Probe start/end match the sequence boundaries.
			probe:  seq('A', 'C'),
			expect: true,
		},

		// The following test probes fully outside the sequence.
		{
			desc:  "it should not match a probe that starts/ends after the sequence (adjacent)",
			start: index('A'),
			end:   index('K'),
			// Probe start index is adjacent to the sequence boundary.
			// Probe end any valid index.
			probe:  seq('L', 'P'),
			expect: false,
		}, {
			desc:  "it should not match a probe that starts/ends after the sequence",
			start: index('A'),
			end:   index('K'),
			// Probe start index is any index after the sequence and not adjacent to the sequence boundary.
			// Probe end is any valid index.
			probe:  seq('N', 'P'),
			expect: false,
		}, {
			desc: "it should not match a probe that starts/ends before the sequence (adjacent)",
			// start/end is any valid index.
			start: index('L'),
			end:   index('P'),
			// Probe start is any valid index.
			// Probe end is adjacent to the sequence boundary.
			probe:  seq('A', 'K'),
			expect: false,
		}, {
			desc:  "it should not match a probe that starts/ends before the sequence",
			start: index('N'),
			end:   index('P'),
			// Probe start index is any index after the sequence and not adjacent to the sequence boundary.
			// Probe end is any valid index.
			probe:  seq('A', 'K'),
			expect: false,
		},

		// The following test probes that cross the sequence boundary.
		{
			// TODO: it should not accept a probe that starts before and ends at the start of the sequence.
			// TODO: it should not accept a probe that starts before and ends within the sequence.

			// The entire probe must be fully contained within the sequence.
			desc: "it should not match a probe that starts before and ends within the sequence",
			// start/end is any valid index
			start: index('B'),
			end:   index('D'),
			// Probe start is any valid index before the sequence.
			// Probe end is any index within the sequence.
			probe:  seq('A', 'C'),
			expect: false,
		}, {
			desc: "it should not match a probe that starts before and ends on the sequence start boundary",
			// start/end is any valid index
			start: index('N'),
			end:   index('P'),
			// Probe start index is any valid index before the sequence.
			// Probe end is any valid index.
			probe:  seq('A', 'N'),
			expect: false,
		}, {
			desc: "it should not match a probe that starts before and ends on the sequence end boundary",
			// start/end is any valid index
			start: index('N'),
			end:   index('P'),
			// Probe start index is any valid index before the sequence.
			// Probe end is any valid index.
			probe:  seq('A', 'P'),
			expect: false,
		}, {
			desc: "it should not match a probe that starts within the sequence and ends after",
			// start/end is any valid index
			start: index('N'),
			end:   index('P'),
			// Probe start index is any valid index before the sequence.
			// Probe end is any valid index.
			probe:  seq('O', 'Z'),
			expect: false,
		}, {
			desc: "it should not match a probe that starts on the sequence start boundary and ends after",
			// start/end is any valid index
			start: index('N'),
			end:   index('P'),
			// Probe start index is any valid index before the sequence.
			// Probe end is any valid index.
			probe:  seq('N', 'Z'),
			expect: false,
		}, {
			desc: "it should not match a probe that starts on the sequence end boundary and ends after",
			// start/end is any valid index
			start: index('N'),
			end:   index('P'),
			// Probe start index is any valid index before the sequence.
			// Probe end is any valid index.
			probe:  seq('P', 'Z'),
			expect: false,
		}, {
			// TODO: This should not happen unless the user selects an invalid min/max range.
			// Warn the user the range is too small for the sequences they are using.
			desc: "it should not match a probe that fully envelopes a sequence.",
			// start/end is any valid index.
			start: index('B'),
			end:   index('Y'),
			// probe is any range that envelops the start/end index.
			probe:  seq('A', 'Z'),
			expect: false,
		},
	}

	for i, tt := range testtable {
		defer func(i int) {
			if r := recover(); r != nil {
				t.Errorf("Test #%d sequenceContainsProbe(%q,%q,true) Probe: %q. PANIC: %s. %s.",
					i, sequence[tt.start], sequence[tt.end], tt.probe, r, tt.desc)
			}
		}(i)

		// Cannot use ProbeFlag.Set() because this test is using the an invalid Nucleotide character set
		// for the sequences which will cause the parser to return an error.
		probes := ProbeFlag{Probe{
			Sequence:    tt.probe,
			Sequences:   []text{tt.probe},
			RcSequences: []text{tt.probe},
		}}
		//if err := probe.Set(string(tt.probe)); err != nil {
		//	t.Fatalf("Test #%d sequenceContainsProbe(%q,%q,true) ProbeFlag.Set() returned an error: %s\n",
		//		i, sequence[tt.start], sequence[tt.end], err)
		//}

		m := NewContigMatch(contig, probes, PrimerList{})

		if result := m.sequenceContainsProbe(tt.start, tt.end, true); result != tt.expect {
			t.Errorf("Test #%d sequenceContainsProbe(%q,%q,true) %s. Probe: %q Result: %t Want: %t.",
				i, sequence[tt.start], sequence[tt.end], tt.desc, tt.probe, result, tt.expect)
		}
	}
}
