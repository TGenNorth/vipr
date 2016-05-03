package main

import (
	"bytes"
	"testing"
)

func TestExpandDegeneratePrimer(t *testing.T) {
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
		in  []string
		out []string
	}{
		{[]string{"TGTTAGGTAATCCAACTMGCACYT", "GATGGAGGACTCGTWYGCTTGT"}, []string{""}},
	}
	for _, tt := range testtable {
		primer := NewPrimer([]byte(tt.in[0]), []byte(tt.in[1]))
		if !bytes.Equal([]byte(tt.in[0]), primer.forward) {
			t.Errorf("NewPrimer(%s, %s) => primer.forward %s, expected %s\n", tt.in[0], tt.in[1], primer.forward, []byte(tt.in[0]))
		}
		if !bytes.Equal([]byte(tt.in[1]), primer.reverse) {
			t.Errorf("NewPrimer(%s, %s) => primer.forward %s, expected %s\n", tt.in[0], tt.in[1], primer.reverse, []byte(tt.in[1]))
		}
		/*if len(primers) != len(tt.out) {
			t.Fatalf("expandPrimerDegeneracies(%s) => %d primers %s, expected %d %s\n", tt.in, len(primers), primers, len(tt.out), tt.out)
		}
		for i := range tt.out {
			if !bytes.Equal(primers[i], tt.out[i]) {
				t.Errorf("expandPrimerDegeneracies(%s) => primer #%d %s, expected %s\n", tt.in, i, primers[i], tt.out[i])
			}
		}*/
	}
}
