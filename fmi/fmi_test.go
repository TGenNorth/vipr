package fmi

import (
	"index/suffixarray"
	"io/ioutil"
	"log"
	"os"
	"reflect"
	"runtime/debug"
	"testing"
)

const PRIMER = "AACTTTYRRCAAYGGATCWCT"

func TestNew(t *testing.T) {
	var testtable = []struct {
		in  []byte
		out Index
	}{
		{
			// It should build an empty index given an empty input
			in:  nil,
			out: Index{},
		}, {
			in:  []byte(""),
			out: Index{},
		}, {
			// It should should append a special character at the end
			// if the user did not provide one.
			in: []byte("banana"),
			out: Index{
				data: []byte("banana$"),
				sa:   []int{6, 5, 3, 1, 0, 4, 2},
			},
		}, {
			// It should build a suffixarray
			in: []byte("banana$"),
			out: Index{
				data: []byte("banana$"),
				sa:   []int{6, 5, 3, 1, 0, 4, 2},
			},
		},
	}
	for i, tt := range testtable {
		defer func(i int) {
			if err := recover(); err != nil {
				t.Errorf("panic on test case %d\n%s", i, debug.Stack())
			}
		}(i)
		idx := New(tt.in)
		if !reflect.DeepEqual(*idx, tt.out) {
			t.Errorf("New(%q) => %#v, expected %#v", tt.in, idx, tt.out)
		}
	}
}

func TestIndexLookup(t *testing.T) {
	index := New([]byte("banana$"))
	var testtable = []struct {
		index    *Index
		query    []byte
		mismatch int
		out      []int
	}{
		{
			// It should find all positions given an empty query
			index:    index,
			query:    []byte(""),
			mismatch: 0,
			out:      []int{6, 5, 3, 1, 0, 4, 2},
		}, {
			// It should find all exact matches
			index:    index,
			query:    []byte("banana"),
			mismatch: 0,
			out:      []int{0},
		}, {
			index:    index,
			query:    []byte("na"),
			mismatch: 0,
			out:      []int{4, 2},
		}, {
			index:    index,
			query:    []byte("a"),
			mismatch: 0,
			out:      []int{5, 3, 1},
		}, {
			// It should still find exact values when allowing for mismatches
			index:    index,
			query:    []byte(""),
			mismatch: 1,
			out:      []int{6, 5, 3, 1, 0, 4, 2},
		}, {
			index:    index,
			query:    []byte("banana"),
			mismatch: 1,
			out:      []int{1},
		}, {
			index:    index,
			query:    []byte("na"),
			mismatch: 1,
			out:      []int{4, 2, 4, 2, 3, 1, 5, 3, 1},
		}, {
			index:    index,
			query:    []byte("a"),
			mismatch: 1,
			out:      []int{5, 3, 1, 3, 1, 0, 4, 2, 6, 5, 3, 1, 0, 4, 2},
		}, {
			// Assuming an insertion of the middle 'a', the query anna should match the substring anana
			index:    index,
			query:    []byte("anna"),
			mismatch: 1,
			// deletion:  a_na  1,3
			// deletion:  an_a  1,3
			// insertion: anAna 1
			out: []int{1, 3, 1, 3, 1},
		}, {
			index:    index,
			query:    []byte("banana"),
			mismatch: 2,
			out:      []int{1},
		}, {
			// TODO: test short suffix in the middle of a suffix range (allowed by insertions) creates gap in range
			// Assuming an insertion of the middle 'a', the query anna should match the substring anana
			/*index:    index,
				query:    []byte("a"),
				mismatch: 3,
				out:      []int{1},
			}, {*/
			// It should behave the same for negative mismatches as exact matches (zero)
			index:    index,
			query:    []byte(""),
			mismatch: -1,
			out:      []int{6, 5, 3, 1, 0, 4, 2},
		}, {
			index:    index,
			query:    []byte("a"),
			mismatch: -1,
			out:      []int{5, 3, 1},
		}, {
			index:    index,
			query:    []byte("a"),
			mismatch: -1,
			out:      []int{5, 3, 1},
		}, {
			// It should not match unless the query is present
			index:    index,
			query:    []byte("x"),
			mismatch: 0,
			out:      []int{},
		}, {
			// It should not match partial queries
			index:    index,
			query:    []byte("bananaa"),
			mismatch: 0,
			out:      []int{},
		}, {
			index:    index,
			query:    []byte("nanaa"),
			mismatch: 0,
			out:      []int{},
		},
	}
	for i, tt := range testtable {
		defer func(i int) {
			if err := recover(); err != nil {
				t.Errorf("panic on test case %d\n%s", i, debug.Stack())
			}
		}(i)
		// TODO: is the result sorted?
		// TODO: is the result unique?
		result := tt.index.Lookup(tt.query, tt.mismatch)
		if !reflect.DeepEqual(result, tt.out) {
			t.Errorf("Test #%d index.Lookup(%q, %d) => %#v, expected %#v", i, tt.query, tt.mismatch, result, tt.out)
		}
	}
}

func TestIndexBoundsInsertAnnaShouldMatchAnana(t *testing.T) {
	index := New([]byte("banana$"))
	var testtable = []struct {
		nt  byte
		iq  indexQuery
		out indexQuery
	}{
		{
			nt: 'a',
			iq: indexQuery{
				mismatch: 1,
				depth:    0,
				query:    []byte("anna"),
				sa:       []int{6, 5, 3, 1, 0, 4, 2},
			},
			out: indexQuery{
				mismatch: 1,
				depth:    1,
				query:    []byte("nna"),
				sa:       []int{3, 1},
			},
		}, {
			nt: 'n',
			iq: indexQuery{
				mismatch: 1,
				depth:    1,
				query:    []byte("nna"),
				sa:       []int{3, 1},
			},
			out: indexQuery{
				mismatch: 1,
				depth:    2,
				query:    []byte("na"),
				sa:       []int{3, 1},
			},
		}, {
			nt: 'n',
			iq: indexQuery{
				mismatch: 0,
				depth:    3,
				query:    []byte("na"),
				sa:       []int{3, 1},
			},
			out: indexQuery{
				mismatch: 0,
				depth:    4,
				query:    []byte("a"),
				sa:       []int{1},
			},
		}, {
			nt: 'a',
			iq: indexQuery{
				mismatch: 0,
				depth:    4,
				query:    []byte("a"),
				sa:       []int{1},
			},
			out: indexQuery{
				mismatch: 0,
				depth:    5,
				query:    []byte(""),
				sa:       []int{1},
			},
		},
	}

	test := func(i int, tt struct {
		nt  byte
		iq  indexQuery
		out indexQuery
	}) {
		defer func(i int) {
			if err := recover(); err != nil {
				t.Errorf("panic on test case %d\n%s", i, debug.Stack())
			}
		}(i)
		result := index.bounds(tt.nt, tt.iq)
		if !reflect.DeepEqual(result, tt.out) {
			t.Errorf("Test #%d index.bounds(%q, %#v) => %#v, expected %#v", i, tt.nt, tt.iq, result, tt.out)
		}
	}

	for i, tt := range testtable {
		test(i, tt)
	}
}

func TestIndexBoundsDeleteAnaanaShouldMatchAnana(t *testing.T) {
	index := New([]byte("banana$"))
	var testtable = []struct {
		nt  byte
		iq  indexQuery
		out indexQuery
	}{
		{
			nt: 'a',
			iq: indexQuery{
				mismatch: 1,
				depth:    0,
				query:    []byte("anaana"),
				sa:       []int{6, 5, 3, 1, 0, 4, 2},
			},
			out: indexQuery{
				mismatch: 1,
				depth:    1,
				query:    []byte("naana"),
				sa:       []int{1},
			},
		}, {
			nt: 'n',
			iq: indexQuery{
				mismatch: 1,
				depth:    1,
				query:    []byte("naana"),
				sa:       []int{1},
			},
			out: indexQuery{
				mismatch: 1,
				depth:    2,
				query:    []byte("aana"),
				sa:       []int{1},
			},
		}, {
			nt: 'a',
			iq: indexQuery{
				mismatch: 0,
				depth:    2,
				query:    []byte("ana"),
				sa:       []int{1},
			},
			out: indexQuery{
				mismatch: 0,
				depth:    3,
				query:    []byte("na"),
				sa:       []int{1},
			},
		}, {
			nt: 'n',
			iq: indexQuery{
				mismatch: 0,
				depth:    3,
				query:    []byte("na"),
				sa:       []int{1},
			},
			out: indexQuery{
				mismatch: 0,
				depth:    4,
				query:    []byte("a"),
				sa:       []int{1},
			},
		},
	}

	test := func(i int, tt struct {
		nt  byte
		iq  indexQuery
		out indexQuery
	}) {
		defer func(i int) {
			if err := recover(); err != nil {
				t.Errorf("panic on test case %d\n%s", i, debug.Stack())
			}
		}(i)
		result := index.bounds(tt.nt, tt.iq)
		if !reflect.DeepEqual(result, tt.out) {
			t.Errorf("Test #%d index.bounds(%q, %#v) => %#v, expected %#v", i, tt.nt, tt.iq, result, tt.out)
		}
	}

	for i, tt := range testtable {
		test(i, tt)
	}
}

func TestIndexBounds(t *testing.T) {
	index := New([]byte("banana$"))
	var testtable = []struct {
		nt  byte
		iq  indexQuery
		out indexQuery
	}{
		// TODO: test mismatch exceeds query length
		// TODO: empty input struct
		{ //}, {
			nt: 'a',
			iq: indexQuery{
				mismatch: 0,
				depth:    0,
				query:    []byte(""),
				sa:       []int{},
			},
			out: indexQuery{
				mismatch: 0,
				depth:    0,
				query:    []byte(""),
				sa:       []int{},
			},
		}, {
			nt: 'a',
			iq: indexQuery{
				mismatch: 0,
				depth:    0,
				query:    []byte("ana"),
				sa:       []int{6, 5, 3, 1, 0, 4, 2},
			},
			out: indexQuery{
				mismatch: 0,
				depth:    1,
				query:    []byte("na"),
				sa:       []int{3, 1},
			},
		}, {
			nt: 'n',
			iq: indexQuery{
				mismatch: 1,
				depth:    1,
				query:    []byte("na"),
				sa:       []int{3, 1},
			},
			out: indexQuery{
				mismatch: 1,
				depth:    2,
				query:    []byte("a"),
				sa:       []int{3, 1},
			},
		},
	}

	test := func(i int, tt struct {
		nt  byte
		iq  indexQuery
		out indexQuery
	}) {
		defer func(i int) {
			if err := recover(); err != nil {
				t.Errorf("panic on test case %d\n%s", i, debug.Stack())
			}
		}(i)
		result := index.bounds(tt.nt, tt.iq)
		if !reflect.DeepEqual(result, tt.out) {
			t.Errorf("Test #%d index.bounds(%q, %#v) => %#v, expected %#v", i, tt.nt, tt.iq, result, tt.out)
		}
	}

	for i, tt := range testtable {
		test(i, tt)
	}
}

var benchmarkIndex *Index
var benchmarkContig []byte

func NewBenchmarkContig() []byte {
	if benchmarkContig != nil {
		return benchmarkContig
	}
	file, err := os.Open("testdata/contig.go")
	if err != nil {
		log.Fatal(err)
	}
	defer file.Close()
	data, err := ioutil.ReadFile("testdata/contig.go")
	if err != nil {
		log.Fatal(err)
	}
	benchmarkContig = data
	return benchmarkContig
}

func NewBenchmarkIndex() *Index {
	if benchmarkIndex != nil {
		return benchmarkIndex
	}

	benchmarkIndex = New(NewBenchmarkContig())
	return benchmarkIndex
}

func benchmarkSuffixarraySearch(n int, b *testing.B) {
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		suffixarray.New([]byte("banana"))
	}
}

func BenchmarkBuild(b *testing.B) {
	for i := 0; i < b.N; i++ {
		New(NewBenchmarkContig())
	}
}

func BenchmarkBuildSuffixarray(b *testing.B) {
	for i := 0; i < b.N; i++ {
		suffixarray.New(NewBenchmarkContig())
	}
}

func BenchmarkSearchSuffixarray(b *testing.B) {
	index := suffixarray.New(NewBenchmarkContig())
	primer := NewBenchmarkContig()[:22]
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		index.Lookup([]byte(primer), -1)
	}
}

func BenchmarkSearch0(b *testing.B) {
	index := NewBenchmarkIndex()
	primer := NewBenchmarkContig()[:22]
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		index.Lookup([]byte(primer), -1)
	}
}

func BenchmarkSearch1(b *testing.B) {
	index := NewBenchmarkIndex()
	primer := NewBenchmarkContig()[:22]
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		index.Lookup([]byte(primer), 1)
	}
}

func BenchmarkSearch2(b *testing.B) {
	index := NewBenchmarkIndex()
	primer := NewBenchmarkContig()[:22]
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		index.Lookup([]byte(primer), 2)
	}
}

func BenchmarkSearch3(b *testing.B) {
	index := NewBenchmarkIndex()
	primer := NewBenchmarkContig()[:22]
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		index.Lookup([]byte(primer), 3)
	}
}
