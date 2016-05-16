package fmi

import (
	"bytes"
	"log"
	"sort"
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
N/- any Nucleotide (not a gap)
*/
const NUCLEOTIDE_CODE = "ACGTUMRWSYKVHDBN-"

//const NUCLEOTIDE_CODE = "ABCDGHKMNRSTUVWY-"

type Index struct {
	data []byte
	sa   []int
}

func New(data []byte) *Index {
	// TODO: reenable special char
	//log.Println("reenable append special char")
	if data[len(data)-1] != '$' {
		data = append(data, '$')
	}
	return &Index{
		data: data,
		sa:   qsufsort(data),
	}
}

type indexQuery struct {
	// TODO: reduce mismatch to a smaller data type to conserve stack space. Be sure to pack struct.
	mismatch       int
	depth          int
	query          []byte
	sa             []int
	shortestSuffix int
}

/*func (q *indexQuery) minSuffixLength(i int) int {
}

func (q *indexQuery) needTrimSuffixArray() bool {
	if q.depth > q.minSuffixLength {
		// It doesn't matter what the current minSuffixLength is,
		// it is still greater than the current depth.
		return false
	}
	// The depth has exceeded the minSuffixLength since the last time minSuffixLength was checked.
	//if len(q.sa) < q.lastSaLength {
	// Scan the suffixarray.
	for i := range q.sa {
		if r := dataLength - sa[i] - 1; r < q.minSuffixLength {
			q.minSuffixLength = r
		}
	}
	// Did the minSuffixLength change since the last check?
	if q.depth > q.minSuffixLength {
		return false
	}
	//}
	return q.depth > q.minSuffixLength
}*/

/*type indexStack []indexQuery

func (s *indexStack) push(q indexQuery) {
	s = append(s, q)
}

func (s *indexQuery) pop() indexQuery {
	q := s[len(s)-1]
	s = s[:len(s)-1]
	return q
}*/

func (x *Index) Lookup(query []byte, mismatch int) []int {
	// Consider negative mismatch an exact match
	if mismatch < 0 {
		mismatch = 0
	}

	if mismatch == 0 {
		// TODO: replace with Index.exact() unless there is a measurable performance gain in this code duplication.
		// find matching suffix index range [i:j]
		// find the first index where s would be the prefix
		i := sort.Search(len(x.sa), func(i int) bool {
			return bytes.Compare(x.data[x.sa[i]:], query) >= 0
		})
		// starting at i, find the first index at which s is not a prefix
		j := i + sort.Search(len(x.sa)-i, func(j int) bool {
			return !bytes.HasPrefix(x.data[x.sa[j+i]:], query)
		})
		return x.sa[i:j]
	}

	// TODO: Does declaring these variables here instead of the start of the function have a measurable performance gain?
	var stack []indexQuery
	var locations []int

	// Initialize the stack with a indexQuery considering all positions.
	stack = append(stack, indexQuery{
		depth:    0,
		mismatch: mismatch,
		query:    query,
		sa:       x.sa,
	})
	for len(stack) > 0 {
		// stack.pop()
		iq := stack[len(stack)-1]
		stack = stack[:len(stack)-1]
		//log.Printf("%#v\n", iq)

		// TODO: can we optimize for single element remaining in the suffixarray?

		// TODO: enable
		/*if len(iq.sa) == 0 {
			// short-circuit - query not found
			continue
		}*/

		if len(iq.query) == 0 {
			// query found
			// Note: An empty initial query will match all positions.
			locations = append(locations, iq.sa...)
			continue
		}

		if iq.mismatch == 0 {
			// short-circuit, no need to compare character-by-character for an exact match.
			locations = append(locations, x.exact(iq)...)
			continue
		}

		// deletion
		// query:  gtgaatcatcraatytty
		// suffix: gt_aatcatcraatytty
		stack = append(stack, indexQuery{
			// Assuming there was as deletion, the query character is discarded
			// as it would not match anything in the suffixarray.
			query:    iq.query[1:],
			depth:    iq.depth,
			mismatch: iq.mismatch - 1,
			sa:       iq.sa,
		})

		// insertion
		// query:  gtgaatcatcraatytty
		// suffix: gtGgaatcatcraatytty
		// Cannot call x.bounds on the indexQuery because x.bounds consumes a query position,
		// thus we could not represent sequential insertions.
		//stack = append(stack, x.bounds(iq.query[0], indexQuery{
		var sa []int
		for i := range sa {
			if len(x.data)-iq.sa[i] > len(iq.query)+iq.depth+1-iq.mismatch-1 {
				sa = append(sa, iq.sa[i])
			} else {
				if len(x.data)-iq.sa[i]+iq.depth > 50 {
					log.Printf("%d %s... discard %d\n", len(iq.query), iq.query[:50], iq.sa[i])
				} else {
					log.Printf("%d %s discard %d\n", len(iq.query), iq.query, iq.sa[i])
				}
			}

		}
		stack = append(stack, indexQuery{
			query: iq.query,
			// Assuming there was an insertion, depth is increased to skip over the insertion
			depth:    iq.depth + 1,
			mismatch: iq.mismatch - 1,
			sa:       sa,
		})

		// match
		// query:  gtgaatcatcraatytty
		// suffix: gtgaatcatcraatytty
		stack = append(stack, x.bounds(iq.query[0], indexQuery{
			query:    iq.query,
			depth:    iq.depth,
			mismatch: iq.mismatch,
			sa:       iq.sa,
		}))

		// mismatch
		// query:  gtgaatcatcraatytty
		// suffix: Htgaatcatcraatytty
		// TODO: reduce the alphabet
		for _, nt := range "GATCN" {
			/*for _, nt := range NUCLEOTIDE_CODE {*/
			if byte(nt) == iq.query[0] {
				continue
			}
			stack = append(stack, x.bounds(byte(nt), indexQuery{
				query:    iq.query,
				depth:    iq.depth + 1,
				mismatch: iq.mismatch - 1,
				sa:       iq.sa,
			}))
		}
	}

	return locations
}

/*func (q *indexQuery) sweepForGaps() {
	var remaining int
	if q.minSuffixLength > len(iq.query)+iq.depth-iq.mismatch {
		return
	}
	dataLength := len(x.data)
	for i := range q.sa {
		remaining = dataLength - q.sa[i]
		if dataLength-q.sa[i] < q.minSuffixLength {
			q.minSuffixLength = remaining
		}
	}
	remaining = len(x.data) - 1 - idx
}*/

func (x *Index) bounds(nt byte, iq indexQuery) indexQuery {
	if len(iq.query) == 0 {
		// TODO: In practice, an indexQuery passed to this function should not have an empty query or zero mismatches.
		// This if block was added because a unit test for an empty query would panic with an index out of bounds.
		// Should panic be a valid response?
		return iq
	}

	// find matching suffix index range [i:j]
	// It is assumed the indexQuery suffixarray is sorted and all prefixes less than the current depth match

	// find the first index where s would be the prefix
	i := sort.Search(len(iq.sa), func(i int) bool {
		idx := iq.sa[i] + iq.depth
		/*defer func(i, idx int, iq indexQuery) {
			if err := recover(); err != nil {
				log.Fatalf("i:%d panic data:%d idx:%d %s %#v err:%s\n", i, len(x.data), idx, x.data, iq, err)
			}
		}(i, idx, iq)*/
		// The suffix, not including the special character, must be at least as long as the current query
		// The suffix must be at least as long as the entire query (allowing for mismatches)
		// The suffix must contain the current query
		//return len(x.data) > idx && len(x.data)-idx >= len(iq.query)-iq.mismatch && x.data[idx] >= nt
		//return len(x.data) > idx && len(x.data)-idx-1 > len(iq.query)-iq.mismatch && x.data[idx] >= nt
		//return len(x.data) > idx && x.data[idx] >= nt
		return len(x.data) > idx+len(iq.query)-iq.mismatch && x.data[idx] >= nt
		//return x.data[idx] >= nt
	})

	// starting at i, find the first index at which s is not a prefix
	j := i + sort.Search(len(iq.sa)-i, func(j int) bool {
		idx := iq.sa[j+i] + iq.depth
		defer func(iq indexQuery, i, j, idx int) {
			if err := recover(); err != nil {
				var sa []int
				for i := range iq.sa[i : j+i+1] {
					//if len(x.data)-iq.sa[i] > len(iq.query)+iq.depth+1-iq.mismatch {
					sa = append(sa, len(x.data[iq.sa[j+i]:]))
					//}
				}
				log.Fatalf("Index.bounds\nquery %s\nsuffix %s\n%v\nlen(data) %d len(sa) %d sa[i+j] %d depth %d i %d j %d idx %d\n%s", iq.query, x.data[iq.sa[j+i]:], sa, len(x.data), len(iq.sa), iq.sa[i+j], iq.depth, i, j, idx, err)
			}
		}(iq, i, j, idx)
		return x.data[idx] > nt
	})

	return indexQuery{
		mismatch: iq.mismatch,
		query:    iq.query[1:],
		sa:       iq.sa[i:j],
		depth:    iq.depth + 1,
	}
}

func (x *Index) exact(iq indexQuery) []int {
	// find matching suffix index range [i:j]
	// find the first index where s would be the prefix
	i := sort.Search(len(iq.sa), func(i int) bool {
		idx := iq.sa[i] + iq.depth
		return idx < len(x.data) && bytes.Compare(x.data[idx:], iq.query) >= 0
	})
	// starting at i, find the first index at which s is not a prefix
	j := i + sort.Search(len(iq.sa)-i, func(j int) bool {
		idx := iq.sa[j+i] + iq.depth
		defer func(iq indexQuery, i, idx int) {
			if err := recover(); err != nil {
				log.Fatalf("Index.bounds\n%v\nlen(data) %d len(sa) %d depth %d i %d j %d idx %d\n%s", iq.sa, len(x.data), len(iq.sa), iq.depth, i, j, idx, err)
			}
		}(iq, i, idx)
		return !bytes.HasPrefix(x.data[idx:], iq.query)
	})

	return iq.sa[i:j]
}

// Locate locates the pattern
/*func (i *Index) Locate(query []byte, mismatches int) ([]int, error) {
	var stack []int
	var locations []int
	//letters := byteutil.Alphabet(query)
	//for _, letter := range letters {
	for _, letter := range NUCELOTIDE_CODE {
		if _, ok := fmi.CountOfLetters[letter]; !ok {
			return locations, nil
		}
	}

	//if fmi.SuffixArray == nil {
	//	return nil, errors.New("SuffixArray is nil, you should call TransformForLocate instead of Transform")
	//}

	n := len(fmi.BWT)
	var matches stack.Stack
	type Match struct {
		query      []byte
		start, end int
		mismatches int
	}
	matches.Put(Match{query, 0, n - 1, mismatches})
	for !matches.Empty() {
		match := matches.Pop().(Match)
		query = match.query[0 : len(query)-1]
		last := match.query[len(query)-1]
		var letters []byte
		if mismatches == 0 {
			letters = []byte{last}
		} else {
			letters = fmi.Alphabet
		}
		for _, c := range letters {
			start := fmi.C[c] + fmi.Occ[c][match.start-2] + 1
			end := fmi.C[c] + fmi.Occ[c][match.end-1]
			if start <= end {
				if len(query) == 0 {
					for _, i := range fmi.SuffixArray[start : end+1] {
						locations = append(locations, i)
					}
				} else {
					mm := match.mismatches
					if c != last {
						if match.mismatches-1 > 0 {
							mm = match.mismatches - 1
						} else {
							mm = 0
						}
					}
					matches.Put(Match{query, start, end, mm})
				}
			}
		}
	}
	sort.Ints(locations)
	return locations, nil
}*/
