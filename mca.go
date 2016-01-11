package main

import (
    "fmt"
    //"os"
    //"strings"
    "sort"
    "math"
)

type mcaType struct {
    seq_a []int32
    seq_b []int32
    swap bool
    nr int
    nc int
    cost func (x int, y int) int
    c []IntSl
    y []IntSl
}

const INF = math.MaxInt32

func Min(x, y int) int {
    if x < y {
        return x
    }
    return y
}

func Max(x, y int) int {
    if x > y {
        return x
    }
    return y
}

type IntSl struct {
    sl []int
}

func (self IntSl) Init(v int) {
    for i := range self.sl {
        self.sl[i] = v
    }
}

func (self IntSl) Max() (int, int) {
    max := - INF
    max_i := -1
    for i := range self.sl {
        if self.sl[i] > max {
            max = self.sl[i]
            max_i = i
        }
    }
    return max, max_i
}

func (self IntSl) BinaryIndexOf(v int) int {
    i := sort.SearchInts(self.sl, v)
    if i < len(self.sl) && self.sl[i] == v {
        return i
    }
    return -1
}

func (selfPtr *IntSl) extend(ext int) {
    (*selfPtr).sl = append((*selfPtr).sl, ext)
}

// ----------------------------------------------------------------------------------------
// --------------------------- cost matrix computation
// ----------------------------------------------------------------------------------------

func MakeMca(a string, b string, minMatch int) (mca mcaType) {
    mca.seq_a, mca.seq_b = []rune(a), []rune(b)
    mca.swap = len(mca.seq_a) > len(mca.seq_b)
    if(mca.swap) {
        mca.seq_a, mca.seq_b = mca.seq_b, mca.seq_a
    }

    mca.nr = len(mca.seq_a)
    mca.nc = len(mca.seq_b)

    mca.c = make([]IntSl, mca.nr)
    mca.y = make([]IntSl, mca.nr)

    mca.cost = func (x int, y int) int {
        i_y := mca.y[x].BinaryIndexOf(y)
        if i_y != -1 {
            return mca.c[x].sl[i_y]
        }
        return 0
    }

    // create index for bigrams of characters in b
    var y_index [256][256]IntSl
    for i := 0; i < len(mca.seq_b) -1; i++ {
        y_index[mca.seq_b[i] % 256][mca.seq_b[i+1] % 256].extend(i)
    }

    // Length of corresponding pieces in a and b is calculated. Efficient implementation by
    // doing lookahead for starting correspondences and grabbing those values otherwise

    c_cur_x := make([]int, mca.nc)  // materialized current row, 0 initialized
    c_prev_x := make([]int, mca.nc) // copy of last row, 0 initialized

    var y_indices []int          // positions matching current letter in a
    var prev_y_indices []int     // copy of last y_indices

    for x := range mca.seq_a {
        for pos := range prev_y_indices {
            c_prev_x[ prev_y_indices[pos] ] = 0
        }

        tmp_ref := c_prev_x
        c_prev_x = c_cur_x
        c_cur_x = tmp_ref

        prev_y_indices = y_indices;

        if x < len(mca.seq_a)-1 {
            y_indices = y_index[mca.seq_a[x] % 256][mca.seq_a[x+1] % 256].sl   // empty if character not in y_dic / b
        } else {
            y_indices = make([]int, 0)
        }

        for pos, prev_pos := 0, 0 ; pos < len(y_indices) || prev_pos < len(prev_y_indices) ; {
            cur_y_ind, prev_y_ind := INF, INF
            if pos < len(y_indices) { cur_y_ind = y_indices[pos] }
            if prev_pos < len(prev_y_indices) { prev_y_ind = prev_y_indices[prev_pos] }

            cxy := 0
            var y_merge int  // CHECK: will this re-initialize ?

            if prev_y_ind < cur_y_ind {
                // seq_a[x-1] == seq_b[prev_y_ind], seq_a[x] == seq_b[prev_y_ind+1]
                y_merge = prev_y_ind + 1
                prev_pos ++
                if prev_y_ind == cur_y_ind - 1 {
                    pos ++
                }
                if y_merge < len(mca.seq_b) && mca.seq_a[x] == mca.seq_b[y_merge] {
                    cxy = c_prev_x[prev_y_ind]
                }
            } else {
                // start new correspondence?
                y_merge = cur_y_ind
                pos ++
                bound := Min(mca.nr - x, mca.nc - y_merge)
                for cxy < bound && mca.seq_a[x+cxy] == mca.seq_b[y_merge+cxy] {
                    cxy++
                }
            }

            if cxy >= minMatch {    // on suboptimal matchings remove constraint for tests
                mca.c[x].extend(cxy)
                mca.y[x].extend(y_merge)
            }
            c_cur_x[y_merge] = cxy
        }
    }
    return mca
}

func main() {
    a := "ksajldkajsldkjalskd"
    b := "ksajldkajkkkkksldkjooooalsppkd"

    mca := MakeMca(a, b, 2)
    fmt.Printf("%+v\n", mca)
}
