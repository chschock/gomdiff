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

func Init(sl []int, v int) {
    for i := range sl {
        sl[i] = v
    }
}

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

func (self.IntSl) Len() {
    return len(self.sl)
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

type solType struct {
    xy []int
    yx []int
    cost int
}

func max_cost_assignment(mca) solType {
    var lx, ly = make([]int, mca.nr), make([]int, mca.nc+mca.nr)
    var xy, yx = make([]int, mca.nr), make([]int, mca.nc+mca.nr)
    var S, T = make([]bool, mca.nr), make([]bool, mca.nc+mca.nr)
    var slack = make([]int, mca.nc+mca.nr)
    var slackx = make([]int, mca.nc+mca.nr)
    var rev_path = make([]int, mca.nr);

    // debug
    var mark_pos = make([]int, mca.nr)  // 0 initialized

    var y_merge = make([]int, mca.nc+mca.nr), prev_y_merge = make([]int, mca.nc+mca.nr)
    var y_merge_back = 0, prev_y_merge_back

    merge_relevant_columns := func (x int) {
        if q_back > 1 {
            y_merge, prev_y_merge = prev_y_merge, y_merge
            prev_y_merge_back = y_merge_back
            y_merge_back = 0

            y_next := INF;
            if mca.y[x].len() > 0 { y_next = Min(y_next, mca.y[x].sl[0]) }
            if prev_y_merge_back > 0 { y_next = Min(y_next, prev_y_merge[0]) }

            var i_y = 0, prev_i_y = 0
            for var y_cur; y_next < INF ; {
                y_cur, y_next = y_next, INF
                y_merge[y_merge_back++] = y_cur

                for i_y < mca.y[x].len() && mca.y[x].sl[i_y] <= y_cur { i_y ++ }
                if  i_y < mca.y[x].len() {
                    y_next = Min(y_next, mca.y[x].sl[i_y])
                }

                for prev_i_y < prev_y_merge_back && prev_y_merge[prev_i_y] <= y_cur { prev_i_y ++ }
                if  prev_i_y < prev_y_merge_back {
                    y_next = Min(y_next, prev_y_merge[prev_i_y])
                }
            }

        } else {
            y_merge_back = 0
            for i_y := 0; i_y < mca.y[x].len(); i_y++ {
                y_merge[y_merge_back++] = mca.y[x].sl[i_y]
            }
        }
    }

    add_to_tree := func (x int, greedy_break bool)
    {
        if xy[x] != -1 {   // not for root of path
            T[xy[x]] = true
        }

        q[q_back++] = x
        S[x] = true

        // update slack and slackx
        for i_y := 0; i_y < mca.y[x].len(); i_y ++ {
            y := mca.y[x].sl[i_y]
            cxy := mca.c[x].sl[i_y]

            if lx[x] + ly[y] - cxy < slack[y] {
                slack[y] = lx[x] + ly[y] - cxy
                slackx[y] = x
            }

            if greedy_break && slack[y] == 0 && !T[y] {
                if slack0_pivot == -1 { slack0_pivot = y }
                // condition leads immediately to a match in grow_eq_subgraph
                if yx[y] == -1 {
                    slack0_pivot = y
                    break
                }
            }
        }
        // check_slack(slack, slackx, mat, lx, ly);

        merge_relevant_columns(x)
    }

    augment := func (x int, y int) {
        x_aug, y_aug = x, y
        found_augmentation = true
    }

    // try to grow a path in the equality subgraph
    explore_eq_subgraph := func () {
        if debug { console.log('explore equality subgraph') }

        for q_back > q_front && !found_augmentation {
            x = q[q_front++]

            // plausible heuristic: first try x_aug+1 with y_aug+1
            if y_aug {
                var y = y_aug + 1;
                if (!T[y] &&  yx[y] == -1 && mat.cost(x, y) == lx[x] + ly[y])
                {
                    augment(x, y);
                    stat_upper++; // stat
                    continue;
                }
            }

            // sparse heuristic: try only y with positive cost
            for (var i_y = 0; i_y < mat.y[x].length; ++i_y)
            {
                var y = mat.y[x][i_y];
                if (mat.c[x][i_y] == lx[x] + ly[y] && !T[y])
                {
                    if (yx[y] == -1)  // if y is unmatched
                    {
                        augment(x, y);
                        stat_mid++; // stat
                        break;
                    }

                    rev_path[yx[y]] = x;
                    add_to_tree(yx[y], false);
                }
            }
        }
    }

    // improve labels lx, ly to grow the equality subgraph
    function improve_labelling()
    {
        if (debug) console.log('improve labelling');
        stat_improve++; // stat

        var delta = MSI;

        // for (var y=0; y<mat.nc+mat.nr; y++)
        //     if (slack[y] < MSI && y_merge.indexOf(y) == -1) // && mat.cost(slackx[y], y) > 0)
        //         console.log(y_merge + ' y ' + y + ' slack[y] ' + slack[y] + ' cost ' + mat.cost(slackx[y], y));

        // for (var y = 0; y < mat.nc+mat.nr; y++) {
        for (var i_y = 0; i_y < y_merge_back; i_y++) {
            var y = y_merge[i_y];
            if (!T[y] && slack[y] < delta)
                delta = slack[y];
        }

        if (debug && delta > 1000000000) alert('no minimal delta found');

        // adapt node labels: q holds all nodes in S
        lx[q[0]] -= delta;
        for (var i = 1; i < q_back; i++) {
            var x = q[i];
            lx[x] -= delta;
            ly[xy[x]] += delta;
        }

        slack0_pivot = -1;
        // for (var y = 0; y < mat.nc+mat.nr; y++) {
        for (var i_y = 0; i_y < y_merge_back; i_y++) {
            var y = y_merge[i_y];
            if (!T[y] && slack[y] < MSI) {
                slack[y] -= delta;
                if (slack[y] == 0) {
                    // store the first with ![T[y] && slack[y]==0, it exists
                    if (slack0_pivot == -1) slack0_pivot = y;
                    // condition leads immediately to a match in grow_eq_subgraph
                    if (yx[y] == -1) {
                        slack0_pivot = y;
                        break;
                    }
                }
            }
        }
    }

    // grow path in the equality subgraph
    function grow_eq_subgraph()
    {
        if debug { console.log('grow equality subgraph') }
        // search all y: this ensures path augmentation
        // for (var y = slack0_pivot; y < mat.nc; ++y)
        while (slack0_pivot != -1)
        {
            var y = slack0_pivot;
            slack0_pivot = -1;

            // if (! (!T[y] && slack[y] == 0)) break;
            if (yx[y] == -1)    // if y is unmatched
            {
                augment(slackx[y], y);
                stat_lower++; // stat
                break;
            }
            else
            {
                stat_marks++; // stat
                mark_pos[match_cnt] ++;
                // if (!S[yx[y]])     // never false because: S[xy[y]] => T[y]
                rev_path[yx[y]] = slackx[y];
                add_to_tree(yx[y], true);
                // else T[y] = true; // is true anyways, as yx[y] cannot be the path root
            }
        }
    }

}

func main() {
    a := "ksajldkajsldkjalskd"
    b := "ksajldkajkkkkksldkjooooalsppkd"

    mca := MakeMca(a, b, 2)
    fmt.Printf("%+v\n", mca)
}
