package main

import (
    "fmt"
    "sort"
    "math"
)

const INF = math.MaxInt32

var debug = false
var sort_by_cost = true
var skip_preproc = false
var early_stop = true

type Fp float32

type McaType struct {
    SeqA []int32    `json:"seq_a"`
    SeqB []int32    `json:"seq_b"`
    Swap bool       `json:"swap"`
    NR int          `json:"nr"`
    NC int          `json:"nc"`
    C []IntSl       `json:"c"`
    Y []IntSl       `json:"y"`
    cost func (x int, y int) int
    minMatch int
}

type SolType struct {
    Cost int    `json:"cost"`
    XY []int    `json:"xy"`
    YX []int    `json:"yx"`
}

func (mca McaType) String() string {
    return fmt.Sprintf("%v\n%v\nnr: %d, nc: %d, swapped: %v\nc: %v\ny: %v\n",
        mca.SeqA, mca.SeqB, mca.NR, mca.NC, mca.Swap, mca.C, mca.Y)
}

func Init(sl []int, v int) {
    for i := range sl {
        sl[i] = v
    }
}

func InitBool(sl []bool, v bool) {
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

// int slice with augmented functionality

type IntSl []int

func (self IntSl) Max() (int, int) {
    max := - INF
    max_i := -1
    for i := range self {
        if self[i] > max {
            max = self[i]
            max_i = i
        }
    }
    return max, max_i
}

func (self IntSl) BinaryIndexOf(v int) int {
    i := sort.SearchInts(self, v)
    if i < len(self) && self[i] == v {
        return i
    }
    return -1
}

func (selfPtr *IntSl) extend(ext int) {
    *selfPtr = append(*selfPtr, ext)
}

// aggregate int slice with custom order function

type CompareFunType func(i, j int) bool

type PriorityIntSl struct {
    IntSl
    less CompareFunType
}

func (self PriorityIntSl) Len() int { return len(self.IntSl) }
func (p PriorityIntSl) Less(i, j int) bool { return p.less(i, j) }
func (p PriorityIntSl) Swap(i, j int) { p.IntSl[i], p.IntSl[j] = p.IntSl[j], p.IntSl[i] }

// ----------------------------------------------------------------------------------------
// --------------------------- cost matrix computation
// ----------------------------------------------------------------------------------------

func MakeMca(a string, b string, minMatch int) (mca McaType) {
    mca.SeqA, mca.SeqB = []rune(a), []rune(b)
    mca.Swap = len(mca.SeqA) > len(mca.SeqB)
    if(mca.Swap) {
        mca.SeqA, mca.SeqB = mca.SeqB, mca.SeqA
    }

    mca.NR = len(mca.SeqA)
    mca.NC = len(mca.SeqB)

    mca.C = make([]IntSl, mca.NR)
    mca.Y = make([]IntSl, mca.NR)

    mca.cost = func (x int, y int) int {
        i_y := mca.Y[x].BinaryIndexOf(y)
        if i_y != -1 {
            return mca.C[x][i_y]
        }
        return 0
    }
    mca.minMatch = minMatch

    // create index for bigrams of characters in b
    var y_index [256][256]IntSl
    for i := 0; i < len(mca.SeqB) -1; i ++ {
        y_index[mca.SeqB[i] % 256][mca.SeqB[i+1] % 256].extend(i)
    }

    // Length of corresponding pieces in a and b is calculated. Efficient implementation by
    // doing lookahead for starting correspondences and grabbing those values otherwise

    c_cur_x := make([]int, mca.NC)  // materialized current row, 0 initialized
    c_prev_x := make([]int, mca.NC) // copy of last row, 0 initialized

    var y_indices []int          // positions matching current letter in a
    var prev_y_indices []int     // copy of last y_indices

    for x := range mca.SeqA {
        for pos := range prev_y_indices {
            c_prev_x[ prev_y_indices[pos] ] = 0
        }

        tmp_ref := c_prev_x
        c_prev_x = c_cur_x
        c_cur_x = tmp_ref

        prev_y_indices = y_indices

        if x < len(mca.SeqA)-1 {
            y_indices = y_index[mca.SeqA[x] % 256][mca.SeqA[x+1] % 256]   // empty if character not in y_dic / b
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
                // SeqA[x-1] == SeqB[prev_y_ind], SeqA[x] == SeqB[prev_y_ind+1]
                y_merge = prev_y_ind + 1
                prev_pos ++
                if prev_y_ind == cur_y_ind - 1 {
                    pos ++
                }
                if y_merge < len(mca.SeqB) && mca.SeqA[x] == mca.SeqB[y_merge] {
                    cxy = c_prev_x[prev_y_ind]
                }
            } else {
                // start new correspondence?
                y_merge = cur_y_ind
                pos ++
                bound := Min(mca.NR - x, mca.NC - y_merge)
                for cxy < bound && mca.SeqA[x+cxy] == mca.SeqB[y_merge+cxy] {
                    cxy ++
                }
            }

            if cxy >= minMatch {    // on suboptimal matchings remove constraint for tests
                mca.C[x].extend(cxy)
                mca.Y[x].extend(y_merge)
            }
            c_cur_x[y_merge] = cxy
        }
    }
    return mca
}

// Calculate Mca solution

func SolveMca(mca McaType) SolType {
    var lx, ly = make([]int, mca.NR), make([]int, mca.NC+mca.NR)
    var xy, yx = make([]int, mca.NR), make([]int, mca.NC+mca.NR)
    var S, T = make([]bool, mca.NR), make([]bool, mca.NC+mca.NR)
    var slack = make([]int, mca.NC+mca.NR)
    var slackx = make([]int, mca.NC+mca.NR)
    var rev_path = make([]int, mca.NR)

    // debug
    var mark_pos = make([]int, mca.NR)  // 0 initialized
    // stat
    var stat_upper, stat_mid, stat_lower, stat_marks, stat_improve int
    var match_cnt = 0

    var row_order = make([]int, mca.NR)
    var q = make([]int, mca.NR) // this is a queue, we maintain it manually with q_front, q_back
    var q_front, q_back int
    var slack0_pivot int
    var found_augmentation bool
    var x_aug, y_aug int
    var x int

    var y_merge, prev_y_merge = make([]int, mca.NC+mca.NR), make([]int, mca.NC+mca.NR)
    var y_merge_back, prev_y_merge_back = 0, 0

    merge_relevant_columns := func (x int) {
        if q_back > 1 {
            y_merge, prev_y_merge = prev_y_merge, y_merge
            prev_y_merge_back = y_merge_back
            y_merge_back = 0

            y_next := INF
            if len(mca.Y[x]) > 0 { y_next = Min(y_next, mca.Y[x][0]) }
            if prev_y_merge_back > 0 { y_next = Min(y_next, prev_y_merge[0]) }

            var i_y, prev_i_y = 0, 0
            for y_cur := 0 ; y_next < INF ; {
                y_cur, y_next = y_next, INF
                y_merge[y_merge_back] = y_cur
                y_merge_back ++

                for i_y < len(mca.Y[x]) && mca.Y[x][i_y] <= y_cur { i_y ++ }
                if  i_y < len(mca.Y[x]) {
                    y_next = Min(y_next, mca.Y[x][i_y])
                }

                for prev_i_y < prev_y_merge_back && prev_y_merge[prev_i_y] <= y_cur { prev_i_y ++ }
                if  prev_i_y < prev_y_merge_back {
                    y_next = Min(y_next, prev_y_merge[prev_i_y])
                }
            }
        } else {
            y_merge_back = 0
            for i_y := 0; i_y < len(mca.Y[x]); i_y ++ {
                y_merge[y_merge_back] = mca.Y[x][i_y]
                y_merge_back ++
            }
        }
    }

    add_to_tree := func (x int, greedy_break bool) {
        if xy[x] != -1 {   // not for root of path
            T[xy[x]] = true
        }

        q[q_back] = x
        q_back ++
        S[x] = true

        // update slack and slackx
        for i_y := 0; i_y < len(mca.Y[x]); i_y ++ {
            y := mca.Y[x][i_y]
            cxy := mca.C[x][i_y]

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
        // check_slack(slack, slackx, mca, lx, ly)

        merge_relevant_columns(x)
    }

    augment := func (x int, y int) {
        x_aug, y_aug = x, y
        found_augmentation = true
    }

    // try to grow a path in the equality subgraph
    explore_eq_subgraph := func () {
        if debug { fmt.Println("explore equality subgraph") }

        for q_back > q_front && !found_augmentation {
            x = q[q_front]
            q_front ++

            // plausible heuristic: first try x_aug+1 with y_aug+1
            if match_cnt > 0 {
                y := y_aug + 1
                if !T[y] && yx[y] == -1 && mca.cost(x, y) == lx[x] + ly[y] {
                    augment(x, y)
                    stat_upper ++ // stat
                    continue
                }
            }

            // sparse heuristic: try only y with positive cost
            for i_y := 0; i_y < len(mca.Y[x]); i_y ++ {
                y := mca.Y[x][i_y]
                if mca.C[x][i_y] == lx[x] + ly[y] && !T[y] {
                    if yx[y] == -1 { // if y is unmatched
                        augment(x, y)
                        stat_mid ++ // stat
                        break
                    }

                    rev_path[yx[y]] = x
                    add_to_tree(yx[y], false)
                }
            }
        }
    }

    // improve labels lx, ly to grow the equality subgraph
    improve_labelling := func() {
        if debug { fmt.Println("improve labelling") }
        stat_improve ++ // stat

        delta := INF

        // for (var y=0; y<mca.NC+mca.NR; y++)
        //     if (slack[y] < INF && y_merge.indexOf(y) == -1) // && mca.cost(slackx[y], y) > 0)
        //         fmt.Println(y_merge + ' y ' + y + ' slack[y] ' + slack[y] + ' cost ' + mca.cost(slackx[y], y))

        // for (var y = 0; y < mca.NC+mca.NR; y++) {
        for i_y := 0; i_y < y_merge_back; i_y ++ {
            y := y_merge[i_y]
            if !T[y] && slack[y] < delta {
                delta = slack[y]
            }
        }

        if debug && delta > 1000000000 { fmt.Println("no minimal delta found") }

        // adapt node labels: q holds all nodes in S
        lx[q[0]] -= delta
        for i := 1; i < q_back; i ++ {
            x := q[i]
            lx[x] -= delta
            ly[xy[x]] += delta
        }

        slack0_pivot = -1
        // for (var y = 0; y < mca.NC+mca.NR; y++) {
        for i_y := 0; i_y < y_merge_back; i_y ++ {
            y := y_merge[i_y]
            if !T[y] && slack[y] < INF {
                slack[y] -= delta
                if slack[y] == 0 {
                    // store the first with ![T[y] && slack[y]==0, it exists
                    if slack0_pivot == -1 {
                        slack0_pivot = y
                    }
                    // condition leads immediately to a match in grow_eq_subgraph
                    if (yx[y] == -1) {
                        slack0_pivot = y
                        break
                    }
                }
            }
        }
    }

    // grow path in the equality subgraph
    grow_eq_subgraph := func() {
        if debug { fmt.Println("grow equality subgraph") }
        // search all y: this ensures path augmentation
        // for (var y = slack0_pivot; y < mca.NC; ++y)
        for slack0_pivot != -1 {
            y := slack0_pivot
            slack0_pivot = -1

            // if (! (!T[y] && slack[y] == 0)) break
            if yx[y] == -1 {    // if y is unmatched
                augment(slackx[y], y)
                stat_lower ++ // stat
                break
            } else {
                stat_marks ++ // stat
                mark_pos[match_cnt] ++
                // if (!S[yx[y]])     // never false because: S[xy[y]] => T[y]
                rev_path[yx[y]] = slackx[y]
                add_to_tree(yx[y], true)
                // else T[y] = true; // is true anyways, as yx[y] cannot be the path root
            }
        }
    }

    // start with nothing matched
    Init(xy, -1)
    Init(yx, -1)

    // prepare sparse problem
    for x := 0; x < mca.NR; x ++ {
        mca.Y[x].extend(mca.NC + x)
        mca.C[x].extend(0)
    }

    // setup lx, ly feasible
    Init(ly, 0)
    for x := range mca.C {
        lx[x], _ = mca.C[x].Max()
    }

    for preprocessed := false; ; preprocessed = true {
        // prepare sorting criteria and order slice
        rowmax, rowmax_i := make([]int, mca.NR), make([]int, mca.NR)
        for x := range mca.C {
            if xy[x] != -1 {
                rowmax[x], rowmax_i[x] = 0, -1  // is already greedy matched
            } else {
                rowmax[x], rowmax_i[x] = mca.C[x].Max()
            }
            row_order[x] = x
        }

        /* Algorithm performance strongly depends on the successive choice of x. We seperate
         * sparse rows from the rest and sort them by decreasing cost, to prefer longer matches.
         * The rest is sorted by decreasing sparsity, to match rows with many options at the end.
         * Zeros are put at the very end. Cost refers to row maxima - the most probable match.
        */

        non_sparse := func(sl IntSl) bool { return len(sl) > 5 }   // threshold to considered row dense

        var tmp_order PriorityIntSl
        tmp_order.IntSl = row_order
        tmp_order.less = func(i, j int) bool {
            switch {
                case rowmax[i] == 0 && rowmax[j] == 0: return i < j // natural order (empty rows)
                case rowmax[i] == 0: return false                   // non-matching last
                case rowmax[j] == 0: return true                    // non-matching last
                case non_sparse(mca.Y[i]) && non_sparse(mca.Y[j]):
                    return len(mca.Y[i]) < len(mca.Y[j])            // natural order (dense rows)
                case non_sparse(mca.Y[i]): return false             // sparse rows first
                case non_sparse(mca.Y[j]): return true              // sparse rows first
                default: return rowmax[j] < rowmax[i]               // natural order (no rows)
            }
        }

        sort.Sort(tmp_order)

        // sort once again after greedy match
        if preprocessed || skip_preproc {
            break
        }

        // greedy match what is possible
        for i_x := range row_order {
            x := row_order[i_x]
            if rowmax[x] == 0 {
                continue
            }
            y := mca.Y[x][rowmax_i[x]]
            if xy[x] == -1 && yx[y] == -1 {
                xy[x], yx[y] = y, x
                match_cnt ++
            }
        }
    }

    fmt.Printf("preprocessor matched %4.2f %%\n", 100.0 * Fp(match_cnt) / Fp(mca.NR))

    InitBool(S,false)
    InitBool(T,false)
    Init(rev_path, -1)
    Init(slack, INF)

    // main loop to grow the matching
    for ; match_cnt < mca.NR; match_cnt ++ {
        q_front, q_back = 0, 0

        // select unmatched x as path root
        for i_x := range mca.C {
            x := i_x
            if sort_by_cost { x = row_order[i_x] }

            if xy[x] == -1 {   // if x is unmatched
                add_to_tree(x, false)
                break
            }
        }

        if debug && q_back == q_front {
            fmt.Println("unexpected empty queue")
            break
        }

        found_augmentation = false
        // augment path to unmatched y
        for !found_augmentation {
            slack0_pivot = -1

            explore_eq_subgraph()

            if found_augmentation { break }

            improve_labelling()

            q_front = q_back   // don't loose old entries of q; need them for sparse resetting

            grow_eq_subgraph()
        }

        // flip edges along augmenting path, thereby growing the matching by one
        if debug { fmt.Println("flipping edges") }
        for cx, cy, ty := x_aug, y_aug, 0; cx != -1 ; cx, cy = rev_path[cx], ty {
            ty = xy[cx]
            yx[cy] = cx
            xy[cx] = cy
        }

        if debug { fmt.Println("flipped edges") }

        // the following early stopping criterion is correct if:
        // there has to be some free x, with some y*, with better cost(x, y*) than cost(yx[y], y)

        if early_stop && match_cnt % 10 == 0 {
            only_crap_left := true
            for x := range mca.C {
                if xy[x] != -1 { continue }
                for i_y := range mca.Y[x] {
                    y := mca.Y[x][i_y]
                    if yx[y] == -1 || mca.C[x][i_y] > mca.cost(yx[y], y) {
                        only_crap_left = false
                        break
                    }
                }
                if ! only_crap_left { break }
            }
            if only_crap_left { break }
        }

        // correct, because matching is already flipped
        for i := 0; i < q_back; i ++ {
            x := q[i]
            S[x] = false
            T[xy[x]] = false
            rev_path[x] = -1
        }

        for i_y := 0; i_y < y_merge_back; i_y ++ {
            slack[y_merge[i_y]] = INF
        }

        // correctness check, on suboptimal matchings uncomment to check if sparse resets work
        if debug {
            if (match_cnt % 10 == 0) {
                for x := 0; x < mca.NR; x ++ {
                    if S[x] != false || rev_path[x] != -1 || xy[x] != -1 && x != yx[xy[x]] { fmt.Println("f*** x") }
                }
                for y := 0; y < mca.NC; y ++ {
                    if T[y] != false || yx[y] != -1 && y != xy[yx[y]] { fmt.Println("f*** y") }
                }
            }
            // S.init(false);  T.init(false);  rev_path.init(-1)
        }
    } // for match_cnt

    // if debug {
    //     // log matches
    //     for(var x=1; x<mca.NR; x++) {
    //         if (mca.cost(x, xy[x]) != mca.cost(x-1, xy[x]-1)) {
    //             var c = mca.cost(x,xy[x])
    //             fmt.Println(''.concat('cost(', x, ',', xy[x], ') = ', c, '\n-> ', mca.a.slice(x, x+c), '\n-> ', mca.b.slice(xy[x], xy[x]+c)))
    //             fmt.Println('+++++++++++++++++++++++++++++++++++++++++++')
    //         }
    //     }

    //     // log highest marked positions
    //     var mark_pos_id = Array(mca.NR)
    //     for (var i=0; i<mark_pos_id.length; i++) mark_pos_id[i] = i
    //     mark_pos_id.sort(function(x,y) {return mark_pos[y]-mark_pos[x];})
    //     for (var i=0; i<10; i++) fmt.Println(mark_pos_id[i] + ': ' + mark_pos[mark_pos_id[i]])
    // }

    var sol = SolType{Cost:0, XY:xy, YX:yx}
    total_match := 0
    for x := 0; x < mca.NR; x ++ {
        if xy[x] == -1 {
            continue
        }
        if c := mca.cost(x, xy[x]); c > 0 {
            sol.Cost += c
            total_match ++
        }
    }

    fmt.Printf("upper: %d, mid: %d, improve %d, lower: %d, marks: %d\n",
        stat_upper, stat_mid, stat_improve, stat_lower, stat_marks)
    fmt.Printf("early stopping at %4.2f %%, saved %d cycles, stopped with %d matches\n",
        100.0 * Fp(match_cnt) / Fp(mca.NR), mca.NR - match_cnt, match_cnt)
    fmt.Printf("total cost: %d, total match: %d, nr: %d, nc: %d,  minMatch: %d, swap: %t\n",
        sol.Cost, total_match, mca.NR, mca.NC, mca.minMatch, mca.Swap)

    return sol

}

