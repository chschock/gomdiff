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

func (self IntSl) Len() int {
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
    for i := 0; i < len(mca.seq_b) -1; i ++ {
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

        prev_y_indices = y_indices

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
                    cxy ++
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

var debug = false
var sort_by_cost = false
var minMatch = 3

type solType struct {
    xy []int
    yx []int
    cost int
}

func SolveMca(mca mcaType) solType {
    var lx, ly = make([]int, mca.nr), make([]int, mca.nc+mca.nr)
    var xy, yx = make([]int, mca.nr), make([]int, mca.nc+mca.nr)
    var S, T = make([]bool, mca.nr), make([]bool, mca.nc+mca.nr)
    var slack = make([]int, mca.nc+mca.nr)
    var slackx = make([]int, mca.nc+mca.nr)
    var rev_path = make([]int, mca.nr)

    // debug
    var mark_pos = make([]int, mca.nr)  // 0 initialized
    // stat
    var stat_upper, stat_mid, stat_lower, stat_marks, stat_improve int
    var match_cnt = 0

    var row_order = make([]int, mca.nr)
    var q = make([]int, mca.nr) // this is a queue, we maintain it manually with q_front, q_back
    var q_front, q_back int
    var slack0_pivot int
    var found_augmentation bool
    var x_aug, y_aug int
    var x int

    var y_merge, prev_y_merge = make([]int, mca.nc+mca.nr), make([]int, mca.nc+mca.nr)
    var y_merge_back, prev_y_merge_back = 0, 0

    merge_relevant_columns := func (x int) {
        if q_back > 1 {
            y_merge, prev_y_merge = prev_y_merge, y_merge
            prev_y_merge_back = y_merge_back
            y_merge_back = 0

            y_next := INF
            if mca.y[x].Len() > 0 { y_next = Min(y_next, mca.y[x].sl[0]) }
            if prev_y_merge_back > 0 { y_next = Min(y_next, prev_y_merge[0]) }

            var i_y, prev_i_y = 0, 0
            for y_cur := 0 ; y_next < INF ; {
                y_cur, y_next = y_next, INF
                y_merge[y_merge_back] = y_cur
                y_merge_back ++

                for i_y < mca.y[x].Len() && mca.y[x].sl[i_y] <= y_cur { i_y ++ }
                if  i_y < mca.y[x].Len() {
                    y_next = Min(y_next, mca.y[x].sl[i_y])
                }

                for prev_i_y < prev_y_merge_back && prev_y_merge[prev_i_y] <= y_cur { prev_i_y ++ }
                if  prev_i_y < prev_y_merge_back {
                    y_next = Min(y_next, prev_y_merge[prev_i_y])
                }
            }
        } else {
            y_merge_back = 0
            for i_y := 0; i_y < mca.y[x].Len(); i_y ++ {
                y_merge[y_merge_back] = mca.y[x].sl[i_y]
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
        for i_y := 0; i_y < mca.y[x].Len(); i_y ++ {
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
            for i_y := 0; i_y < mca.y[x].Len(); i_y ++ {
                y := mca.y[x].sl[i_y]
                if mca.c[x].sl[i_y] == lx[x] + ly[y] && !T[y] {
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

        // for (var y=0; y<mca.nc+mca.nr; y++)
        //     if (slack[y] < INF && y_merge.indexOf(y) == -1) // && mca.cost(slackx[y], y) > 0)
        //         fmt.Println(y_merge + ' y ' + y + ' slack[y] ' + slack[y] + ' cost ' + mca.cost(slackx[y], y))

        // for (var y = 0; y < mca.nc+mca.nr; y++) {
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
        // for (var y = 0; y < mca.nc+mca.nr; y++) {
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
        // for (var y = slack0_pivot; y < mca.nc; ++y)
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
    for x := 0; x < mca.nr; x ++ {
        mca.y[x].extend(mca.nc + x)
        mca.c[x].extend(0)
    }

    // setup lx, ly feasible
    lx = make([]int, mca.nr)
    Init(ly, 0)
    c_pp := make([]IntSl, mca.nr)

    for x := range mca.c {
        c_pp[x].sl = make([]int, mca.c[x].Len())   // copy of cost for preprocessing
        copy(c_pp[x].sl, mca.c[x].sl)
        lx[x], _ = mca.c[x].Max()
    }


    // for testing
    for x := range mca.c {
        row_order[x] = x
    }


    fmt.Println("preprocessor matched %4.2f %%", 100.0 * match_cnt / mca.nr)

    InitBool(S,false)
    InitBool(T,false)
    Init(rev_path, -1)
    Init(slack, INF)


    // main loop to grow the matching
    for ; match_cnt < mca.nr; match_cnt ++ {
        q_front, q_back = 0, 0

        // select unmatched x as path root
        for x_o := range mca.c {
            x := x_o
            if sort_by_cost { x = row_order[x_o] }

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

        // if (early_stop && match_cnt % 10 == 0) {
        //     var only_crap_left = true
        //     for (var x=0; x<mca.nr; x++) {
        //         if (xy[x] != -1) continue
        //         for (var i_y=0; i_y<mca.y[x].length - 1; i_y++) { // ignore the 0-cost node
        //             var y = mca.y[x][i_y]
        //             if (yx[y] == -1 || mca.c[x][i_y] > mca.cost(yx[y], y)) {
        //                 only_crap_left = false
        //                 break
        //             }
        //         }
        //         if (!only_crap_left) break
        //     }
        //     if (only_crap_left) break
        // }

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
                for x := 0; x < mca.nr; x ++ {
                    if S[x] != false || rev_path[x] != -1 || xy[x] != -1 && x != yx[xy[x]] { fmt.Println("f*** x") }
                }
                for y := 0; y < mca.nc; y ++ {
                    if T[y] != false || yx[y] != -1 && y != xy[yx[y]] { fmt.Println("f*** y") }
                }
            }
            // S.init(false);  T.init(false);  rev_path.init(-1)
        }
    } // for match_cnt

    // if debug {
    //     // log matches
    //     for(var x=1; x<mca.nr; x++) {
    //         if (mca.cost(x, xy[x]) != mca.cost(x-1, xy[x]-1)) {
    //             var c = mca.cost(x,xy[x]);
    //             fmt.Println(''.concat('cost(', x, ',', xy[x], ') = ', c, '\n-> ', mca.a.slice(x, x+c), '\n-> ', mca.b.slice(xy[x], xy[x]+c)));
    //             fmt.Println('+++++++++++++++++++++++++++++++++++++++++++');
    //         }
    //     }

    //     // log highest marked positions
    //     var mark_pos_id = Array(mca.nr);
    //     for (var i=0; i<mark_pos_id.length; i++) mark_pos_id[i] = i;
    //     mark_pos_id.sort(function(x,y) {return mark_pos[y]-mark_pos[x];});
    //     for (var i=0; i<10; i++) fmt.Println(mark_pos_id[i] + ': ' + mark_pos[mark_pos_id[i]]);
    // }

    var sol = solType{xy:xy, yx:yx, cost:0}
    total_match := 0
    for x := 0; x < mca.nr; x ++ {
        if xy[x] == -1 {
            continue
        }
        if c := mca.cost(x, xy[x]); c > 0 {
            sol.cost += c
            total_match ++
        }
    }

    fmt.Println("upper: %d, mid: %d, improve %d, lower: %d, marks: %d",
        stat_upper, stat_mid, stat_improve, stat_lower, stat_marks)
    fmt.Println("early stopping at %4.2f %%, saved %d cycles, stopped with %d matches",
        100.0 * match_cnt / mca.nr, mca.nr - match_cnt, match_cnt)
    fmt.Println("total cost: %d, total match: %d, nr: %d, nc: %d,  minMatch: %d, swap %d",
        sol.cost, total_match, mca.nr, mca.nc, minMatch, mca.swap)

    return sol;

}

func main() {
    a := "ksajldkajsldkjalskd"
    b := "ksajldkajkkkkksldkjooooalsppkd"

    mca := MakeMca(a, b, 2)
    fmt.Printf("%+v\n", mca)
    sol := SolveMca(mca)
    fmt.Printf("%+v\n", sol)
}
