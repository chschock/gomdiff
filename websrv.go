package main

import (
    "fmt"
    "log"
    "strconv"
    "encoding/json"
    "time"
    "net/http"
    "github.com/gorilla/mux"
    "github.com/davecheney/profile"
)

const gen_profile_data = false

type Output struct {
    Sol SolType     `json:"sol"`
    Mca McaType     `json:"mca"`
}

func about(w http.ResponseWriter, r * http.Request){
    fmt.Fprintf(w, "MatchDiff Webservice API\ncall solve(a, b, minMatch)")
}

func solve(w http.ResponseWriter, r * http.Request){
    w.Header().Set("Server","MatchDiff Webservice")
    w.Header().Set("Content-Type", "text/json")

    a := r.URL.Query().Get("a")
    b := r.URL.Query().Get("b")
    minMatch, _ := strconv.Atoi(r.URL.Query().Get("minMatch"))

    // fmt.Printf("a: %+v\nb: %+v\n minMatch: %+v\n", a, b, minMatch)

    start := time.Now()

    mca := MakeMca(a, b, minMatch)
    fmt.Printf("MakeMca took %s sec\n", time.Since(start))

    sol := SolveMca(mca)
    fmt.Printf("MakeMca + SolveMca took %s sec\n", time.Since(start))

    json.NewEncoder(w).Encode(Output{Mca: mca, Sol: sol})
    fmt.Printf("MakeMca + SolveMca + Write took %s sec\n", time.Since(start))
}


func main() {
    if gen_profile_data {
        defer profile.Start(profile.CPUProfile).Stop()
    }

    router := mux.NewRouter().StrictSlash(true)
    router.HandleFunc("/", about)
    router.HandleFunc("/solve", solve)

    log.Fatal(http.ListenAndServe(":8080", router))
}