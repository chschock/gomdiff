package main

import (
    "bytes"
    "strconv"
)

type IntSlSl []IntSl

func (self IntSlSl) MarshalJSON() ([]byte, error) {
    var buf bytes.Buffer
    buf.WriteRune('[')
    for r, row := range self {
        if r > 0 { buf.WriteRune(',') }
        buf.WriteRune('[')
        for c, col := range row {
            if c > 0 { buf.WriteRune(',') }
            buf.WriteString(strconv.Itoa(col))
        }
        buf.WriteRune(']')
    }
    buf.WriteRune(']')
    return buf.Bytes(), nil
}
