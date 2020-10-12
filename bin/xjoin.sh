#!/bin/bash 

# xjoin take as inputs 3 files containing index number when they mapped 
# it returns a file containing: read ids and their mapped indexes
# ex: AACCRT:15000 C56 D01 B65

xjoin() {
    local f
    local srt="cat"

    # If there is less than 3 args, join them in "join"
    if [ "$#" -lt 3 ]; then
        join <($srt "$1") <($srt "$2")
    # else, join the three in join
    else
        f=$1
        shift
        join <($srt "$f") <(xjoin "$@")
    fi
}

if [ "$#" -eq "0" ]; then
    exit 1
elif [ "$#" -eq 1 ]; then
    cat $1
else
    xjoin "$@" | cat
fi
