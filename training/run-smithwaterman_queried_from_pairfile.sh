#!/bin/bash
GAPOPEN=$1
GAPEXTEND=$2
PAIRFILE=$3
OUTFILE=$4

if [ "$#" -ne 4 ]; then
    echo "Illegal number of parameters"
    exit
fi

wc $PAIRFILE -l

mkdir -p tmp/alignments_queried_from_pairfile

open_sem(){
    mkfifo pipe-$$
    exec 3<>pipe-$$
    rm pipe-$$
    local i=$1
    for((;i>0;i--)); do
        printf %s 000 >&3
    done
}

# run the given command asynchronously and pop/push tokens
run_with_lock(){
    local x
    # this read waits until there is something to read
    read -u 3 -n 3 x && ((0==x)) || exit $x
    (
     ( "$@"; )
    # push the return code of the command to the semaphore
    printf '%.3d' $? >&3
    )&
}

task(){
    ssw/ssw_test -o "$GAPOPEN" -e "$GAPEXTEND" -a tmp/s.mat -p \
        tmp/target.fasta "$1" \
        2> /dev/null \
    | awk '/^target/{target=$2}
            /^query/{query=$2}
            /^optimal_alignment_score/{score=$2; print query,target,score}' \
    | awk 'NR==FNR {a[$1 " " $2]=1; next} (($1 " " $2) in a) {print $1 " " $2 " " $3}' "$PAIRFILE" - \
    >  "tmp/alignments_queried_from_pairfile/${1##*/}"".m8";
}


cp tmp/sub_score.mat tmp/s.mat  # must be shorter than 16 characters

N=64
open_sem $N
for thing in tmp/splits/*; do
    run_with_lock task "$thing"
done

(while [[ $(pidof ssw_test) ]]; do sleep 1; done)

cat tmp/alignments_queried_from_pairfile/*.m8 \
    | sort -k1,1 -k3,3nr \
    > "tmp/$OUTFILE"
rm -r tmp/alignments_queried_from_pairfile

wc "tmp/$OUTFILE" -l