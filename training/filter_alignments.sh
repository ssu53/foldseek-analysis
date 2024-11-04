# Filter the lines of alignments output files to only those pairs contained in INDEX_FILE

# INDEX_FILE='foo.out'
# INDEX_FILE='../metrics/pairfile_random_subset.out'
# INDEX_FILE='data/tmaln-06.out'
INDEX_FILE='../metrics//pairfile_val.out'

for data_file in tmp21_v6_encodings/alignments/*; do
    echo $data_file
    awk 'NR==FNR {a[$1 " " $2]=1; next} (($1 " " $2) in a) {print $1 " " $2 " " $3}' \
        $INDEX_FILE $data_file
done

# awk 'NR==FNR {a[$1 " " $2]=1; next} (($1 " " $2) in a) {print $1 " " $2 " " $3}' ../metrics/pairfile_random_subset.out tmp21_v6_encodings/alignments/split_00.fasta.m8