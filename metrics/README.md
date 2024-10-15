# Compiling TMalign

Our TMalign software is Version 20220412, compiled from the .cpp source as suggested by the authors ([link](https://seq2fun.dcmb.med.umich.edu//TM-align/)).


> Click TMalign.cpp (last update: 2022/4/12) and readme.c++.txt to download the newest version of the TM-align source code in C++. You can compile the program in your Linux computer by (you can ignore the '-static' option for some machines, such as Mac, which does not support static build):
>
> ```
> g++ -static -O3 -ffast-math -lm -o TMalign TMalign.cpp
> ```

The -static flag did not compile in default Sherlock environment. We removed it and ran
```
g++ -O3 -ffast-math -lm -o TMalign TMalign.cpp
```

We also downloaded the `TMalign_cpp` executable directly from ([link](https://seq2fun.dcmb.med.umich.edu//TM-align/)). Though it is "strongly recommended to download the TM-align source code and compile it on your machine" for speed, a quick test running the first 500 pairs of `training/data/tmaln-06.out` showed negligible speed difference between `TMalign_cpp` and our (non-statically) compiled `TMalign` (18.68 vs. 19.00 respectively). 



# Running TMalign and computing extra metrics

## Run TMalign on a protein pairs. 

By default the pair files is `training/data/tmaln-06_500.out`.

Output of TMalign is saved to local file `tmalign.out`.

Transformations (rotation matrix and translation vector to superpose first structure onto second structure) are saved to local dir `rot_mats` by default.

```
./run_tmalign.sh
```

## Parse the TMalign output file into tsv

Parses fields from `tmalign.out` into tab-separated `tmalign.csv`, including computing the CIGAR string.

```
python parse_tmalign_output.py
```


## Enrich the TMalign tabulted results with other metrics

Reads 'tmalign.csv' and enriches it with
- LDDT (Approximate implementation from Foldseek. Computed between aligned residues only)
- Chamfer distance (Using superposition found by TMalign, saved to `rot_mats`.)
- Earth movers distance (Using superposition found by TMalign, saved to `rot_mats`.)

Saves updated table with these metrics to `tmalign_extra.csv`.

```
python compute_extra_metrics.py
```
