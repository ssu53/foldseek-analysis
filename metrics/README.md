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


TMalign metrics and other metrics are computed for a list of protein pairs. 
The pairs are sampled uniformly from the files contained in `ALIGNMENT_DIR`, which is the `alignments/` output directory produced by running `learnAlphabet.sh`.


```
./run_tmalign.sh ALIGNMENT_DIR
```

This script creates a datestamped directory within `metrics/` to hold outputs.

The slowest step by far is computing the EMD. 

High level overview of the step sbelow.


## Sample list of pairs from alignments

There are likely too many alignments to compute metrics for all of them. Sample uniformly over those files:

```
python get_sampled_pairfile.py
```

This can be skipped if already have a pair file, e.g. `training/data/tmaln-06.out`.

## Run TMalign on protein pairs

Run TMalign for pairs in a pairfile.

```
TMalign_cpp protein1.pdb protein2.pdb -m rot_mats/protein1-protein2.txt
```

Raw outputs of TMalign is dumped to `tmalign.out`.

Transformations (rotation matrix and translation vector to superpose first structure onto second structure) are saved to `rot_mats/`.



## Parse the TMalign output file into tsv

Parses fields from `tmalign.out` into tab-separated `tmalign.csv`, including computing the CIGAR string.

```
python parse_tmalign_output.py
```


## Enrich the TMalign tabulated results with other metrics

Reads 'tmalign.csv' and enriches it with
- LDDT (Approximate implementation from Foldseek. Computed between aligned residues only)
- Chamfer distance (Using superposition found by TMalign, saved to `rot_mats`.)
- Earth movers distance (Using superposition found by TMalign, saved to `rot_mats`.)

Saves updated table with these metrics to `tmalign_extra.csv`.

```
python compute_extra_metrics.py
```
