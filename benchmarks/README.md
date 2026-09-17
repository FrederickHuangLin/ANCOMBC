# benchmarks

Verification and timing material for ANCOMBC 2.15.1. `optimization_report.tex` and `optimization_report.pdf` are the report; `timings_2.13.2_vs_2.15.1.csv` is the back-to-back comparison of the 17 harness workloads.

`harness/`: `workloads.R` (workload definitions), `run.R` (timing and peak memory), `compare.R` (result equality), `edge_cases.R` (nine edge cases), `repro_pseudo0.R` (the `pseudo = 0` reproduction). All paths are taken from command-line arguments; with `lib_a` and `lib_b` the library directories of the two installations to compare:

    Rscript harness/run.R lib_a out_a all 3
    Rscript harness/run.R lib_b out_b all 3
    Rscript harness/compare.R out_a out_b 1e-8
