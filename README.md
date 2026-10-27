# RIQKDNumerics

Code and data for the key rates of

> A. Mishra, M. Masini, M. Balanzó-Juandó, A. Tavakoli and S. Pironio,
> *Prepare-and-measure quantum key distribution from a bound on information leakage*.

The code computes lower bounds on the asymptotic key rate of prepare-and-measure QKD protocols under the
leakage bound (LB, called restricted information, RI, in earlier versions and in the name of this package). The
assumption is that the probability of guessing Alice's input from everything that leaves her laboratory is at
most G. It also computes the entropy certified by a CHSH-type value
alone.

The bound on the conditional entropy H(A|E) is computed in three steps:
- the Koßmann–Schwonnek relative-entropy bound on a fixed partition of [0, 1] (paper, Sec. V.A and App. C);
- one semidefinite program per node of the partition, from the one-factor moment relaxation of Sec. V.B;
- the dual of each program, solved with MOSEK.

The programs are generated with [PCPOP.jl](https://github.com/abh1mishra/PCPOP.jl). The rate per key round is
K = H(A|E) − H(A|B). Thresholds are the zero crossings of K, located by bisection.

The package covers the key-rate results of Sec. VI and Apps. E–F, and the CHSH relation of Sec. VI.D. It does not
include the linear programs of Sec. III and App. A (Table II), or the device-independent comparison of
Sec. VI.C. Section and table numbers refer to the arXiv version of the paper.

**Units of G.** The code works with the guessing probability P_g = G/n, where n is the number of preparations,
and its `G`, `dG` and `G_honest` are the paper's G, δG and G_hon divided by n. The honest value is `G` = 2/n for
every protocol, i.e. G = 2 in the paper, and a file at `dG0.01` of a four-state protocol is the paper's δG = 0.04,
while `dG0.006667` of (6,3) and `dG0.013333` of (3,3) are both δG = 0.04.

## Installation

Requirements:
- Julia 1.11 (the results were computed with 1.11.7, and this package was tested with 1.11.3);
- a MOSEK licence in `~/mosek/mosek.lic` (free for academic use, <https://www.mosek.com/products/academic-licenses/>).

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

This installs PCPOP.jl at commit `813b644` and the other packages at the versions of `Manifest.toml`. When
Mosek.jl is built, it downloads the MOSEK 11.2 binaries. To use an existing MOSEK 11.2 installation instead, set
`MOSEKBINDIR` to its `bin` directory before instantiating. The results were computed with MOSEK 11.2.4.

## Usage

All commands are run from this directory. By default the scripts write to `results/`, which holds the results of
the paper. The examples below therefore write to a new directory, `rerun/`. Without it, `run_job.jl` would
overwrite a shipped result, and the batch scripts would skip every job, since all shipped jobs are done.

```bash
# one key-rate job: a curve, and the thresholds of the von Neumann and min-entropy rates
julia --project=. scripts/run_job.jl prot=43 rel=lean:3 dG=0.00 sweep=v out=rerun   # → rerun/43__lean_3__dG0.00__v__g26.json

# all the key-rate jobs of the paper, or all its CHSH scans, with N parallel workers, written under rerun/
scripts/run_batch.sh jobs/keyrates.txt N rerun
scripts/run_chsh_batch.sh jobs/chsh.txt N rerun

# the largest relaxation dG of the bound at which the ideal behaviour (v = eta = 1) still certifies key
julia --project=. scripts/dG_scan.jl prot=43 rel=lean:3 out=rerun/dG  # → rerun/dG/dG__43__lean_3__v1.0__eta1.0__g26.json

# one CHSH scan
julia --project=. scripts/chsh_scan.jl G=0.5 out=rerun/chsh      # → rerun/chsh/chsh__lean_4__G0.5__g26.json

# a single point, e.g. to compare relaxations or partitions
julia --project=. scripts/point.jl prot=43 rel=lean:3 v=0.852 grid=fine120

# checks of the relaxation (App. F) and residual correction of the dual solutions (App. E), a few minutes each
julia --project=. scripts/validate.jl
julia --project=. scripts/residual_check.jl
```

The first lines of each script document all its options.
- **Output root.** The third argument of the batch scripts replaces the leading `results` of each job's output
  directory. For example, `results/curves` becomes `rerun/curves`.
- **Relaunching.** The batch scripts skip every job whose result file is marked `"status": "done"`, so an
  interrupted batch can be relaunched. An unfinished job restarts from the beginning.
- **Run times.** Most jobs take minutes on one core. The slowest take up to about ten hours; see the header of
  `jobs/keyrates.txt`.

## Layout

| path | content |
|---|---|
| `src/RIQKDNumerics.jl` | module |
| `src/protocols.jl` | protocols, honest statistics, H(A|B), honest value of G |
| `src/relaxation.jl` | construction of the semidefinite programs (state-weighted and, for validation only, tracial) |
| `src/entropy.jl` | partitions, entropy bound, min-entropy, key rate at one point |
| `src/utils.jl` | Shannon entropies and the coefficients of the entropy bound |
| `scripts/run_job.jl`, `scripts/run_batch.sh` | key-rate curves and thresholds |
| `scripts/dG_scan.jl` | largest relaxation of the bound with a positive rate |
| `scripts/chsh_scan.jl`, `scripts/run_chsh_batch.sh` | entropy as a function of the CHSH value |
| `scripts/point.jl` | a single point |
| `scripts/validate.jl`, `scripts/residual_check.jl` | checks of App. F and App. E |
| `jobs/keyrates.txt`, `jobs/dG_scan.txt`, `jobs/chsh.txt` | all the jobs whose results appear in the paper, and those of the (4,2) protocol of its earlier versions |
| `results/`, `results/curves/`, `results/dG/`, `results/chsh/` | their results |

## Protocols

The names used by the code, with their states, Bob's measurements and key-generating rounds:

| name | paper | states | measurements | key |
|---|---|---|---|---|
| `42` | (4,2), earlier versions only | \|0⟩, \|1⟩, \|+⟩, \|−⟩ | Z, (Z+X)/√2 | Z |
| `43` | (4,3) | as above | Z, (Z±X)/√2 | Z |
| `44` | (4,4) | as above | Z, (Z±X)/√2, X | Z and X |
| `6sym` | (6,3) | eigenstates of three observables at Bloch angles 0°, 60°, 120° | the same three observables | all three bases |
| `3st_ZBB` | (3,3) | \|0⟩, \|1⟩, \|+⟩ | Z, (Z±X)/√2 | Z |
| `63` | (6,3) variant, Table V | eigenstates of Z, X, (Z+X)/√2 | Z, X, (Z+X)/√2 | all three bases |
| `BB84+1` | variant, Table V | BB84 states | Z, X, (Z+X)/√2 | Z and X |
| `3st_ZB`, `3st_ZXB` | variants, Table V | \|0⟩, \|1⟩, \|+⟩ | {Z, (Z+X)/√2}, {Z, X, (Z+X)/√2} | Z |
| `3st_ZX` | three-state protocol, Sec. III.A and Table V | \|0⟩, \|1⟩, \|+⟩ | Z, X | Z (LB-classical: no key) |
| `42@φ` | tilted (4,2), earlier versions only | as `42` | Z, cos φ Z + sin φ X | Z |
| `43@φ` | tilted (4,3), App. F | as `43` | Z, cos φ Z ± sin φ X | Z |
| `4chsh` | CHSH relation, Sec. VI.D | BB84 states | (Z±X)/√2 | Z-basis bit |

The paper labels the states x = (θ, a) by key basis θ and raw key symbol a. For the four-state and three-state
protocols the code numbers them 1–4 in the order |0⟩, |1⟩, |+⟩, |−⟩, that is, (1,0), (1,1), (2,0), (2,1).

Noise model (Sec. VI.A): visibility `v` (depolarisation) and detection efficiency `eta`. The no-click event is
binned into Bob's outcome 0 in the semidefinite programs and kept as a separate outcome in H(A|B). `G` is, in the
units of the code (the paper's values divided by n):
- the guessing probability of the target states, plus `dG`;
- with `dG=source` (visibility sweeps only), that of the depolarised states, the paper's G = 1 + v.

## Relaxations and partitions

In the relaxation strings, `B` and `E` are the PCPOP letters for Bob's projectors M_{0|y} and Eve's projectors P.

| `rel` | words W in M and P (App. E) |
|---|---|
| `lean:L` | all words of length at most L, e.g. `lean:3` |
| `lean:1+B*E` | 1, M, P, MP, equivalent to the tracial level 1+B·E (Table IV) |
| `lean:2+B*B*E` | all words of length at most two, plus the products M M′ P |
| `lean:2s` | all words of length at most two, except the products of Eve's projectors belonging to different key bases (a shorthand of this code) |
| `lean:2s+B*B*E`, `lean:2+B*B*E+B*B*B` | as above, plus the listed PCPOP terms |
| `tracial:LEVEL[:keepfree]` | the tracial relaxation (validation only, `scripts/validate.jl`) |

`lean` is the one-factor relaxation of Sec. V.B (Proposition 3). It imposes that the matrices [Tr(u† F v)] for u, v ∈ W be
positive semidefinite, for F = ρ_x and for F = σ − ρ_x/n.

The partitions (`grid`) are:
- `g26`: the 26-node partition of App. E, the default, used for all results except the comparisons of partitions;
- `uniform10`: t_k = k/10;
- `fine120`: 120 nodes, for convergence checks.

## Result files

A key-rate result file contains:
- the job parameters `job`;
- the protocol and the relaxation, with `nwords` words and `nmoments` moments;
- the honest value `G_honest` and the partition `grid`;
- `curve`, `bisection` and `bisection_minentropy`: the evaluated points.

Each point has the fields:
- the noise parameter `x`, and `v`, `η`, `G`;
- `H` = H(A|E), `Hmin`, and `Pg`, Eve's guessing probability of the key bit;
- `HAB` = H(A|B), the rates per key round `K` and `Kmin`, and the rates per signal `R` and `Rmin`;
- the values of the individual programs, `node_values`;
- `nonoptimal`, the termination and primal statuses of the entropy programs that were not solved to optimality;
- `pg_status`, the status of the guessing-probability program;
- `seconds`, the solution time of the point.

The points of `bisection_minentropy` only evaluate the min-entropy, so their `H`, `K` and `R` are `null`.

The file also records:
- `threshold` and `threshold_minentropy`, with their brackets;
- `no_key_at_ideal_point`, which is true if no evaluated point certified a positive rate;
- `bisection_stopped_vN` or `bisection_stopped_min`, present only if a program failed during a bisection (none did
  in the shipped files).

A `dG_scan.jl` result file lists the evaluated points (fields as above, plus `dG` and `minonly`, true for the
points that only evaluate the min-entropy), and the largest relaxations with a positive rate, `dG_max` (von Neumann)
and `dG_max_minentropy`.

A CHSH result file lists, for each value `S` of the CHSH value, the bounds `H`, `Hmin` and `Pg`. With them are the
device-independent values `H_DI` and `Pg_DI`, which are meaningful only up to S = 2√2 (beyond it they are clamped
to 1 and 1/2). A program that failed is recorded as `null`.

## Where the numbers of the paper come from

The files are in `results/` unless stated otherwise. `vN` is the von Neumann threshold (`threshold`) and `min` the
min-entropy threshold (`threshold_minentropy`).

| paper | result files |
|---|---|
| Table III, LB columns | `43__lean_3`, `44__lean_2+B_B_E`, `6sym__lean_2s+B_B_E` and `3st_ZBB__lean_3`, `v` and `eta`: `threshold`. The columns δG = 0, 0.04, 0.08, 0.12 are the files at `dG0.00`, `dG0.01`, `dG0.02`, `dG0.03` for (4,3) and (4,4), at `dG0.00`, `dG0.006667`, `dG0.013333`, `dG0.02` for (6,3), and at `dG0.00`, `dG0.013333`, `dG0.026667`, `dG0.04` for (3,3) |
| Fig. 2, curves | the same files at `dG0.00` (fields `x` and `K` of `curve` and `bisection`), and `results/curves/` for (4,4) and (6,3) |
| Fig. 3 and Sec. VI.B, thresholds versus δG | the files of (4,3), (4,4) and (6,3) at all values of `dG` (up to `dG0.17`), `v` and `eta`: `threshold`; the end points: `results/dG/`, field `dG_max` (and `dG_max_minentropy`) |
| Secs. IV.D and VI.B, error budget δG = 0.02 for the four-state protocols and 0.03 for (6,3) | `43__lean_3`, `44__lean_2+B_B_E` and `6sym__lean_2s+B_B_E` at `dG0.005` |
| Sec. VI.B, ideal rate, and rates of (4,3) and (3,3) at v = 0.95 and η = 0.95 | the Table III files at `dG0.00`, field `curve` |
| Sec. VI.C, rates at v = 0.99 | (4,4) and (6,3): `44__lean_2` and `6sym__lean_2s` at v = 0.99 |
| App. E, word sets | (4,4): `44__lean_2` vs `44__lean_2+B_B_E` (0.034 bits at v = 0.86), and `44__lean_2+B_B_E+B_B_B`; (6,3): `6sym__lean_2s` vs `6sym__lean_2s+B_B_E` (products M M′ P) and `6sym__lean_2` vs `6sym__lean_2s` (products of Eve's projectors of different bases); (4,3): `43__lean_4` vs `43__lean_3` |
| App. E, size of the programs | field `nmoments` |
| Sec. VI.C, (4,4) with M M′ M″ | `44__lean_2+B_B_E+B_B_B` |
| Sec. VI.A and App. F, noise in the source | `43__lean_3`, `43__lean_4`, `44__lean_2`, `63__lean_2s`, all `dGsource__v` |
| Table IV | (4,3): `lean_1+B_E` with `uniform10` and `g26`, `lean_2`, `lean_3` |
| Table V | `vN` and `min` of the main protocols (as for Table III at `dG0.00`) and of `63__lean_2s+B_B_E`, `BB84+1__lean_2+B_B_E`, `3st_ZB__lean_3`, `3st_ZXB__lean_3`, `3st_ZX__lean_3`; the entries marked † from `44__lean_2`, `6sym__lean_2s` and `63__lean_2s` |
| App. F, relaxed bound | `results/dG/` |
| App. F, tilted measurements | `43@φ__lean_3` |
| Sec. VI.D, Fig. 4 | `results/chsh/chsh__lean_4__G*__g26` and `chsh__lean_4__source__g26` |
| App. F, entropy versus the CHSH value | `results/chsh/`: `lean_3`, `lean_5`, and `lean_4` with `fine120` |
| App. F, validation of the relaxation | `scripts/validate.jl` |
| App. E, residual correction | `scripts/residual_check.jl` |

Some numbers of the paper come from single points or from post-processing of the files, and are not stored as
such.
- **Partitions (App. E).** The comparisons with the 120-node and uniform partitions near the thresholds can be
  repeated with `scripts/point.jl`, for instance for (4,3) at v = 0.852 with `grid=g26`, `fine120` and `uniform10`.
- **Timings and length-four words (App. E).** The times per program were measured with `scripts/point.jl` on a
  single job. The field `seconds` of the result files is the time of a whole point (27 programs), measured with
  several jobs running in parallel. The effect of words of length four for (4,3) can also be checked with
  `scripts/point.jl`.
- **Thresholds from F(S, G) (Sec. VI.D).** They follow from the CHSH scans. The rate of (4,3) from the CHSH value
  alone is F(S, G) − H(A|B). For a visibility v, S = 2√2 v and H(A|B) = h((1 − v)/2); for an efficiency η, S = 2√2 η
  and H(A|B) = 1 − η. F is non-decreasing in S and non-increasing in G, so its values on the grid of S bound it
  rigorously between grid points.

## Numerical notes

- **Values used.** The value of each program is the dual objective at MOSEK's solution. It is used only when MOSEK
  reports a feasible or nearly feasible point, and is then a lower bound up to the solver tolerances.
  `scripts/residual_check.jl` removes that dependence at given points.
- **Honest G.** The honest value of G is computed by semidefinite programming. It is below 2/n by at most 10⁻⁸,
  which makes the assumption very slightly stronger than the honest one. The effect on the rates is negligible.
- **Thresholds.** They are bisected to 5·10⁻⁴ after a sweep in steps of 0.01, and assume that the rate decreases
  with the noise. The jobs with `xmax` < 1 start their sweep just above the von Neumann threshold. A min-entropy
  threshold is found only when the sweep crosses it; otherwise it is `null`.
- **Largest δG.** `dG_scan.jl` treats rates below `Ktol` = 10⁻⁴ bits as zero: where the exact value of a
  relaxation is zero, MOSEK returns about 10⁻⁷ bits. For the same reason, threshold jobs at a δG beyond the end
  point return a meaningless `threshold` close to 1; such runs are not shipped.
- **3st_ZX.** This protocol admits an exact classical attack. Its computed entropy at the ideal point is about
  10⁻⁶ instead of zero, so the `threshold` recorded for it (≈ 0.9998) is not meaningful.
- **Output files.** The name of a result file does not contain `xmax`, `xmin` or `step`. Runs of the same job with
  different ranges must therefore be written to different directories (`out=`), as for `results/curves/`.
- **Reproducibility.** In our tests, re-running a job on another machine gave the same thresholds, and point values
  that agree to about 10⁻⁸. The exceptions are ill-conditioned programs, for which MOSEK reports slow progress: at the ideal
  point, and for the small word set `lean:1+B*E` near it, values can differ by up to about 10⁻⁴.

## Provenance of the results

The results were computed on 28–29 September 2026, those of Fig. 3 on 30 September–1 October 2026, and the files
at `dG0.006667`, `dG0.013333`, `dG0.026667` and `dG0.04` (the columns of Table III for (6,3) and (3,3)) on 7–8 October
2026, with:
- Julia 1.11.7;
- MOSEK 11.2.4;
- PCPOP.jl at commit `813b644`;
- a preliminary version of this code for the results of 28–29 September, and the present version for the later
  ones (the files at `dG0.005`, at `dG0.04`–`dG0.17` and at the fractional values of `dG`, and `results/dG/`).

That version read the coefficients of the entropy bound and the Shannon entropies from `utils.jl` of A. Mishra's
[SDIQKD](https://github.com/abh1mishra/SDIQKD) repository, at commit `a5fa6f7`. `src/utils.jl` computes them with
the same expressions, and the values are bit-for-bit identical.

The metadata of the result files come from that version:
- the key-rate files record `mosek`, `pcpop_commit` and `sdiqkd_commit`, whereas files written by the present
  version record `versions`;
- the CHSH files have no version fields.

The protocol descriptions and the label of the uniform partition were updated to the names used here. The numerical
content of the files is unchanged.
