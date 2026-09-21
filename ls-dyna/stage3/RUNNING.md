# Running Stage 3

Nine decks, one folder each. The folders matter: LS-DYNA writes `d3plot`,
`glstat`, `spcforc`, `nodout` and `rcforc` into whatever directory it is
started in, always under those same names, so nine runs in one folder would
overwrite each other. The extractor looks for the same folder names.

```
ls-dyna/stage3/
    run_all.bat                       Windows
    run_all.sh                        Linux / macOS
    stage3_lance_only/
        stage3_lance_only.k
    stage3_extract_mu000/
        stage3_extract_mu000.k
    ...
```

## Run them

Start the solver **inside** each folder. In LS-Run or LS-PrePost, set the
working directory to the case folder and open the `.k` file there.

From a command line, both loops do it in one go:

```
run_all.bat "C:\path\to\lsdyna_dp.exe" 4
```

```
./run_all.sh /path/to/ls-dyna 4
```

Second argument is the core count. Every deck is small — 1200 to 2064
elements — and finishes in a couple of minutes. `insert_mu020` is the long
one: 31 ms of simulated time against 5.5 ms for the rest.

A note on output size. `spcforc` is written every 0.5 microseconds, which is
about 30 MB per extraction run and 90 MB for the insertion one. That is
deliberate and it is half the fix for the first set of runs — see below. Keep
the run folders until the post-processor's section 0 comes back clean; after
that nothing downstream needs them.

## Then

```
python scripts/extract_stage3.py ls-dyna/stage3
python scripts/postprocess_stage3.py
```

The first walks the folders, reads `spcforc`, `nodout` and `glstat` from each,
and writes `results/stage3_results.csv` plus the curve history. The second
compares them against the closed form and draws `results/figures/`.

**Read section 0 of the post-processor first.** It says whether each run's
force columns may be quoted, against three gates:

| Gate | Threshold | If it fires |
|---|---|---|
| samples per cycle of the ring | ≥ 8, where the ripple is over 5 % of the level | halve `force` in `make_stage3_models.py` and rerun that deck |
| movement across averaging windows of 50–400 µs | ≤ 1.5° | no rerun helps; the estimator is at its limit for that run |
| force route against shape route, worst station | ≤ 1.5° | same |

Only the first gate is fixed by running again. The other two are properties of
how noisy that particular run's contact force is, and a run that fails them
still reports its displacement results, which do not touch the force output.

Pass the extractor whatever folder actually holds the runs — if the decks were
copied somewhere else, give it that path instead.

## The nine

| folder | what it is for |
|---|---|
| `stage3_lance_only` | regression. Tip lifted 0.6 mm, nothing else touching the lance. The root reaction must come back at the Stage 1 value plus the 1.06 % the tooth stiffens the beam by. If this drifts, nothing below it means anything |
| `stage3_extract_mu000` | frictionless extraction. W/P is the tangent of the effective face angle with no friction term in it, so it is the one measurement in the stage that assumes nothing |
| `stage3_extract_mu020` | the baseline. µ = 0.20 |
| `stage3_extract_mu030` | µ = 0.30. The wedge relation has to hold at both |
| `stage3_insert_mu020` | the same terminal driven the other way, over the 30° lead-in |
| `stage3_extract_tpa_g010` | TPA at 0.10 mm clearance. Must block |
| `stage3_extract_tpa_g085` | TPA at 0.85 mm, above the 0.757 mm limit. Must not block, and must not change the force either — a contact that is defined but never reached has to cost exactly nothing |
| `stage3_extract_mu020_fine` | baseline with the element over the tooth halved |
| `stage3_extract_mu020_slow` | baseline with the ramp doubled |

## If a run does not finish

The four failures this model already went through, and what they looked like:

- **Runs, writes zeros.** Card-perfect deck, no force anywhere. Two causes
  found: a part ID used twice, and `*MAT_RIGID` CON1 set to 6 instead of 5,
  which constrains the direction the terminal is driven in. The generator now
  refuses both before it writes the file.
- **Terminates on element distortion during insertion.** One-way contact let
  the terminal's leading corner sit inside a tooth face between nodes. Fixed
  with a two-way contact and a chamfer on the nose.
- **Blocked run keeps driving.** The terminal is rigid and displacement
  controlled, so once the TPA stops the lance the stroke has nowhere to go.
  The blocked decks are stroked to the blocking travel plus 0.05 mm.

## The fifth failure, and why the output rates are what they are

The first complete set of nine ran cleanly. Energy balance 0.99983 to 1.00000,
hourglass exactly zero, every run finished its full stroke, the TPA control
matched the baseline to every decimal place. Every force number in it was
still worthless.

Everything was written every 10 microseconds. The contact force rings far
above 50 kHz, so at that rate it was aliased — folded down into the data,
where no filter can separate it from the answer. It did not look like a
failure. It looked like a noisy curve.

Two things gave it away. The measured angle moved 12 degrees depending on how
wide a filter was used, while the same angle taken from the deflected shape
moved less than one. And the apparent ringing frequency came out 44 kHz in one
run and 16 kHz in the same model run at half the loading rate, which no real
structural mode can do — only an alias changes frequency when you change the
sample rate.

So forces are now written at 2 MHz and displacements stay at the old rate,
because a displacement is the double integral of the acceleration and the ring
is already small in it. That is also why the first set's rotation measurement
survived intact while its forces did not.

**Resolving the ring turned out to be necessary and not sufficient.** At
0.5 µs the ring resolves at 12 to 19 samples per cycle and the answer still
moved several degrees with the filter. The ripple is 55 to 72 % of the force's
own level, and W and P are two projections of the same noisy contact impulse,
so their ratio *at an instant* scatters by degrees however it is filtered.

Stage 2 had already settled this and its rule is the one that applies: take
the ratio of the MEANS over a window, not the mean of a ratio, and not a
point. Applied at stations through the stroke, each compared against the
rotation over the same window, the two routes agree to under 1 degree the
whole way down the sweep.

Three gates now guard it, and the lesson is the gates rather than any one
number: this failed silently three times, so something has to look for it.
