# Provenance of the test fixtures

The five files in this directory are excerpts of the output of one
Nbody6++GPU-Beijing run, kept as fixtures for the reader tests
(`test/test_io.jl`, "Real-output fixtures"). They exercise the readers on the
engine's real column layouts and number formats, which synthetic files cannot.

| File | Content | Excerpt |
|---|---|---|
| `out1000` | captured stdout: input echo, PHYSICAL SCALING line, ADJUST records | first 370 lines, seven adjustments |
| `lagr.7` | Lagrangian radii table with its header line | header and the first epochs |
| `esc.11` | escaper records with the direction columns | 88 escapers at t = 7 to 7.5 N-body units |
| `sev.83_0` | single-star evolution snapshot at t = 0, 15-field `hrplot.F` layout | first 40 stars |
| `bev.82_0` | regularised-binary snapshot at t = 0, 32-field `hrplot.F` layout | first 40 pairs |

Origin. Produced on 2026-07-31 with the package's own install and run
pipeline (package commit `757d1bd`), engine Nbody6PPGPU-beijing at commit
`618d7a4` (upstream v2026.07), the commit `install.ref` defaults to. The
generating parameters are those the engine echoes at the top of `out1000`:
N = 3948, NRAND = 3948, NNBOPT = 63, ETAI = ETAR = 0.02, RS0 = 0.5,
DTADJ = DELTAT = 0.5, TCRIT = 10, QE = 1.0, RBAR = 5.5 pc, ZMBAR = 7.0,
with the KZ options printed below them; the PHYSICAL SCALING line gives
R* = 5.504 pc, M* = 27 696 M☉, T* = 1.159 Myr. The run directory itself was
not retained (run data of that date were purged), so the fixtures are
reproducible only in layout, not bit for bit: a new run with these
parameters and this engine commit writes files of the same structure with
different realisations.

Licence. The files are simulation output produced by the package author and
are released with the package under its MIT licence.
