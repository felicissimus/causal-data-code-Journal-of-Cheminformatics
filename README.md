# Causal analysis for Diels-Alder reactions

## Scope and contents

This package accompanies **Causal Associations in Diels-Alder Reactions** for the Journal of Cheminformatics. **Additional file 2** is
`causal-data-code-Journal-of-Cheminformatics.zip`, containing the
top-level folder `causal-JACS-submit-data-code/`. Additional file 1 is the
separate supporting-information document.

The package analyzes precomputed descriptors. It contains no raw trajectories, coordinates, electronic-structure outputs, or upstream descriptor-generation
scripts. It does not run molecular dynamics or generate manuscript figures.

| File | Purpose |
| --- | --- |
| `README.md` | Package documentation |
| `input.json` | Datasets, variables, units, DAG, models, refuters, output location |
| `analysis_input.py` | `Input`: configuration and CSV loading/validation |
| `causal_inference.py` | Entry point; identification, estimation, refutation |
| `analysis_output.py` | `Output`: result objects and file export |
| `da-01-BU.csv`, `da-02-BU.csv`, `da-04-BU.csv` | Three descriptor datasets |
| `LICENSE-CODE.txt`, `LICENSE-DATA.txt` | License notices |
| `SHA256SUMS` | Release file-integrity manifest |

## Data and preprocessing

Each CSV has **4,402 rows and eight columns**. All three datasets have `temperature_K=300` in the configuration, not in the CSVs. No trajectory ID,
timestamp, or time-step column is supplied.

| Column | Role | Configured unit |
| --- | --- | --- |
| `Interaction_Energy` | Outcome | `eV/atom` |
| `omega_k_R1` | Treatment | `eV` |
| `Charge_Transfer` | Treatment | `e` |
| `alpha_C_C` | Treatment | `angstrom` |
| `beta_C_C` | Treatment | `angstrom` |
| `torsion` | Treatment | `degree` |
| `Pyramidal_R1` | Treatment | `degree` |
| `Delta_phi` | Treatment | `degree` |

The loader uses pandas `float_precision="round_trip"`, selects `graph.nodes`
in configured order, ignores extra columns, and retains row order. It requires
at least `max(4, cv)` rows and rejects missing columns, nonnumeric/nonfinite
values, and constant features. It performs no imputation, row filtering,
rescaling, standardization, angle wrapping, trajectory grouping, or descriptor
recalculation. Unit strings label outputs; no unit conversion is performed.

## Environment and execution

Dependencies: DoWhy, EconML, LightGBM, NetworkX, NumPy, pandas, and scikit-learn.
The original README reported **Python 3.10, DoWhy 0.12, EconML 0.16.0, and
LightGBM 4.6.0** as its tested environment. These are historical reports, not
newly verified versions. Other dependency versions were unspecified; no
environment lockfile or verified numerical reference output is supplied.

From the extracted package directory in PowerShell, with Python 3.10 selected as `python`:

```powershell
python --version
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install "dowhy==0.12" "econml==0.16.0" "lightgbm==4.6.0" networkx numpy pandas scikit-learn
$env:PYTHONHASHSEED = "42"
.\.venv\Scripts\python.exe -B .\causal_inference.py
.\.venv\Scripts\python.exe -m pip freeze
```

Unpinned dependencies remain resolver-selected; the installation command does
not reconstruct the complete historical environment. Retain `pip freeze`
output for the environment used. These commands are instructions, not a claim
of a new successful run.

The entry point reads `input.json` beside itself; there is no alternative-config
command-line option. Dataset/output paths resolve relative to the JSON file.
The `Input` API can load another JSON path.

## Configuration, workflow, and seeds

The supplied `input.json` is the authoritative DAG/configuration source:
eight nodes, 12 directed edges, three datasets, and seven treatments in the
orders above. The default outcome is `Interaction_Energy`; output goes to
`results`. Settings described here are those currently configured in code.

- Estimator: `backdoor.econml.dml.LinearDML`; `discrete_treatment=false`,
  `cv=2`, `random_state=42`; `fit_params={}`.
- Both separate LightGBM nuisance models: `n_estimators=50`, `max_depth=4`,
  `learning_rate=0.1`, `random_state=42`, `n_jobs=1`,
  `force_row_wise=true`, `verbose=-1`.
- Estimation: `control_value=0`, `treatment_value=1`, `target_units="ate"`.
- Placebo: `placebo_treatment_refuter`, `placebo_type="permute"`;
  random-cause refuter: `random_common_cause`. Each uses
  `num_simulations=100` and `n_jobs=1`.
- Neither refuter specifies `random_state`; installed-library defaults apply.
  No global Python/NumPy RNG seed is set. `PYTHONHASHSEED` is separately set
  before launch and recorded; it does not seed the refuters.
- The entry point defaults `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, and
  `MKL_NUM_THREADS` to `1` only when unset.

Configuration validation checks DAG structure, variables, units, dataset paths,
names, continuous-treatment/all-row settings, finite nonzero contrast,
`cv>=2`, and at least two simulations per refuter.

Each dataset is loaded once. For each treatment, a fresh `CausalModel`
identifies a backdoor estimand using DoWhy defaults, estimates with fresh
nuisance learners, then runs the placebo and random-cause refuters. No manual
adjustment or effect-modifier list is passed. The default loop has 21 cases.

## Outputs

Each run creates `results/run_YYYYMMDD_HHMMSS_microseconds/` containing:

- `run_settings.json`: resolved configuration/paths, hash seed, adjustment-selection
  description, and interpretation text; not installed dependency versions.
- `DAG.graphml` and `DAG_edges.csv` (`cause,effect`).
- `causal_effects_and_refutations.csv`: one row per saved case, 21 on full success.
- Per-case `<dataset>_300K_<treatment>_CATE.csv` and matching `_report.txt`.

Summary columns, in order:
```text
system, temperature_K, treatment, outcome, effect_modifiers_X, adjustment_W,
mean_marginal_effect, effect_unit, placebo_effect, placebo_p_value,
random_cause_effect, random_cause_p_value
```
CATE columns: `input_row_zero_based`, `CATE_per_treatment_unit`; 4,402 rows
per supplied dataset. A scalar CATE is repeated across rows. The exporter
checks finite CATEs and agreement of their mean with the effect
(`rtol=1e-8`, `atol=1e-10`).

CSV effects are divided by the treatment-control contrast; reports retain
original DoWhy result representations and list X/W. Units are
`(<outcome unit>) per <treatment unit>`. CSV schemas contain no SE/CI columns.
Errors stop execution; partial outputs may remain. Existing per-case files
are not overwritten.

## Quick checks

These checks compile syntax and use the actual input loader without fitting
estimators or creating result directories. Run from the package directory:
```powershell
.\.venv\Scripts\python.exe -B -c "from pathlib import Path; [compile(Path(f).read_text(encoding='utf-8-sig'), f, 'exec') for f in ('analysis_input.py','analysis_output.py','causal_inference.py')]; print('Syntax OK')"
.\.venv\Scripts\python.exe -B -c "from analysis_input import Input; x=Input('input.json'); [print(d['name'], x.read_data(d['path']).shape) for d in x.settings['datasets']]"
```
Expected loader shapes: `(4402, 8)` for each dataset. These checks do not
establish estimator compatibility or numerical agreement.

## Licenses

`LICENSE-CODE.txt` supplies the **MIT** license for the three Python scripts,
`input.json`, and README. `LICENSE-DATA.txt` supplies the **CC BY 4.0** notice
and legal-terms URL for the three CSVs. Follow those files for attribution and
license terms. Third-party dependencies retain their own licenses.
