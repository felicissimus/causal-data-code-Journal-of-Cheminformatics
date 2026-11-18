"""Output object: retain DoWhy results and resolve them into saved outputs."""
from datetime import datetime
import json
import os
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd


class Output:
    """Hold the original estimand, estimate, placebo, and random-refuter objects."""

    def __init__(self, estimand, estimate, placebo, random_refute):
        self.estimand = estimand
        self.estimate = estimate
        self.placebo = placebo
        self.random_refute = random_refute

    @staticmethod
    def prepare_output(settings):
        """Create a fresh run folder and record the resolved input and DAG."""
        output_dir = Path(settings["output_dir"]) / datetime.now().strftime("run_%Y%m%d_%H%M%S_%f")
        output_dir.mkdir(parents=True, exist_ok=False)
        graph = nx.DiGraph()
        graph.add_nodes_from(settings["graph"]["nodes"])
        graph.add_edges_from(settings["graph"]["edges"])
        nx.write_graphml(graph, output_dir / "DAG.graphml")
        pd.DataFrame(graph.edges, columns=["cause", "effect"]).to_csv(output_dir / "DAG_edges.csv", index=False)
        recorded = dict(settings,
                        adjustment_selection="Automatic DoWhy default backdoor set; no overrides",
                        python_hash_seed=os.environ.get("PYTHONHASHSEED"),
                        interpretation="Mean marginal slopes; no SE/CI. Refuters do not account for trajectory autocorrelation.")
        (output_dir / "run_settings.json").write_text(json.dumps(recorded, indent=2) + "\n", encoding="utf-8")
        return output_dir

    def resolve_output(self, dataset, treatment, data, settings, output_dir):
        """Save CATEs and reports, append one summary row, and return that row."""
        estimand = self.estimand
        estimate = self.estimate
        placebo = self.placebo
        random_refute = self.random_refute
        controls = estimand.get_backdoor_variables()
        modifiers = estimate.params["effect_modifiers"]

        # LinearDML contrasts are linear in treatment; report slopes per unit.
        contrast = estimate.treatment_value - estimate.control_value
        contrast = float(np.asarray(contrast).item())
        effect = float(np.asarray(estimate.value).item()) / contrast
        cate = np.asarray(estimate.cate_estimates, dtype=float).reshape(-1) / contrast
        if cate.size == 1:
            cate = np.full(len(data), cate.item())
        if cate.size != len(data) or not np.isfinite(cate).all():
            raise ValueError("Invalid row-wise CATE output.")
        np.testing.assert_allclose(cate.mean(), effect, rtol=1e-8, atol=1e-10)
        outcome = settings["outcome"]
        units = settings["units"]
        row = dict(system=dataset["name"], temperature_K=dataset.get("temperature_K"),
                   treatment=treatment, outcome=outcome,
                   effect_modifiers_X=", ".join(modifiers), adjustment_W=", ".join(controls),
                   mean_marginal_effect=effect, effect_unit=f"({units[outcome]}) per {units[treatment]}",
                   placebo_effect=float(np.asarray(placebo.new_effect).item()) / contrast,
                   placebo_p_value=float(np.asarray(placebo.refutation_result["p_value"]).item()),
                   random_cause_effect=float(np.asarray(random_refute.new_effect).item()) / contrast,
                   random_cause_p_value=float(np.asarray(random_refute.refutation_result["p_value"]).item()))

        output_dir = Path(output_dir)
        temperature = dataset.get("temperature_K")
        label = dataset["name"] if temperature is None else f"{dataset['name']}_{temperature}K"
        prefix = output_dir / f"{label}_{treatment}"
        cate_file = Path(f"{prefix}_CATE.csv")
        report_file = Path(f"{prefix}_report.txt")
        if cate_file.exists() or report_file.exists():
            raise FileExistsError(f"Results already exist for {label}, {treatment}.")
        pd.DataFrame({"input_row_zero_based": np.arange(len(data)),
                      "CATE_per_treatment_unit": cate}).to_csv(cate_file, index=False)
        reports = [f"X = {modifiers}\nW = {controls}", str(estimand), str(estimate),
                   str(placebo), str(random_refute)]
        report_file.write_text("\n\n".join(reports) + "\n", encoding="utf-8")
        summary_file = output_dir / "causal_effects_and_refutations.csv"
        pd.DataFrame([row]).to_csv(summary_file, mode="a", header=not summary_file.exists(), index=False)
        print(f"{dataset['name']}, {treatment}: X={modifiers}, W={controls}", flush=True)
        print(f"  Mean marginal effect = {effect:+.8g}", flush=True)
        return row
