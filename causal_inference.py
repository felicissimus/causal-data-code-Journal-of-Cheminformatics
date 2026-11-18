"""Run the configured workflow: python causal_inference.py"""
import os
from pathlib import Path

for option in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(option, "1")

import networkx as nx
from dowhy import CausalModel
from lightgbm import LGBMRegressor
from analysis_input import Input
from analysis_output import Output

INPUT_FILE = Path(__file__).with_name("input.json")


def causal_analysis(model, inputs):
    """Identify, estimate, and refute one effect using an Input object.

    Return an Output object holding all four original DoWhy results.
    """
    settings = inputs.settings
    # Select the adjustment set automatically from the model's supplied DAG.
    estimand = model.identify_effect()
    if estimand.estimands.get("backdoor") is None:
        raise ValueError("No backdoor effect is identified by this DAG.")

    # Create fresh nuisance learners for each treatment/system.
    dml = settings["dml"]
    init_params = dict(
        dml["init_params"],
        model_y=LGBMRegressor(**dml["model_y"]),
        model_t=LGBMRegressor(**dml["model_t"]),
    )
    estimate = model.estimate_effect(
        estimand,
        method_name=dml["method_name"],
        method_params={"init_params": init_params, "fit_params": dml["fit_params"]},
        **dml["estimate_params"],
    )

    # Omitted refuter seeds use DoWhy defaults, separately from the DML seed.
    placebo = model.refute_estimate(
        estimand, estimate, **settings["refuters"]["placebo"],
    )
    random_refute = model.refute_estimate(
        estimand, estimate, **settings["refuters"]["random_refute"],
    )
    return Output(estimand, estimate, placebo, random_refute)


if __name__ == "__main__":
    inputs = Input(INPUT_FILE)
    output_dir = Output.prepare_output(inputs.settings)
    graph = nx.DiGraph()
    graph.add_nodes_from(inputs.settings["graph"]["nodes"])
    graph.add_edges_from(inputs.settings["graph"]["edges"])

    for dataset in inputs.settings["datasets"]:
        data = inputs.read_data(dataset["path"])
        for treatment in inputs.settings["treatments"]:
            model = CausalModel(
                data=data,
                graph=graph.copy(),
                treatment=treatment,
                outcome=inputs.settings["outcome"],
            )
            output = causal_analysis(model, inputs)
            output.resolve_output(dataset, treatment, data, inputs.settings, output_dir)

    print(f"Results saved to {output_dir}")
