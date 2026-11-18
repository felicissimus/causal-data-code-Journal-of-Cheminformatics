"""Input object: load and validate the configuration and data."""
import json
import math
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd


class Input:
    """Own the analysis settings and the methods that read and use them."""

    def __init__(self, filename):
        self.filename = Path(filename).expanduser().resolve()
        self.read_input()

    def read_input(self):
        """Load and validate one continuous-treatment LinearDML configuration."""
        filename = self.filename
        settings = json.loads(filename.read_text(encoding="utf-8-sig"))
        required = {"outcome", "treatments", "datasets", "units", "graph", "dml", "refuters", "output_dir"}
        missing = required - settings.keys()
        if missing:
            raise ValueError(f"Missing input fields: {sorted(missing)}")

        nodes = settings["graph"]["nodes"]
        treatments = settings["treatments"]
        if not nodes or len(nodes) != len(set(nodes)):
            raise ValueError("Graph nodes must be nonempty and unique.")
        if not treatments or len(treatments) != len(set(treatments)):
            raise ValueError("Treatments must be nonempty and unique.")
        if settings["outcome"] not in nodes or not set(treatments).issubset(nodes):
            raise ValueError("The outcome and all treatments must appear in graph.nodes.")
        if settings["outcome"] in treatments:
            raise ValueError("The outcome cannot also be a treatment.")
        if not {settings["outcome"], *treatments}.issubset(settings["units"]):
            raise ValueError("Provide units for the outcome and every treatment.")
        graph = nx.DiGraph()
        graph.add_nodes_from(nodes)
        for edge in settings["graph"]["edges"]:
            if len(edge) != 2 or not set(edge).issubset(nodes):
                raise ValueError(f"Invalid graph edge: {edge}")
            graph.add_edge(*edge)
        if not nx.is_directed_acyclic_graph(graph):
            raise ValueError("The graph must not contain directed cycles.")

        dml = settings["dml"]
        if dml["method_name"] != "backdoor.econml.dml.LinearDML":
            raise ValueError("This workflow supports backdoor.econml.dml.LinearDML.")
        if dml["init_params"].get("discrete_treatment", False):
            raise ValueError("This workflow reports slopes for continuous treatments only.")
        if {"model_y", "model_t"} & dml["init_params"].keys():
            raise ValueError("Place nuisance-model parameters in dml.model_y and dml.model_t.")
        cv = dml["init_params"].get("cv", 2)
        if type(cv) is not int or cv < 2:
            raise ValueError("dml.init_params.cv must be an integer of at least two.")
        estimate = dml["estimate_params"]
        if estimate.get("target_units", "ate") != "ate":
            raise ValueError("This workflow summarizes all input rows; use target_units='ate'.")
        contrast = estimate.get("treatment_value", 1) - estimate.get("control_value", 0)
        if not math.isfinite(contrast) or contrast == 0:
            raise ValueError("The treatment-control contrast must be finite and nonzero.")
        for key, method in [("placebo", "placebo_treatment_refuter"), ("random_refute", "random_common_cause")]:
            refuter = settings["refuters"][key]
            if refuter["method_name"] != method:
                raise ValueError(f"refuters.{key}.method_name must be {method}.")
            if type(refuter["num_simulations"]) is not int or refuter["num_simulations"] < 2:
                raise ValueError("Use at least two simulations per refuter.")

        datasets = settings["datasets"]
        if not datasets or len({item["name"] for item in datasets}) != len(datasets):
            raise ValueError("Provide at least one dataset, with unique names.")
        for name in [*treatments, *(item["name"] for item in datasets)]:
            if not name or name in (".", "..") or any(c in name for c in '<>:"/\\|?*'):
                raise ValueError(f"Name cannot be used in output filenames: {name!r}")
        for item in datasets:
            temperature = item.get("temperature_K")
            if temperature is not None and (type(temperature) not in (int, float) or not math.isfinite(temperature)):
                raise ValueError("temperature_K must be a finite number or null.")
            item["path"] = str((filename.parent / Path(item["path"]).expanduser()).resolve())
            if not Path(item["path"]).is_file():
                raise FileNotFoundError(item["path"])
        settings["output_dir"] = str((filename.parent / Path(settings["output_dir"]).expanduser()).resolve())
        settings["input_file"] = str(filename)
        self.settings = settings
        return settings

    def read_data(self, filename):
        """Select graph columns in their configured order and validate the data."""
        settings = self.settings
        data = pd.read_csv(filename, float_precision="round_trip")
        columns = settings["graph"]["nodes"]
        missing = set(columns) - set(data.columns)
        if missing:
            raise ValueError(f"{filename}: missing graph columns {sorted(missing)}")
        data = data.loc[:, columns]
        minimum_rows = max(4, settings["dml"]["init_params"].get("cv", 2))
        if len(data) < minimum_rows:
            raise ValueError(f"{filename}: need at least {minimum_rows} rows for this configuration.")
        if not np.isfinite(data.to_numpy(dtype=float)).all() or data.nunique().le(1).any():
            raise ValueError(f"{filename}: missing, nonfinite, or constant feature values.")
        return data

