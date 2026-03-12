import os
import sys
import yaml
import pandas as pd

os.makedirs("logs", exist_ok=True)
os.makedirs("plots", exist_ok=True)

# checking env yaml files
assert os.path.isfile("envs/main.yaml"), "Missing envs/main.yaml (conda env)"

# checking config exists
if not os.path.isfile("config.yaml"):
    raise FileNotFoundError("Missing config.yaml")

try:
    with open("config.yaml", "r") as handle:
        data = yaml.safe_load(handle)
except yaml.YAMLError as e:
    raise ValueError(f"config.yaml is not valid YAML: {e}")

if not isinstance(data, dict):
    raise ValueError("config.yaml is empty or has wrong structure (expected key:value pairs).")

# check required keys
required = ["ref_organism", "ref_organism_annotation", "ref_enh", "sample_list"]
missing_k = [k for k in required if k not in data]
if missing_k:
    raise KeyError(f"Missing keys in config.yaml: {missing_k}")
missing_v = [k for k in required if not data.get(k)]
if missing_v:
    raise ValueError(f"Empty values for required keys in config.yaml: {missing_v}")

print("You are going to run the pipeline with following input parameters:")
print(f"Reference organism assembly: {data['ref_organism']}")
print(f"Reference annotation for gene-enhancer pairing: {data['ref_organism_annotation']}")

# check ref enhancers file
ref_enh_path = data["ref_enh"]
if not os.path.isfile(ref_enh_path):
    raise FileNotFoundError(f"Reference enhancers file not found: {ref_enh_path}")

ref_enhancers = pd.read_csv(ref_enh_path, sep="\t", header=None, usecols=[0, 1, 2])
if ref_enhancers.shape[1] != 3:
    raise ValueError("Reference enhancer regions are not BED-like (expected 3 columns: chr, start, end, without a header).")

# check species list
sample_list_path = data["sample_list"]
if not os.path.isfile(sample_list_path):
    raise FileNotFoundError(f"sample_list file not found: {sample_list_path}")

species_df = pd.read_csv(sample_list_path, sep=";", header=0)
expected_cols = {"Species", "Latin", "Assembly", "Annotation"}
if not expected_cols.issubset(set(species_df.columns)):
    raise ValueError(
        f"Species list columns look wrong. Expected {expected_cols}, got {list(species_df.columns)}"
    )

print("[OK] Preflight checks passed.")