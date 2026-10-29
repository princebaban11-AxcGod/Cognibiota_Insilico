"""
Cognibiota - in silico screening sketch (COBRApy)

What it does, per strain model:
  1. Loads the metabolic model (SBML file, e.g. from AGORA/VMH or CarveMe).
  2. Runs FBA under a "relaxed" and a "stressed" condition  -> growth rate.
  3. Fixes growth at >= 90% of its optimum and runs FVA on the product
     exchanges (GABA, acetate, propionate, butyrate) -> max possible output.
  4. Computes resilience (stressed/relaxed growth) and retention
     (stressed/relaxed product output) and ranks the strains.

Install:   pip install cobra pandas
Run demo:  python cognibiota_insilico.py          (uses the built-in E. coli
           core model, only to prove the pipeline runs - it makes no GABA)
Real use:  put your SBML files in STRAINS below.

IMPORTANT: stress values are illustrative PROXIES, not measured numbers.
Reaction IDs differ between AGORA and BiGG/CarveMe, so the script tries
several common IDs; check the printed warnings for anything it could not find.
"""

import math

import cobra
import pandas as pd
from cobra.flux_analysis import flux_variability_analysis

# ----------------------------------------------------------------------
# 1. Your strains: name -> path to SBML file (None = built-in demo model)
# ----------------------------------------------------------------------
STRAINS = {
    "L_brevis_ATCC367": "Models/Lactobacillus_brevis_ATCC_367_AGORA2_v2.01.xml",
}

# ----------------------------------------------------------------------
# 2. Conditions (proxies). Uptakes are in mmol/gDW/h.
#    Stress = less carbon/glutamate + extra ATP spent on maintenance.
# ----------------------------------------------------------------------
CONDITIONS = {
    "relaxed":  {"glc_uptake": 10.0, "glu_uptake": 10.0, "extra_atpm": 0.0},
    "stressed": {"glc_uptake": 3.0,  "glu_uptake": 3.0,  "extra_atpm": 5.0},
}

# Candidate reaction IDs (AGORA style first, then BiGG style)
GLUCOSE_IDS = ["EX_glc_D_e_", "EX_glc_D(e)", "EX_glc__D_e", "EX_glc_D_e"]
GLUTAMATE_IDS = ["EX_glu_L_e_", "EX_glu_L(e)", "EX_glu__L_e", "EX_glu_L_e"]
ATPM_IDS = ["DM_atp_c_", "ATPM", "DM_atp(c)"]

# Products to track: label -> candidate exchange IDs
PRODUCTS = {
    "GABA":       ["EX_4abut_e_", "EX_4abut(e)", "EX_4abut_e", "EX_gaba_e"],
    "acetate":    ["EX_ac_e_", "EX_ac(e)", "EX_ac_e"],
    "propionate": ["EX_ppa_e_", "EX_ppa(e)", "EX_ppa_e"],
    "butyrate":   ["EX_but_e_", "EX_but(e)", "EX_but_e"],
}

GROWTH_FRACTION = 0.9  # keep growth >= 90% of optimum while maximising product
UPTAKE_CAP = 5.0   # max uptake for any nutrient not set explicitly (mmol/gDW/h)

def first_present(model, ids):
    """Return the first reaction from `ids` that exists in the model, else None."""
    for rid in ids:
        if rid in model.reactions:
            return model.reactions.get_by_id(rid)
    return None


def apply_conditions(model, cond):
    """Change bounds (call inside `with model:` so it is undone afterwards)."""
    # Limit every nutrient so nothing is unlimited
    for r in model.exchanges:
        if r.lower_bound < -UPTAKE_CAP:
            r.lower_bound = -UPTAKE_CAP

    glc = first_present(model, GLUCOSE_IDS)
    glu = first_present(model, GLUTAMATE_IDS)
    atpm = first_present(model, ATPM_IDS)

    if glc is None:
        print("[warn] glucose exchange not found")
    else:
        glc.lower_bound = -cond["glc_uptake"]      # negative = uptake
    if glu is None:
        print("[warn] glutamate exchange not found")
    else:
        glu.lower_bound = -cond["glu_uptake"]
    if atpm is not None:
        atpm.lower_bound = atpm.lower_bound + cond["extra_atpm"]

def simulate(model, cond):
    """FBA for growth, then FVA for max product output at >=90% growth."""
    result = {}
    with model:
        apply_conditions(model, cond)

        growth = model.slim_optimize()
        if math.isnan(growth) or growth < 1e-6:
            result["growth"] = 0.0
            return result          # no growth -> products not meaningful
        result["growth"] = growth

        found = {label: first_present(model, ids) for label, ids in PRODUCTS.items()}
        rxns = [r for r in found.values() if r is not None]
        if rxns:
            fva = flux_variability_analysis(
                model, reaction_list=rxns, fraction_of_optimum=GROWTH_FRACTION
            )
            for label, rxn in found.items():
                if rxn is not None:
                    result[f"{label}_max"] = max(fva.loc[rxn.id, "maximum"], 0.0)
    return result


def load(name, path):
    if path is None:
        return cobra.io.load_model("textbook")   # built-in demo model
    return cobra.io.read_sbml_model(path)


def main():
    rows = []
    for name, path in STRAINS.items():
        model = load(name, path)

        # Warn about anything we could not find
        for label, ids in PRODUCTS.items():
            if first_present(model, ids) is None:
                print(f"[warn] {name}: no exchange found for {label}")
        if first_present(model, ATPM_IDS) is None:
            print(f"[warn] {name}: no ATP maintenance reaction found "
                  "(stress will only act through nutrient limits)")

        row = {"strain": name}
        for cname, cond in CONDITIONS.items():
            for key, val in simulate(model, cond).items():
                row[f"{key}_{cname}"] = val
        rows.append(row)

    df = pd.DataFrame(rows).set_index("strain")

    # Resilience = stressed growth / relaxed growth
    df["resilience"] = df["growth_stressed"] / df["growth_relaxed"]

    # Retention = stressed output / relaxed output, per product
    for label in PRODUCTS:
        s, r = f"{label}_max_stressed", f"{label}_max_relaxed"
        if s in df and r in df:
            df[f"{label}_retention"] = df[s] / df[r].replace(0, float("nan"))

    df = df.sort_values("resilience", ascending=False)
    pd.set_option("display.width", 200)
    print(df.round(3).T)
    df.to_csv("insilico_results.csv")
    print("\nSaved: insilico_results.csv")


if __name__ == "__main__":
    main()
