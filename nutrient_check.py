print("script started")

import os
import cobra

BASE = os.path.dirname(os.path.abspath(__file__))
MODEL_FILE = os.path.join(BASE, "Models", "Lactobacillus_brevis_ATCC_367_AGORA2_v2.01.xml")

model = cobra.io.read_sbml_model(MODEL_FILE)
print("model loaded")

CAP = 5.0   # same cap used in the main script

with model:
    for r in model.exchanges:
        if r.lower_bound < -CAP:
            r.lower_bound = -CAP

    sol = model.optimize()
    print("Growth:", sol.objective_value)

    flux = sol.fluxes[[r.id for r in model.exchanges]]
    uptakes = flux[flux < -1e-6].sort_values()

    print("\nNutrients being taken up (negative = uptake):")
    print(uptakes.head(40).to_string())

    print("\nNutrients at the cap (these are limiting growth):")
    print(uptakes[uptakes <= -CAP + 1e-6].to_string())