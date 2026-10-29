import os
BASE = os.path.dirname(os.path.abspath(__file__))
MODEL_FILE = os.path.join(BASE, "Models", "Lactobacillus_brevis_subsp_gravesensis_ATCC_27305_AGORA2_v2.01.xml ")
print("Looking for:", MODEL_FILE)
print("Exists?", os.path.exists(MODEL_FILE))
print("Folder contents:", os.listdir(os.path.join(BASE, "Models")))
import cobra
from cobra.flux_analysis import flux_variability_analysis

model = cobra.io.read_sbml_model(MODEL_FILE)
print("Loaded:", model.id, "|", len(model.reactions), "reactions")
print("Growth with default settings:", model.slim_optimize())

# 1. Look for glutamate decarboxylase (GAD), the GABA-making enzyme
print("\nPossible GAD reactions:")
found_gad = False
for r in model.reactions:
    text = (r.id + " " + r.name).lower()
    if "gludc" in text or "glutamate decarboxylase" in text:
        print("  ", r.id, "|", r.name, "|", r.reaction)
        found_gad = True
if not found_gad:
    print("   none found")

# 2. Look for a GABA exchange reaction (how GABA leaves the cell)
gaba_ex = [r for r in model.exchanges
           if "4abut" in r.id.lower() or "gaba" in r.id.lower()]
print("\nGABA exchange reactions:", [r.id for r in gaba_ex])

# 3. Can it make GABA while still growing?
if gaba_ex:
    ex = gaba_ex[0]
    with model:
        for g in model.exchanges:
            if g.id.lower().startswith(("ex_glu_l", "ex_glu__l")):
                g.lower_bound = -10          # open glutamate uptake
        growth = model.slim_optimize()
        print("\nGrowth with glutamate open:", growth)
        if growth == growth and growth > 1e-6:   # not NaN and above zero
            fva = flux_variability_analysis(
                model, reaction_list=[ex], fraction_of_optimum=0.9)
            print("Max GABA output at >=90% growth:", fva.loc[ex.id, "maximum"])
        else:
            print("Model did not grow, so GABA cannot be tested.")
else:
    print("\nNo GABA exchange, so this model cannot predict GABA output.")