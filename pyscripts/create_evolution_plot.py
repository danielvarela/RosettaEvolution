

import sys
import pandas as pd
import matplotlib.pyplot as plt


input_file = sys.argv[-1]

with open(input_file, "r") as f:
    lines = [d.strip() for d in  f.readlines()]

rmsds = [l for l in lines if "[RMSDS]" in l]
energies = [l for l in lines if "[POPUL]" in l]
gens = [l for l in lines if "[GEN]" in l]

results = []
for g in gens:
    tokens = g.split(" ")
    print(tokens)
    best_energy = float(tokens[2])
    avg_energy = float(tokens[3])
    results.append({"bst" : best_energy, "avg" : avg_energy })

df = pd.DataFrame(results)
df["gen"] = range(0,len(df))
ax = df.plot.line(x="gen")

ax.get_figure().savefig("evolution.png")

plt.clf()
plt.close()
