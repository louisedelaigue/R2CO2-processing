import copy
import numpy as np, pandas as pd
from matplotlib import pyplot as plt
import itertools
import koolstof as ks, calkulate as calk

# Import logfile and dbs file
logfile = ks.read_logfile(
    "data/logfile.bak",
    methods=[
        "3C standard",
    ],
)

dbs = ks.read_dbs("data/64PE503_SO289_2022.dbs", logfile=logfile)

# Create empty metadata columns
for meta in [
    "salinity",
    "dic_certified",
    "alkalinity_certified",
    "total_phosphate",
    "total_silicate",
    "total_ammonium",
]:
    dbs[meta] = np.nan

# Assign metadata values for CRMs
dbs["crm"] = dbs.bottle.str.startswith("CRM")
crm_batches = [189, 195]

for b in crm_batches:
    if b == 189:
        L = dbs["bottle"].str.startswith("CRM-{}".format(b))
        dbs.loc[L, "dic_certified"] = 2009.48  # micromol/kg-sw
        dbs.loc[L, "alkalinity_certified"] = 2205.26  # micromol/kg-sw
        dbs.loc[L, "salinity"] = 33.494
        dbs.loc[L, "total_phosphate"] = 0.45  # micromol/kg-sw
        dbs.loc[L, "total_silicate"] = 2.1  # micromol/kg-sw
        dbs.loc[L, "total_ammonium"] = 0  # micromol/kg-sw
    
    if b == 195:
        L = dbs["bottle"].str.startswith("CRM-{}".format(b))
        dbs.loc[L, "dic_certified"] = 2024.96  # micromol/kg-sw
        dbs.loc[L, "alkalinity_certified"] = 2213.51  # micromol/kg-sw
        dbs.loc[L, "salinity"] = 33.485
        dbs.loc[L, "total_phosphate"] = 0.49  # micromol/kg-sw
        dbs.loc[L, "total_silicate"] = 3.6  # micromol/kg-sw
        dbs.loc[L, "total_ammonium"] = 0  # micromol/kg-sw
    
# Assign temperature = 25.0 for VINDTA analysis temperature
dbs["temperature_override"] = 25.0

# Assign metadata for junks
dbs['salinity'] = dbs['salinity'].fillna(35)
dbs['total_phosphate'] = dbs['total_phosphate'].fillna(0)
dbs['total_silicate'] = dbs['total_silicate'].fillna(0)
dbs['total_ammonium'] = dbs['total_ammonium'].fillna(0)

# Add optional column "file_good"
dbs['file_good'] = True

# === ALKALINITY
# Assign alkalinity metadata
dbs["analyte_volume"] = 95.939  # TA pipette volume in ml
dbs["file_path"] = "data/64PE503_SO289_2022/"

# Assign TA acid batches
dbs["analysis_batch"] = 0
dbs.loc[dbs['analysis_datetime'].dt.month == 10, 'analysis_batch'] = 1

# Select which TA CRMs to use/avoid for calibration
dbs["reference_good"] = ~np.isnan(dbs.alkalinity_certified)

# Calibrate and solve alkalinity and plot calibration
calk.io.get_VINDTA_filenames(dbs)
calk.dataset.calibrate(dbs)
calk.dataset.solve(dbs)
calk.plot.titrant_molinity(dbs, figure_fname="figs/titrant_molinity.png", show_bad=False)
calk.plot.alkalinity_offset(dbs, figure_fname="figs/alkalinity_offset.png", show_bad=False)

# === DIC
# Add optional column "blank_good
dbs['blank_good'] = True

# Select which DIC CRMs to use/avoid for calibration --- only fresh bottles
dbs["k_dic_good"] = dbs.crm & dbs.bottle.str.endswith("-01")

# Get blanks and apply correction
dbs.get_blank_corrections()
dbs.plot_blanks(figure_path="figs/dic_blanks/")

# Calibrate DIC and plot calibration
dbs.calibrate_dic()
dic_sessions = copy.deepcopy(dbs.sessions)
dbs.plot_k_dic(figure_path="figs/")
dbs.plot_dic_offset(figure_path="figs/")

# === ENTIRE DATASET
# Demote dbs to a standard DataFrame
dbs = pd.DataFrame(dbs)

# === PLOT NUTS FOR EACH ANALYSIS DAY
# Prepare colours and markers
markers = itertools.cycle(("o", "^", "s", "v", "D", "<", ">"))
colors = itertools.cycle(
    (
        "xkcd:purple",
        "xkcd:green",
        "xkcd:blue",
        "xkcd:pink",
        "xkcd:deep blue",
        "xkcd:red",
        "xkcd:teal",
        "xkcd:orange",
        "xkcd:fuchsia",
    )
)

# Only keep nuts bottles
L = dbs["bottle"].str.startswith("NUTS")
nuts = dbs[L]

# Only keep real analysis days
real_days = ["C_Aug27-22_0908", "C_Oct17-22_0910", "C_Oct18-22_0810"]
L = nuts["dic_cell_id"].isin(real_days)
nuts = nuts[L]

# Create an hour column
nuts["hour"] = nuts["analysis_datetime"].dt.hour

# Create figure
fig, ax = plt.subplots(dpi=300, figsize=(6, 4))

# Scatter NUTS DIC
for r in real_days:
    L = nuts["dic_cell_id"] == r
    data = nuts[L]
    m = next(markers)
    c = next(colors)
    ax.scatter(
        x="hour",
        y="dic",
        data=data,
        marker=m,
        color=c,
        alpha=0.3,
        label=r,
    )
    
ax.legend()
ax.grid(alpha=0.3)
ax.set_xlabel("Time (hrs)")
ax.set_ylabel("$DIC$ / μmol · $kg^{-1}$")

# Save plot
plt.tight_layout()
plt.savefig("./figs/day_nuts.png")

