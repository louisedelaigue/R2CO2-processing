import pandas as pd, numpy as np, koolstof as ks, calkulate as calk
import matplotlib.dates as mdates

# Import logfile and dbs file
logfile = ks.read_logfile(
    "data/TA_ONLY/logfile.bak",
    methods=[
        "3C standard AT only",
    ],
)

dbs = ks.read_dbs("data/TA_ONLY/SO289_TA_only.dbs", logfile=logfile)

# Convert datetime to datenum
dbs["analysis_datenum"] = mdates.date2num(dbs["analysis_datetime"])

# Create empty metadata columns
for meta in [
    "salinity",
    "alkalinity_certified",
    "total_phosphate",
    "total_silicate",
    "total_ammonium",
]:
    dbs[meta] = np.nan

# Assign metadata values for CRMs
dbs["crm"] = dbs.bottle.str.startswith("CRM")
crm_batches = [189, 195, 198]

for b in crm_batches:
    if b == 189:
        L = dbs["bottle"].str.startswith("CRM-{}".format(b))
        dbs.loc[L, "alkalinity_certified"] = 2205.26  # micromol/kg-sw
        dbs.loc[L, "salinity"] = 33.494
        dbs.loc[L, "total_phosphate"] = 0.45  # micromol/kg-sw
        dbs.loc[L, "total_silicate"] = 2.1  # micromol/kg-sw
        dbs.loc[L, "total_ammonium"] = 0  # micromol/kg-sw

    if b == 195:
        L = dbs["bottle"].str.startswith("CRM-{}".format(b))
        dbs.loc[L, "alkalinity_certified"] = 2213.51  # micromol/kg-sw
        dbs.loc[L, "salinity"] = 33.485
        dbs.loc[L, "total_phosphate"] = 0.49  # micromol/kg-sw
        dbs.loc[L, "total_silicate"] = 3.6  # micromol/kg-sw
        dbs.loc[L, "total_ammonium"] = 0  # micromol/kg-sw

    if b == 198:
        L = dbs["bottle"].str.startswith("CRM-{}".format(b))
        dbs.loc[L, "alkalinity_certified"] = 2200.67  # micromol/kg-sw
        dbs.loc[L, "salinity"] = 33.504
        dbs.loc[L, "total_phosphate"] = 0.67  # micromol/kg-sw
        dbs.loc[L, "total_silicate"] = 3.8  # micromol/kg-sw
        dbs.loc[L, "total_ammonium"] = 0  # micromol/kg-sw

# Assign temperature = 25.0 for VINDTA analysis temperature
dbs["temperature_override"] = 25.0

# Assign metadata for junks
dbs["salinity"] = dbs["salinity"].fillna(35)
dbs["total_phosphate"] = dbs["total_phosphate"].fillna(0)
dbs["total_silicate"] = dbs["total_silicate"].fillna(0)
dbs["total_ammonium"] = dbs["total_ammonium"].fillna(0)

# Assign alkalinity metadata
dbs["analyte_volume"] = 98.865  # TA pipette volume in ml
dbs["file_path"] = "data/TA_ONLY/SO289_TA_only/"

# Assign TA acid batches
dbs["analysis_batch"] = 0
dbs.loc[
    (dbs["analysis_datetime"].dt.day >= 24) & (dbs["analysis_datetime"].dt.month == 11),
    "analysis_batch",
] = 1
dbs.loc[
    (dbs["analysis_datetime"].dt.day >= 1) & (dbs["analysis_datetime"].dt.month == 12),
    "analysis_batch",
] = 1
dbs.loc[
    (dbs["analysis_datetime"].dt.day >= 3) & (dbs["analysis_datetime"].dt.month == 12),
    "analysis_batch",
] = 2

# Select which TA CRMs to use/avoid for calibration
dbs["reference_good"] = ~np.isnan(dbs.alkalinity_certified)

# Calibrate and solve alkalinity and plot calibration
calk.io.get_VINDTA_filenames(dbs)
calk.dataset.calibrate(dbs)
calk.dataset.solve(dbs)
calk.plot.titrant_molinity(
    dbs, figure_fname="figs/SO289-TA_only/titrant_molinity.png", show_bad=False
)
calk.plot.alkalinity_offset(
    dbs, figure_fname="figs/SO289-TA_only/alkalinity_offset.png", show_bad=False
)

# Demote dbs to a standard DataFrame
dbs = pd.DataFrame(dbs)

# Save to .csv
dbs.to_csv("data/SO289-TA_only_results.csv")