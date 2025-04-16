import pandas as pd, numpy as np, koolstof as ks, calkulate as calk
import matplotlib.dates as mdates

# Import logfile and dbs file765t4rfgedvc 56tyrghujugytferdtyujikhhgfdcvxb nmbvc
logfile = ks.read_logfile(
    "data/TA_ONLY/logfile.bak",
    methods=[
        "3C standard AT only",
    ],
)

dbs = ks.read_dbs("data/TA_ONLY/SO289_TA_only.dbs")

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
        
# Assign metadata values for SO289 bottle samples
so289_metadata = pd.read_csv("data/SO289_CTD_data.csv")
so289_metadata["bottle"] = ["SO289-" + str(s) for s in so289_metadata["bottle"]]
so289_samples = list(so289_metadata["bottle"])

for s in so289_samples:
    dbs.loc[dbs["bottle"] == s, "salinity"] = so289_metadata.loc[
        so289_metadata["bottle"] == s, "salinity"
    ].values
    dbs.loc[dbs["bottle"] == s, "total_phosphate"] = so289_metadata.loc[
        so289_metadata["bottle"] == s, "phosphate"
    ].values
    dbs.loc[dbs["bottle"] == s, "total_silicate"] = so289_metadata.loc[
        so289_metadata["bottle"] == s, "silicate"
    ].values

# Assign metadata for SO289/ UWS samples
so289_metadata_uws = pd.read_csv("data/TA_ONLY/uws_S06_match_samples_SMB_sal_temp.csv")
so289_samples_uws = list(so289_metadata_uws["bottle"])

for s in so289_samples_uws:
    dbs.loc[dbs["bottle"] == s, "salinity"] = so289_metadata_uws.loc[
        so289_metadata_uws["bottle"] == s, "SBE45_sal"
    ].values
    

print("/!\ \n CAREFUL: MISSING METADATA FOR 3 SO289 UWS SAMPLES \n /!\ ")

# Assign metadata for alkalinity experiement samples
meta = pd.read_excel("data/TA_ONLY/SO289_NaOH.xlsx", na_values="<LOD")

# Calculate sample density
meta["density"] = calk.density.seawater_1atm_MP81(25, meta["Salinity"])

# Convert nutrients from nM to umol/kg
meta["ton"] = (meta["TON (nM)"] * 0.001) / meta["density"] # from nM to umol/L, then divide by density to umol/kg
meta["total_phosphate"] = (meta["Phosphate (nM)"] * 0.001) / meta["density"] # from nM to umol/L, then divide by density to umol/kg
meta["total_silicate"] = (meta["Silicate (nM)"] * 0.001) / meta["density"] # from nM to umol/L, then divide by density to umol/kg

# Isolate alkalinity experiment
L = dbs.bottle.str.startswith("AE")

# Create a list of samples
ae_samples = list(dbs[L].bottle.unique())

# Assign experiment number based on sample name
for s in ae_samples:
    if len(s.split("-")) == 3:
        experiment = s.split("-")[1]
        dbs.loc[dbs["bottle"]==s, "exp_number"] = int(experiment[3])
        
    if len(s.split("-")) == 4:
        experiment = s.split("-")[2]
        dbs.loc[dbs["bottle"]==s, "exp_number"] = int(experiment[1])

# Drop experiment 8 because of contamination
L = dbs["exp_number"] == 8
dbs = dbs[~L]

# Create a list of experiments
experiments = list(meta["Exp."].unique())

for e in experiments:  
    dbs.loc[dbs["exp_number"]==e, "sample_date"] = meta.loc[meta["Exp."]==e, "Date"].item()
    dbs.loc[dbs["exp_number"]==e, "time"] = meta.loc[meta["Exp."]==e, "Time"].item()
    dbs.loc[dbs["exp_number"]==e, "latitude"] = meta.loc[meta["Exp."]==e, "Lat"].item()
    dbs.loc[dbs["exp_number"]==e, "longitude"] = meta.loc[meta["Exp."]==e, "Lon"].item()
    dbs.loc[dbs["exp_number"]==e, "temperature"] = meta.loc[meta["Exp."]==e, "SST"].item()
    dbs.loc[dbs["exp_number"]==e, "salinity"] = meta.loc[meta["Exp."]==e, "Salinity"].item()
    dbs.loc[dbs["exp_number"]==e, "ton"] = meta.loc[meta["Exp."]==e, "ton"].item()
    dbs.loc[dbs["exp_number"]==e, "total_phosphate"] = meta.loc[meta["Exp."]==e, "total_phosphate"].item()
    dbs.loc[dbs["exp_number"]==e, "total_silicate"] = meta.loc[meta["Exp."]==e, "total_silicate"].item()

# Assign temperature = 25.0 for VINDTA analysis temperature
dbs["temperature_override"] = 25.0

# Assign metadata for junks, nutswater and storage test samples
L = dbs.bottle.str.startswith(("JUNK", "NUTS", "ST"))
dbs.loc[L, "salinity"] = dbs["salinity"].fillna(35)
dbs.loc[L, "total_phosphate"] = dbs["total_phosphate"].fillna(0)
dbs.loc[L, "total_silicate"] = dbs["total_silicate"].fillna(0)
dbs.loc[L, "total_ammonium"] = dbs["total_ammonium"].fillna(0)

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