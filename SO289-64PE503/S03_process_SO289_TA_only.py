import numpy as np, pandas as pd, koolstof as ks, calkulate as calk
import copy, itertools
from pandas.tseries.offsets import DateOffset
import matplotlib.dates as mdates
from matplotlib import pyplot as plt

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
    "dic_certified",
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

    if b == 198:
        L = dbs["bottle"].str.startswith("CRM-{}".format(b))
        dbs.loc[L, "dic_certified"] = 2033.64  # micromol/kg-sw
        dbs.loc[L, "alkalinity_certified"] = 2200.67  # micromol/kg-sw
        dbs.loc[L, "salinity"] = 33.504
        dbs.loc[L, "total_phosphate"] = 0.67  # micromol/kg-sw
        dbs.loc[L, "total_silicate"] = 3.8  # micromol/kg-sw
        dbs.loc[L, "total_ammonium"] = 0  # micromol/kg-sw