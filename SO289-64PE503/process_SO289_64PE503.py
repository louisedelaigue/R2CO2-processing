import copy
import numpy as np, pandas as pd
from pandas.tseries.offsets import DateOffset
import matplotlib.dates as dates
from matplotlib import pyplot as plt
import matplotlib.dates as mdates
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

# Fix datetime issue in .dbs samples
dbs.loc[dbs["bottle"] == "SO289-41297", "analysis_datetime"] = "2022-10-21 11:58:00"
dbs.loc[dbs["bottle"] == "64PE503-28-5-8", "analysis_datetime"] = "2022-10-21 12:12:00"
dbs.loc[dbs["bottle"] == "64PE503-68-6-5", "analysis_datetime"] = "2022-10-21 12:44:00"
dbs.loc[dbs["bottle"] == "SO289-40698", "analysis_datetime"] = "2022-10-21 13:00:00"
dbs.loc[dbs["bottle"] == "64PE503-46-4-8", "analysis_datetime"] = "2022-10-21 13:30:00"
dbs.loc[dbs["bottle"] == "64PE503-56-4-2", "analysis_datetime"] = "2022-10-21 13:46:00"
dbs.loc[dbs["bottle"] == "SO289-40611", "analysis_datetime"] = "2022-10-21 14:18:00"
dbs.loc[dbs["bottle"] == "SO289-40499", "analysis_datetime"] = "2022-10-21 14:33:00"
dbs.loc[dbs["bottle"] == "64PE503-58-4-5", "analysis_datetime"] = "2022-10-21 14:50:00"
dbs.loc[dbs["bottle"] == "64PE503-66-5-8", "analysis_datetime"] = "2022-10-21 15:09:00"
dbs.loc[dbs["bottle"] == "64PE503-14-4-11", "analysis_datetime"] = "2022-10-21 15:24:00"
dbs.loc[dbs["bottle"] == "64PE503-47-4-11", "analysis_datetime"] = "2022-10-21 15:39:00"
dbs.loc[dbs["bottle"] == "64PE503-26-5-8", "analysis_datetime"] = "2022-10-21 15:55:00"
dbs.loc[dbs["bottle"] == "CRM-189-0526-01", "analysis_datetime"] = "2022-10-21 16:28:00"
dbs.loc[dbs["bottle"] == "CRM-189-0526-02", "analysis_datetime"] = "2022-10-21 17:17:00"

# Fix datetime issue in .dbs for sample "64PE503-12-7-2" // because there's a duplicate
L = (dbs["dic_cell_id"] == "C_Oct21-22_0810") & (dbs["bottle"] == "64PE503-12-7-2")
dbs.loc[L, "analysis_datetime"] = "2022-10-21 13:14:00"

# Fix datetime issue in .dbs nuts
L = (dbs["bottle"] == "NUTSLAB03") & (dbs["analysis_datetime"].isnull())
dbs.loc[L, "analysis_datetime"] = "2022-10-21 12:28:00"
L = (dbs["bottle"] == "NUTSLAB04") & (dbs["analysis_datetime"].isnull())
dbs.loc[L, "analysis_datetime"] = "2022-10-21 14:02:00"
L = (dbs["bottle"] == "NUTSLAB05") & (dbs["analysis_datetime"].isnull())
dbs.loc[L, "analysis_datetime"] = "2022-10-21 16:11:00"

# Convert datetime to datenum
dbs["analysis_datenum"] = dates.date2num(dbs["analysis_datetime"])

# Add a column to locate real analysis days
dbs["real_day"] = True
dbs.loc[dbs["dic_cell_id"] == "C_Aug26-22_0808", "real_day"] = False
dbs.loc[dbs["dic_cell_id"] == "C_Aug28-22_0808", "real_day"] = False
dbs.loc[dbs["dic_cell_id"] == "C_Aug30-22_0808", "real_day"] = False

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
dbs["salinity"] = dbs["salinity"].fillna(35)
dbs["total_phosphate"] = dbs["total_phosphate"].fillna(0)
dbs["total_silicate"] = dbs["total_silicate"].fillna(0)
dbs["total_ammonium"] = dbs["total_ammonium"].fillna(0)

# Add optional column "file_good"
dbs["file_good"] = True

L = (
    (dbs["bottle"] == "NUTSLAB05")
    & (dbs["analysis_datetime"].dt.month == 10)
    & (dbs["analysis_datetime"].dt.day == 28)
)
dbs.loc[L, "file_good"] = False

L = (
    (dbs["bottle"] == "NUTSLAB03")
    & (dbs["analysis_datetime"].dt.month == 10)
    & (dbs["analysis_datetime"].dt.day == 25)
)
dbs.loc[L, "file_good"] = False

# Add a flag column
# where good = 4, questionable = 3, bad = 2
dbs["flag"] = 4

L = (
    (dbs["bottle"] == "NUTSLAB05")
    & (dbs["analysis_datetime"].dt.month == 10)
    & (dbs["analysis_datetime"].dt.day == 28)
)
dbs.loc[L, "flag"] = 2

L = (
    (dbs["bottle"] == "NUTSLAB03")
    & (dbs["analysis_datetime"].dt.month == 10)
    & (dbs["analysis_datetime"].dt.day == 25)
)
dbs.loc[L, "flag"] = 2

dbs.loc[dbs["bottle"] == "64PE503-10-5-2", "flag"] = 3  # weird DIC
dbs.loc[dbs["bottle"] == "SO289-41003", "flag"] = 3  # weird DIC
dbs.loc[dbs["bottle"] == "64PE503-53-9-3", "flag"] = 3  # red tape
dbs.loc[dbs["bottle"] == "64PE503-53-4-5", "flag"] = 3  # red tape
dbs.loc[dbs["bottle"] == "64PE503-53-4-6", "flag"] = 3  # red tape
dbs.loc[dbs["bottle"] == "64PE503-53-4-9", "flag"] = 3  # red tape
dbs.loc[dbs["bottle"] == "64PE503-57-3-2", "flag"] = 3  # red tape
dbs.loc[dbs["bottle"] == "64PE503-53-4-2", "flag"] = 3  # red tape
dbs.loc[dbs["bottle"] == "SO289-40060", "flag"] = 3  # bottle popped, DIC only 1 x rinse
dbs.loc[dbs["bottle"] == "SO289-40975", "flag"] = 3  # not so great TA curve during analysis
dbs.loc[dbs["bottle"] == "SO289-41296", "flag"] = 3  # not so great TA curve during analysis

# Flag any nan
if dbs["counts"].isnull().any():
    print("ERROR: DIC counts include nan values.")

# === ALKALINITY
# Assign alkalinity metadata
dbs["analyte_volume"] = 98.865  # TA pipette volume in ml
dbs["file_path"] = "data/64PE503_SO289_2022/"

# Assign TA acid batches
dbs["analysis_batch"] = 0
dbs.loc[
    (dbs["analysis_datetime"].dt.day >= 11) & (dbs["analysis_datetime"].dt.month == 10),
    "analysis_batch",
] = 1
dbs.loc[
    (dbs["analysis_datetime"].dt.day >= 20) & (dbs["analysis_datetime"].dt.month == 10),
    "analysis_batch",
] = 2
dbs.loc[
    (dbs["analysis_datetime"].dt.day >= 28) & (dbs["analysis_datetime"].dt.month == 10),
    "analysis_batch",
] = 3

# Select which TA CRMs to use/avoid for calibration
dbs["reference_good"] = ~np.isnan(dbs.alkalinity_certified)

# Calibrate and solve alkalinity and plot calibration
calk.io.get_VINDTA_filenames(dbs)
calk.dataset.calibrate(dbs)
calk.dataset.solve(dbs)
calk.plot.titrant_molinity(
    dbs, figure_fname="figs/titrant_molinity.png", show_bad=False
)
calk.plot.alkalinity_offset(
    dbs, figure_fname="figs/alkalinity_offset.png", show_bad=False
)

# === DIC
# Add optional column "blank_good
dbs["blank_good"] = True
dbs.loc[dbs["bottle"] == "CRM-189-0526-02", "blank_good"] = False
dbs.loc[dbs["bottle"] == "64PE503-10-5-2", "blank_good"] = False
dbs.loc[dbs["bottle"] == "SO289-41003", "blank_good"] = False
dbs.loc[dbs["bottle"] == "SO289-41003", "blank_good"] = False

dbs.loc[
    (dbs["bottle"] == "JUNK01") & (dbs["dic_cell_id"] == "C_Oct24-22_0810"),
    "blank_good",
] = False
dbs.loc[
    (dbs["bottle"] == "JUNK04") & (dbs["dic_cell_id"] == "C_Oct28-22_0910"),
    "blank_good",
] = False
dbs.loc[
    (dbs["bottle"] == "NUTSLAB05") & (dbs["dic_cell_id"] == "C_Oct28-22_0910"),
    "blank_good",
] = False

# Select which DIC CRMs to use/avoid for calibration --- only fresh bottles
dbs["k_dic_good"] = dbs.crm & dbs.bottle.str.startswith("CRM")
dbs.loc[
    dbs.crm & dbs.bottle.str.endswith("-02"), "k_dic_good"
] = False  # remove CRMs used twice for DIC calibration

# Get blanks and apply correction
dbs.get_blank_corrections()
dbs.plot_blanks(figure_path="figs/dic_blanks/")

# Calibrate DIC and calibration
dbs.calibrate_dic()
dic_sessions = copy.deepcopy(dbs.sessions)
dbs.plot_k_dic(figure_path="figs/")
dbs.plot_dic_offset(figure_path="figs/")

# === ENTIRE DATASET
# Fix datetime for 31 October 2022 (time change in real life)
L = (dbs["analysis_datetime"].dt.month == 10) & (dbs["analysis_datetime"].dt.day == 31)
dbs.loc[L, "analysis_datetime"] = dbs["analysis_datetime"] - DateOffset(hours=1)

# Correct bottle name typos
dbs.loc[dbs["bottle"] == "CRM-189-0897", "bottle"] = "CRM-189-0897-01"
dbs.loc[dbs["bottle"] == "64PE503-65-2-1", "bottle"] = "64PE503-65-1-2"

# Demote dbs to a standard DataFrame
dbs = pd.DataFrame(dbs)

# Save to .csv
dbs.to_csv("data/SO289-64PE503_results.csv")

# === PLOT NUTS FOR EACH "REAL" ANALYSIS DAY
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

# Only keep nuts bottles that have file_good
L = (dbs["bottle"].str.startswith("NUTS")) & (dbs["file_good"] == True)
nuts = dbs[L]

# Only keep real analysis days
L = nuts["real_day"] == True
nuts = nuts[L]

real_days = list(nuts["dic_cell_id"].unique())

# Create a column with hours and minutes
nuts["datetime"] = nuts["analysis_datetime"].dt.strftime("%H:%M")
nuts["datetime"] = pd.to_datetime(nuts["datetime"])

# Create figure
fig, ax = plt.subplots(dpi=300, figsize=(6, 4))

# Scatter NUTS DIC
for r in real_days:
    L = nuts["dic_cell_id"] == r
    data = nuts[L]
    m = next(markers)
    c = next(colors)
    l = r.split("_")[1].replace("-22", "")
    ax.scatter(
        x="datetime",
        y="dic",
        data=data,
        marker=m,
        color=c,
        alpha=0.3,
        label=l,
    )

    # Fit a polynomial
    a, b, d = np.polyfit(data["analysis_datenum"], data["dic"], 2)
    new_dic = a * (data["analysis_datenum"] ** 2) + b * (data["analysis_datenum"]) + d

    # Plot polynomial
    ax.plot(data["datetime"], new_dic, color=c, alpha=0.3)

myFmt = mdates.DateFormatter("%H")
ax.xaxis.set_major_formatter(myFmt)

ax.legend(loc="upper left", ncol=3)  # bbox_to_anchor=(1, 0.5)
# ax.set_ylim(2000, 2200)
ax.grid(alpha=0.3)
ax.set_xlabel("Time (hrs)")
ax.set_ylabel("$DIC$ / μmol · $kg^{-1}$")


# Save plot
plt.tight_layout()
plt.savefig("./figs/day_nuts.png")
