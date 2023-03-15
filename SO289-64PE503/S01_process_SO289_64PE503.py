import numpy as np, pandas as pd, koolstof as ks, calkulate as calk
from koolstof import vindta as ksv
import copy, itertools
from pandas.tseries.offsets import DateOffset
import matplotlib.dates as mdates
from matplotlib import pyplot as plt

# Import logfile and dbs file
# logfile = ks.read_logfile(
#     "data/logfile.bak",
#     methods=[
#         "3C standard",
#     ],
# )


logfile = ksv.read_logfile("data/logfile.bak", methods="3C standard")

dbs = ks.read_dbs("data/64PE503_SO289_2022.dbs") # no need for logfile anymore?

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

# Fix datetime issue for duplicate samples
L = (dbs["dic_cell_id"] == "C_Oct21-22_0810") & (dbs["bottle"] == "64PE503-12-7-2")
dbs.loc[L, "analysis_datetime"] = "2022-10-21 13:14:00"

L = (dbs["dic_cell_id"] == "C_Nov01-22_0811") & (dbs["bottle"] == "64PE503-66-5-8")
dbs.loc[L, "analysis_datetime"] = "2022-11-01 16:01:00"

L = (dbs["dic_cell_id"] == "C_Nov01-22_0811") & (dbs["bottle"] == "64PE503-41-4-5")
dbs.loc[L, "analysis_datetime"] = "2022-11-01 11:09:00"

# Fix datetime issue in .dbs nuts
L = (dbs["bottle"] == "NUTSLAB03") & (dbs["analysis_datetime"].isnull())
dbs.loc[L, "analysis_datetime"] = "2022-10-21 12:28:00"
L = (dbs["bottle"] == "NUTSLAB04") & (dbs["analysis_datetime"].isnull())
dbs.loc[L, "analysis_datetime"] = "2022-10-21 14:02:00"
L = (dbs["bottle"] == "NUTSLAB05") & (dbs["analysis_datetime"].isnull())
dbs.loc[L, "analysis_datetime"] = "2022-10-21 16:11:00"

# Convert datetime to datenum
dbs["analysis_datenum"] = mdates.date2num(dbs["analysis_datetime"])

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

# Assign metadata values for SO289 cruise samples
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
# where good = 2, questionable = 3, bad = 4
dbs["flag"] = 2

L = (
    (dbs["bottle"] == "NUTSLAB05")
    & (dbs["analysis_datetime"].dt.month == 10)
    & (dbs["analysis_datetime"].dt.day == 28)
)
dbs.loc[L, "flag"] = 4

L = (
    (dbs["bottle"] == "NUTSLAB03")
    & (dbs["analysis_datetime"].dt.month == 10)
    & (dbs["analysis_datetime"].dt.day == 25)
)
dbs.loc[L, "flag"] = 4

# 64PE503 CRUISE
dbs.loc[dbs["bottle"] == "64PE503-10-5-2", "flag"] = 3  # weird DIC
dbs.loc[dbs["bottle"] == "64PE503-53-9-3", "flag"] = 3  # red tape
dbs.loc[dbs["bottle"] == "64PE503-53-4-5", "flag"] = 3  # red tape
dbs.loc[dbs["bottle"] == "64PE503-53-4-6", "flag"] = 3  # red tape
dbs.loc[dbs["bottle"] == "64PE503-53-4-9", "flag"] = 3  # red tape
dbs.loc[dbs["bottle"] == "64PE503-57-3-2", "flag"] = 3  # red tape
dbs.loc[dbs["bottle"] == "64PE503-53-4-2", "flag"] = 3  # red tape
dbs.loc[dbs["bottle"] == "64PE503-42-4-2", "flag"] = 3  # red tape
dbs.loc[dbs["bottle"] == "64PE503-66-5-8", "flag"] = 3  # red tape

# SO289 CRUISE
dbs.loc[dbs["bottle"] == "SO289-41003", "flag"] = 3  # weird DIC
dbs.loc[dbs["bottle"] == "SO289-40060", "flag"] = 3  # bottle popped, DIC only 1 x rinse
dbs.loc[
    dbs["bottle"] == "SO289-40975", "flag"
] = 3  # not so great TA curve during analysis
dbs.loc[
    dbs["bottle"] == "SO289-41296", "flag"
] = 3  # not so great TA curve during analysis
dbs.loc[
    dbs["bottle"] == "SO289-41304", "flag"
] = 3  # not so great TA curve during analysis
dbs.loc[dbs["bottle"] == "SO289-40326", "flag"] = 3  # electrical tape off
dbs.loc[
    dbs["bottle"] == "SO289-40760", "flag"
] = 3  # not so great TA curve during analysis
dbs.loc[
    dbs["bottle"] == "SO289-40832", "flag"
] = 3  # not so great TA curve during analysis
dbs.loc[
    dbs["bottle"] == "SO289-41047", "flag"
] = 3  # not so great TA curve during analysis
dbs.loc[dbs["bottle"] == "SO289-40887", "flag"] = 3  # electrical tape off
dbs.loc[
    dbs["bottle"] == "SO289-41045", "flag"
] = 3  # not so great TA curve during analysis

# Add lab book notes as column
dbs["notes"] = np.nan

dbs.loc[dbs["bottle"] == "64PE503-10-5-2", "notes"] = "weird DIC"
dbs.loc[dbs["bottle"] == "64PE503-53-9-3", "notes"] = "red tape"
dbs.loc[dbs["bottle"] == "64PE503-53-4-5", "notes"] = "red tape"
dbs.loc[dbs["bottle"] == "64PE503-53-4-6", "notes"] = "red tape"
dbs.loc[dbs["bottle"] == "64PE503-53-4-9", "notes"] = "red tape"
dbs.loc[dbs["bottle"] == "64PE503-57-3-2", "notes"] = "red tape"
dbs.loc[dbs["bottle"] == "64PE503-53-4-2", "notes"] = "red tape"
dbs.loc[dbs["bottle"] == "64PE503-42-4-2", "notes"] = "red tape"
dbs.loc[dbs["bottle"] == "64PE503-66-5-8", "notes"] = "red tape"

dbs.loc[dbs["bottle"] == "SO289-41003", "notes"] = "weird DIC"
dbs.loc[dbs["bottle"] == "SO289-40060", "notes"] = "bottle popped, DIC only 1 x rinse"
dbs.loc[
    dbs["bottle"] == "SO289-40975", "notes"
] = "not so great TA curve during analysis"
dbs.loc[
    dbs["bottle"] == "SO289-41296", "notes"
] = "not so great TA curve during analysis"
dbs.loc[
    dbs["bottle"] == "SO289-41304", "notes"
] = "not so great TA curve during analysis"
dbs.loc[dbs["bottle"] == "SO289-40326", "notes"] = "electrical tape off"
dbs.loc[
    dbs["bottle"] == "SO289-40760", "notes"
] = "not so great TA curve during analysis"
dbs.loc[
    dbs["bottle"] == "SO289-40832", "notes"
] = "not so great TA curve during analysis"
dbs.loc[
    dbs["bottle"] == "SO289-41047", "notes"
] = "not so great TA curve during analysis"
dbs.loc[dbs["bottle"] == "SO289-40887", "notes"] = "electrical tape off"
dbs.loc[
    dbs["bottle"] == "SO289-41045", "notes"
] = "not so great TA curve during analysis"

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
    (dbs["analysis_datetime"].dt.day >= 28) & (dbs["analysis_datetime"].dt.month >= 10),
    "analysis_batch",
] = 3

dbs.loc[
    (dbs["analysis_datetime"].dt.day >= 11) & (dbs["analysis_datetime"].dt.month >= 11),
    "analysis_batch",
] = 4

dbs.loc[
    (dbs["analysis_datetime"].dt.day >= 20) & (dbs["analysis_datetime"].dt.month >= 11),
    "analysis_batch",
] = 5

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
dbs["k_dic_good"] = ~dbs.dic_certified.isnull()
dbs.loc[
    dbs.crm & dbs.bottle.str.endswith("-02"), "k_dic_good"
] = False  # remove CRMs used twice for DIC calibration

# Assign analysis temperature, if it's not 25 °C
dbs["temperature_analysis_dic"] = 23.0

# Assign dbs["salinity"] here too, if it's not always 35

sessions = ksv.blank_correction(
    dbs,
    logfile,
    blank_col="blank",
    counts_col="counts",
    runtime_col="run_time",
    session_col="dic_cell_id",
    use_from=6,
)

ksv.plot_increments(dbs, logfile, use_from=6)
ksv.plot_blanks(dbs, sessions)

# Calibrate
ksv.calibrate_dic(dbs, sessions)

# Plot calibration factors
ksv.plot_k_dic(dbs, sessions, show_ignored=True)

# Plot CRM offsets
ksv.plot_dic_offset(dbs, sessions)

# === ENTIRE DATASET
# Fix datetime for 31 October 2022 (time change in real life)
L = (dbs["analysis_datetime"].dt.month == 10) & (dbs["analysis_datetime"].dt.day == 31)
dbs.loc[L, "analysis_datetime"] = dbs["analysis_datetime"] - DateOffset(hours=1)

# Correct bottle name typos
dbs.loc[dbs["bottle"] == "CRM-189-0897", "bottle"] = "CRM-189-0897-01"
dbs.loc[dbs["bottle"] == "64PE503-65-2-1", "bottle"] = "64PE503-65-1-2"
dbs.loc[dbs["bottle"] == "46PE503-20-4-5", "bottle"] = "64PE503-20-4-5"

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

ax.legend(loc="upper left", ncol=4)  # bbox_to_anchor=(1, 0.5)
# ax.set_ylim(2000, 2200)
ax.grid(alpha=0.3)
ax.set_xlabel("Time / hrs")
ax.set_ylabel("$DIC$ / μmol · $kg^{-1}$")

# Save plot
plt.tight_layout()
plt.savefig("./figs/day_nuts.png")
