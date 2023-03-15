import pandas as pd
from scipy.interpolate import PchipInterpolator
from matplotlib import pyplot as plt
import matplotlib.dates as mdates

# Import dataframe
df = pd.read_csv("data/SO289-64PE503_results.csv")

# Only keep flag = 2
L = df["flag"] == 2
df = df[L]

# Only keep NUTS and samples
L = df["bottle"].str.startswith(("CRM-", "JUNK"))
df = df[~L]

# Create a list of analysis days
analysis_days = list(df.loc[df["real_day"] == True, "dic_cell_id"].unique())

# Create df to store corrected values
corrected_dic = pd.DataFrame()

# Calculate DIC offset for each analysis day
for d in analysis_days:

    L = df["dic_cell_id"] == d
    data = df[L].copy()

    # Calculate nuts offset throughout the day
    first_nuts = data.loc[data.loc[data.bottle.str.contains("NUTS")].index[0], "dic"]
    L = data["bottle"].str.startswith("NUTS")
    data.loc[L, "nuts_offset"] = abs(data.loc[L, "dic"] - first_nuts)

    # Use a PCHIP to interpolate the offset throughout the day
    L = data["nuts_offset"].notnull()
    interp_obj = PchipInterpolator(
        data.loc[L, "analysis_datenum"], data.loc[L, "nuts_offset"], extrapolate=False
    )
    data["offset_pchip"] = interp_obj(data["analysis_datenum"])

    # Correct DIC for samples
    data["dic_corrected"] = data["dic"] + data["offset_pchip"]

    # Store results
    corrected_dic = pd.concat([corrected_dic, data])

# Add corrected DIC to original df
samples = list(df["bottle"])

for s in samples:
    df.loc[df["bottle"] == s, "nuts_offset"] = corrected_dic.loc[
        corrected_dic["bottle"] == s, "nuts_offset"
    ]

    df.loc[df["bottle"] == s, "offset_pchip"] = corrected_dic.loc[
        corrected_dic["bottle"] == s, "offset_pchip"
    ]

    df.loc[df["bottle"] == s, "dic_corrected"] = corrected_dic.loc[
        corrected_dic["bottle"] == s, "dic_corrected"
    ]
    
# Save to .csv
df.to_csv("data/SO289-64PE503_results_recalibrated.csv", index=False)

# === PLOT
# Create a column with hours and minutes
df["analysis_datetime"] = pd.to_datetime(df["analysis_datetime"])
df["datetime_for_plotting_only"] = df["analysis_datetime"].dt.strftime("%H:%M")
df["datetime_for_plotting_only"] = pd.to_datetime(df["datetime_for_plotting_only"])

for d in analysis_days:

    # Create figure
    fig, ax = plt.subplots(dpi=300, figsize=(6, 4))

    L = df["dic_cell_id"] == d

    # Scatter original DIC
    ax.scatter(
        x="datetime_for_plotting_only", y="dic", data=df[L], alpha=0.3, label="Initial"
    )

    # Scatter corrected DIC
    ax.scatter(
        x="datetime_for_plotting_only",
        y="dic_corrected",
        data=df[L],
        alpha=0.3,
        label="Corrected",
    )

    # Improve plot
    myFmt = mdates.DateFormatter("%H")
    ax.xaxis.set_major_formatter(myFmt)

    ax.grid(alpha=0.3)
    ax.set_xlabel("Time / hrs")
    ax.set_ylabel("$DIC$ / μmol · $kg^{-1}$")

    ax.set_title(d)

    ax.legend()

    # Save plot
    plt.tight_layout()
    plt.savefig("./figs/drift_correction/correct_DIC_drift_{}.png".format(d))

# ==== PLOT ALL SAMPLES
# Create a column with hours and minutes
df["analysis_datetime"] = pd.to_datetime(df["analysis_datetime"])
df["datetime_for_plotting_only"] = df["analysis_datetime"].dt.strftime("%H:%M")
df["datetime_for_plotting_only"] = pd.to_datetime(df["datetime_for_plotting_only"])

# Create figure
fig, ax = plt.subplots(dpi=300, figsize=(6, 4))

L = df["dic_corrected"].notnull()

# Scatter original DIC
ax.scatter(
    x="datetime_for_plotting_only", y="dic", data=df[L], alpha=0.3, label="Initial"
)

# Scatter corrected DIC
ax.scatter(
    x="datetime_for_plotting_only",
    y="dic_corrected",
    data=df[L],
    alpha=0.3,
    label="Corrected",
)

# Improve plot
myFmt = mdates.DateFormatter("%H")
ax.xaxis.set_major_formatter(myFmt)

ax.grid(alpha=0.3)
ax.set_xlabel("Time / hrs")
ax.set_ylabel("$DIC$ / μmol · $kg^{-1}$")

ax.legend()

# Save plot
plt.tight_layout()
plt.savefig("./figs/drift_correction/correct_DIC_drift_all.png")

# Save nuts as csv
# nuts.to_csv("data/nuts.csv", index=False)
