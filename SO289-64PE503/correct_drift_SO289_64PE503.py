import pandas as pd
from scipy.interpolate import PchipInterpolator
from matplotlib import pyplot as plt
import matplotlib.dates as mdates

# Import dataframe
df = pd.read_csv("data/SO289-64PE503_results.csv")

# For now only keep Oct 31
L = df["dic_cell_id"] == "C_Oct31-22_0910"
df = df[L]

# Calculate nuts offset throughout the day
first_nuts = df.loc[df["bottle"]=="NUTSLAB01", "dic"].values
L= df["bottle"].str.startswith("NUTS")
df.loc[L, "nuts_offset"] = abs(df.loc[L, "dic"] - first_nuts)

# Use a PCHIP to interpolate the offset throughout the day
L = df["nuts_offset"].notnull()
interp_obj = PchipInterpolator(df.loc[L, 'analysis_datenum'], df.loc[L, 'nuts_offset'], extrapolate=False)
df['offset_pchip'] = interp_obj(df['analysis_datenum'])

# Correct DIC for samples
df["dic_corrected"] = df["dic"] + df["offset_pchip"]

# === PLOT
# Create a column with hours and minutes
df["analysis_datetime"] = pd.to_datetime(df["analysis_datetime"])
df["datetime"] = df["analysis_datetime"].dt.strftime("%H:%M")
df["datetime"] = pd.to_datetime(df["datetime"])

# Create figure
fig, ax = plt.subplots(dpi=300, figsize=(6, 4))

L = df["dic_corrected"].notnull()

# Scatter original DIC
ax.scatter(
    x="datetime",
    y="dic",
    data=df[L],
    alpha=0.3,
    label="Initial"
)    

# Scatter corrected DIC
ax.scatter(
    x="datetime",
    y="dic_corrected",
    data=df[L],
    alpha=0.3,
    label="Corrected"
)               

# Improve plot
myFmt = mdates.DateFormatter("%H")
ax.xaxis.set_major_formatter(myFmt)

ax.grid(alpha=0.3)
ax.set_xlabel("Time (hrs)")
ax.set_ylabel("$DIC$ / μmol · $kg^{-1}$")

ax.legend()

# Save plot
plt.tight_layout()
plt.savefig("./figs/correct_drift.png")
