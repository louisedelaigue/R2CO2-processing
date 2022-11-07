import pandas as pd
import PyCO2SYS as pyco2

# Import nuts data
df = pd.read_csv("data/nuts.csv")

# Only keep last three days of nuts
L = df.station <= 3
df = df[L]

# Only keep first nuts of the day
L = df["Unnamed: 0"].isin([423, 453, 484])
df = df[L]

# Average DIC and TA
DIC = df["dic"].mean()
TA = df["alkalinity"].mean()

# Calculate average pH
pH = pyco2.sys(
    TA,
    DIC,
    1,
    2,
    salinity=35,
    temperature=26,
    temperature_out=26,
    pressure=0,
    opt_pH_scale=1,
    total_phosphate=0,
    total_silicate=0,
)["pH_total"]

print(round(pH, 2))