import pandas as pd

# Import dataframe
df = pd.read_csv("data/SO289-64PE503_results.csv")

# Define the dates when the NUTS jug was refilled
refill_dates = ["2021-10-12", "2021-10-31", "2021-11-06", "2021-11-15", "2021-11-20"]
refill_dates = pd.to_datetime(refill_dates)

# Convert 'analysis_datetime' to datetime if it's not already
df['analysis_datetime'] = pd.to_datetime(df['analysis_datetime'])

# Filter for NUTS samples with flag = 2
df_nuts = df[(df['flag'] == 2) & (df['bottle'].str.startswith("NUTS"))]

# Exclude specific samples on specified dates
df_nuts = df_nuts[~((df_nuts['bottle'] == 'NUTSLAB01') & (df_nuts['station'] == 16))]
df_nuts = df_nuts[~((df_nuts['bottle'] == 'NUTSLAB02') & (df_nuts['station'] == 12))]

# Function to determine the jug refill date for each measurement
def assign_refill_date(row, refill_dates):
    for date in refill_dates:
        if row <= date:
            return date
    return refill_dates[-1]

# Assign a jug refill date to each measurement
df_nuts['jug_refill_date'] = df_nuts['analysis_datetime'].apply(lambda x: assign_refill_date(x, refill_dates))

# Group by jug refill date and acid batch
grouped = df_nuts.groupby(['jug_refill_date', 'analysis_batch'])

# Calculate the standard deviation for TA and DIC within each group
std_devs = grouped['alkalinity', 'dic'].std().reset_index()

# Calculate the overall uncertainty as the mean of the standard deviations
overall_ta_uncertainty = std_devs['alkalinity'].mean()
overall_dic_uncertainty = std_devs['dic'].mean()

# Print the overall uncertainties
print(f"Overall Uncertainty for TA: {overall_ta_uncertainty:.3f}")
print(f"Overall Uncertainty for DIC: {overall_dic_uncertainty:.3f}")
