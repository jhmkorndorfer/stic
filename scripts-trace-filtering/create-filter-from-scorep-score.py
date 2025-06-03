import pandas as pd

data = []

with open("new-scorep-score-no-header") as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) < 7:
            continue  # Skip malformed lines
        # First 6 fields are fixed: Type, Samples, Samples_wo_children, Time, Time%, CPU%
        fixed = parts[:6]
        function_name = " ".join(parts[6:])  # Join the rest as function name
        data.append(fixed + [function_name])

# Create DataFrame
df = pd.DataFrame(data, columns=[
    "Type", "buff", "visits",
    "Time", "Time_percent", "time_visit", "Function"
])

# Clean numeric columns
df["buff"] = df["buff"].str.replace(",", "").astype(float)
df["visits"] = df["visits"].str.replace(",", "").astype(float)

for col in ["Time", "Time_percent", "time_visit"]:
    df[col] = pd.to_numeric(df[col], errors="coerce")

# print the whole DataFrame and force all columns and rows to be displayed
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', None)

#Filter df to contail only rows with Time_percent > 1
df = df[df["Time_percent"] > 1]

# print(df.head(200))

for func in df["Function"]:
    print(f"{func:<}")
