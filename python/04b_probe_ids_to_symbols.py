import pandas as pd

# --- 1. LOAD YOUR EXISTING DATA ---
print("Loading existing Top 5000 dataset...")
# Update this if your file is named differently
df = pd.read_csv("data/GSE57338_Top5000_Var.csv", index_col=0)

y = df['Target']
X = df.drop(columns=['Target'])

# --- 2. LOAD THE R MAPPING DICTIONARY ---
print("Loading mapping dictionary...")
map_df = pd.read_csv("data/Probe_to_Symbol_Map.csv")

# Create a dictionary: {'1234567': 'GAPDH', '891011': 'TP53', ...}
mapping_dict = dict(zip(map_df['PROBEID'].astype(str), map_df['SYMBOL']))

# --- 3. FILTER AND RENAME COLUMNS ---
# Only keep the numerical columns that actually have a known gene mapping
valid_cols = [col for col in X.columns if col in mapping_dict]
X = X[valid_cols]

# Rename the numerical columns to Gene Symbols
X.rename(columns=mapping_dict, inplace=True)

# --- 4. AVERAGE ANY DUPLICATE GENES ---
# Sometimes two numeric probes map to the exact same gene. 
# We average them together so we only have one column per gene.
print("Averaging duplicate gene symbols...")
X = X.T.groupby(X.columns).mean().T

# --- 5. SAVE THE FINAL FILE ---
final_df = pd.concat([y, X], axis=1)
output_filename = "data/GSE57338_Top5000_Mapped.csv"
final_df.to_csv(output_filename)

print(f"\nSUCCESS! Your dataset now has {X.shape[1]} unique Gene Symbols.")
print(f"Saved as '{output_filename}'")