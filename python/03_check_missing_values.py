import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# 1. Load the ML-Ready Data
df_full = pd.read_csv("data/GSE57338_Raw_Multiclass.csv", index_col=0)

# 2. ISOLATE GENES (Drop 'Target' temporarily)
# We create a new dataframe 'df_genes' that contains ONLY the numeric gene values
df_genes = df_full.drop(columns=['Target'])

print(f"Analyzing Gene Data Only. Shape: {df_genes.shape} (Samples x Genes)")

# 3. Calculate Missing Values
total_missing = df_genes.isna().sum().sum()
print(f"\nTotal Missing Values in Gene Entries: {total_missing}")

# 4. Detailed Breakdown
# Check if any specific gene is empty
missing_per_gene = df_genes.isna().sum()
bad_genes = missing_per_gene[missing_per_gene > 0]

# Check if any specific sample is empty
missing_per_sample = df_genes.isna().sum(axis=1)
bad_samples = missing_per_sample[missing_per_sample > 0]

if len(bad_genes) > 0:
    print(f"\nFound {len(bad_genes)} genes with missing values.")
    print(bad_genes.head())
else:
    print("No genes have missing values.")

if len(bad_samples) > 0:
    print(f"\nFound {len(bad_samples)} samples with missing values.")
else:
    print("No samples have missing values.")
