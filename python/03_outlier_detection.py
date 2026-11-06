import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram, linkage
from sklearn.preprocessing import StandardScaler

# 1. Load Data
# Usage: python python/03_outlier_detection.py [GSE57338|GSE5406]   (default: GSE57338, the training cohort)
DATASET = sys.argv[1] if len(sys.argv) > 1 else "GSE57338"
df = pd.read_csv(f"data/{DATASET}_Raw_Multiclass.csv", index_col=0)
OUT = f"data/{DATASET}_No_Outliers.csv"

# Separate Target and Features
y = df['Target']
X = df.drop(columns=['Target'])

# 2. Standardize Data (Crucial for distance calculations)
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

print(f"Data Loaded: {X.shape[0]} samples, {X.shape[1]} genes")

# --- PART 1: HIERARCHICAL CLUSTERING (VISUALIZATION) ---
print("\nGenerating Dendrogram...")

# Calculate Linkage (Ward's method minimizes variance within clusters)
linked = linkage(X_scaled, method='ward')

plt.figure(figsize=(12, 6))
plt.title("Sample Dendrogram (Hierarchical Clustering)")
plt.xlabel("Samples")
plt.ylabel("Euclidean Distance")

# Plot Dendrogram
# We color the branches to help see the main clusters
dendrogram(linked,
           orientation='top',
           distance_sort='descending',
           show_leaf_counts=True,
           no_labels=True) # Labels hidden for clarity (300+ is too messy)

# Draw a red line to visualize a potential "Cut" height for outliers
# (Adjust 'y=150' based on the plot height you see)
plt.axhline(y=150, c='r', ls='--', lw=2)
plt.text(0, 155, 'Potential Outlier Threshold', color='red')
plt.show()

# --- PART 2: AUTOMATED OUTLIER REMOVAL (STATISTICAL) ---
# We use "Inter-Sample Correlation".
# Theory: A good sample correlates well with others. A bad sample (outlier) has low correlation.

print("Calculating Sample Connectivity (Z-Scores)...")

# 1. Calculate Correlation Matrix (Sample vs Sample)
# We transpose X.T because .corr() calculates correlation between columns
corr_matrix = X.T.corr()

# 2. Calculate Mean Correlation for each sample
mean_connectivity = corr_matrix.mean()

# 3. Calculate Z-Scores
# (Value - Mean) / Std_Dev
z_scores = (mean_connectivity - mean_connectivity.mean()) / mean_connectivity.std()

# 4. Define Threshold (Standard is Z < -2.0 or -2.5)
# This means the sample is "2 Standard Deviations worse" than the average sample.
outlier_threshold = -2.5

outliers = z_scores[z_scores < outlier_threshold]

print(f"\nFound {len(outliers)} Outliers (Z < {outlier_threshold}):")
if len(outliers) > 0:
    print(outliers.index.tolist())
else:
    print("No statistical outliers found.")

# --- PART 3: REMOVE OUTLIERS & SAVE ---
if len(outliers) > 0:
    # Drop the outlier rows from the original dataframe
    df_clean = df.drop(index=outliers.index)
    print(f"\nRemoved {len(outliers)} samples.")
    print(f"New Shape: {df_clean.shape}")
    
    # Save
    df_clean.to_csv(OUT)
    print(f"Saved clean dataset to: '{OUT}'")
else:
    print("\nNo samples removed. Original data kept.")
    df.to_csv(OUT)  # save anyway so later steps have a consistent file