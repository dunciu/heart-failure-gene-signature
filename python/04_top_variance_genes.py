import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler

# --- STEP 1: LOAD THE CLEAN DATA ---
# Using the multiclass file we just made. 
# (If you have a separate file with outliers removed, update this name!)
df = pd.read_csv("data/GSE57338_Raw_Multiclass.csv", index_col=0)

# Separate Target and Genes
y = df['Target']
X = df.drop(columns=['Target'])

print(f"Original Data: {X.shape[0]} samples, {X.shape[1]} genes")

# --- STEP 2: CALCULATE VARIANCE ---
print("Calculating variance per gene...")
gene_variances = X.var()

# --- STEP 3: SORT AND SELECT TOP 5000 ---
top_n = 5000
sorted_genes = gene_variances.sort_values(ascending=False)
top_genes_list = sorted_genes.head(top_n).index

# Create the filtered dataset
X_top5000 = X[top_genes_list]

print(f"Filtered Data: {X_top5000.shape[0]} samples, {X_top5000.shape[1]} genes")

# --- STEP 4: SAVE THE FILTERED DATASET ---
df_top5000 = pd.concat([y, X_top5000], axis=1)

output_filename = "data/GSE57338_Top5000_Var.csv"   # used by 05 and 06
df_top5000.to_csv(output_filename)
print(f"Saved: '{output_filename}'")

# --- STEP 5: VISUALIZATION CHECK (PCA & t-SNE) ---
print("Running dimensionality reduction for visualization...")

# Standardize (StandardScaler is crucial before PCA/t-SNE)
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X_top5000)

# Run PCA
pca = PCA(n_components=2)
pca_coords = pca.fit_transform(X_scaled)
df_pca = pd.DataFrame(data=pca_coords, columns=['PC1', 'PC2'], index=X.index)
# Convert Target to string so Seaborn treats it as a discrete category, not a gradient
df_pca['Target'] = y.astype(str) 

# Run t-SNE
# GSE57338 has 313 samples, so perplexity=30 is reasonable.
tsne = TSNE(n_components=2, perplexity=30, random_state=42, init='pca', learning_rate='auto')
tsne_coords = tsne.fit_transform(X_scaled)
df_tsne = pd.DataFrame(data=tsne_coords, columns=['tSNE1', 'tSNE2'], index=X.index)
df_tsne['Target'] = y.astype(str)

# Plotting
fig, axes = plt.subplots(1, 2, figsize=(18, 8))

# 0 = non-failing, 1 = heart failure (ischemic + dilated pooled)
palette = sns.color_palette("Set1", n_colors=2)

# PCA Plot
sns.scatterplot(ax=axes[0], x='PC1', y='PC2', hue='Target', data=df_pca, 
                palette=palette, s=80, alpha=0.8, edgecolor='k')
axes[0].set_title(f"PCA (Top {top_n} Variable Genes)", fontsize=14)
axes[0].grid(True, alpha=0.3)

# t-SNE Plot
sns.scatterplot(ax=axes[1], x='tSNE1', y='tSNE2', hue='Target', data=df_tsne, 
                palette=palette, s=80, alpha=0.8, edgecolor='k')
axes[1].set_title(f"t-SNE (Top {top_n} Variable Genes)", fontsize=14)
axes[1].grid(True, alpha=0.3)

plt.tight_layout()
plt.show()