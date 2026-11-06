import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.feature_selection import f_classif, SelectKBest
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

# 1. Load the top-variance genes from step 04
# (originally read 'GSE57338_Raw_NoNorm.csv', a file no script in this repo produces)
df = pd.read_csv("data/GSE57338_Top5000_Var.csv", index_col=0)

y = df['Target']
X = df.drop(columns=['Target'])

# 2. Run ANOVA (Supervised Feature Selection)
# Which genes have the most distinct means between failing and non-failing hearts?
print("Running ANOVA to find the most discriminating genes...")
selector = SelectKBest(f_classif, k=50) # Select Top 50 Genes
X_new = selector.fit_transform(X, y)

# Get the names of the best genes
selected_indices = selector.get_support(indices=True)
selected_genes = X.columns[selected_indices]
X_best = pd.DataFrame(X_new, columns=selected_genes, index=X.index)

print(f"Selected Top 50 Genes: {selected_genes.tolist()[:5]} ...")

# 3. Visualize with PCA/t-SNE using ONLY these 50 genes
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X_best)

# PCA
pca = PCA(n_components=2)
pca_coords = pca.fit_transform(X_scaled)
df_pca = pd.DataFrame(data=pca_coords, columns=['PC1', 'PC2'], index=X.index)
df_pca['Target'] = y

# t-SNE
tsne = TSNE(n_components=2, perplexity=30, random_state=42)
tsne_results = tsne.fit_transform(X_scaled)
df_tsne = pd.DataFrame(data=tsne_results, columns=['tSNE1', 'tSNE2'], index=X.index)
df_tsne['Target'] = y

# --- PLOTTING ---
fig, axes = plt.subplots(1, 2, figsize=(18, 8))
palette = 'viridis'

# PCA Plot
sns.scatterplot(ax=axes[0], x='PC1', y='PC2', hue='Target', data=df_pca, 
                palette=palette, s=100, alpha=0.8, edgecolor='k')
axes[0].set_title("PCA (Top 50 ANOVA Genes)\nFocused on Differences", fontsize=14)
axes[0].grid(True, alpha=0.3)

# t-SNE Plot
sns.scatterplot(ax=axes[1], x='tSNE1', y='tSNE2', hue='Target', data=df_tsne, 
                palette=palette, s=100, alpha=0.8, edgecolor='k')
axes[1].set_title("t-SNE (Top 50 ANOVA Genes)", fontsize=14)
axes[1].grid(True, alpha=0.3)

plt.show()

# Save the best genes for your Machine Learning model later
df_best = pd.concat([y, X_best], axis=1)
df_best.to_csv("data/GSE57338_Top50_ANOVA.csv")
print("Saved Top 50 discriminating genes to 'GSE57338_Top50_ANOVA.csv'")