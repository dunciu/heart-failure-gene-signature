import pandas as pd

# --- STEP 1: LOAD THE DATASETS ---
print("Loading datasets...")
df_train = pd.read_csv("data/GSE57338_Top5000_Var.csv", index_col=0)
df_test = pd.read_csv("data/GSE5406_Raw_Multiclass.csv", index_col=0)

train_genes = list(df_train.columns.drop('Target'))

# --- STEP 2: KEEP ALL 5000 FOR TRAINING ---
# We don't shrink the training set at all.
df_train_final = df_train[['Target'] + train_genes]

# --- STEP 3: PAD THE TESTING SET ---
print("Aligning and padding the external test set...")

# Create a brand new, empty dataframe for the test set
df_test_final = pd.DataFrame(index=df_test.index)

# 1. Copy over the Target column
df_test_final['Target'] = df_test['Target']

# 2. Loop through all 5000 training genes
missing_count = 0
found_count = 0

for gene in train_genes:
    if gene in df_test.columns:
        # If the gene exists in GSE5406 (e.g., the 3400), copy the real data
        df_test_final[gene] = df_test[gene]
        found_count += 1
    else:
        # If the gene is missing, create the column but fill it with zeros
        df_test_final[gene] = 0.0
        missing_count += 1

print(f"\nTraining set features: {len(train_genes)}")
print(f"Testing set features found: {found_count}")
print(f"Testing set features missing (padded with 0): {missing_count}")

# Verify shapes match (Target + 5000 genes = 5001 columns)
print(f"\nFinal Training Shape: {df_train_final.shape}")
print(f"Final Testing Shape:  {df_test_final.shape}")

# --- STEP 4: SAVE THE DATASETS ---
df_train_final.to_csv("data/GSE57338_Train_5000.csv")
df_test_final.to_csv("data/GSE5406_Test_Padded.csv")

print("\nSUCCESS! Saved 'GSE57338_Train_5000.csv' and 'GSE5406_Test_Padded.csv'")