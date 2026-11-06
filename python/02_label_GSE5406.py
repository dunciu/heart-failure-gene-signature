import pandas as pd
import numpy as np

# --- STEP 1: LOAD THE FILES ---
print("Loading data files...")

# Load the Raw Expression Matrix (Genes are rows, Samples are columns)
# This file was created by your R script
df_expression = pd.read_csv("data/GSE5406_Expression_Matrix.csv", index_col=0)

# Load the Metadata (Contains the labels)
# This file was downloaded by your R script
df_meta = pd.read_csv("data/GSE5406_Metadata.csv")

print(f"Expression Data Loaded: {df_expression.shape[0]} genes, {df_expression.shape[1]} samples")
print(f"Metadata Loaded: {df_meta.shape[0]} samples")

# --- STEP 2: TRANSPOSE THE DATA ---
# Machine Learning expects: Rows = Patients, Columns = Genes
df_transposed = df_expression.T

# --- STEP 3: CLEAN UP FILENAMES ---
# Ensure the index exactly matches the GSM IDs in the metadata
def clean_filename(filename):
    # Split by '_' and take the first part (GSMxxxxxx)
    return str(filename).split('_')[0]

# Apply the cleaning function to the index
clean_index = [clean_filename(x) for x in df_transposed.index]
df_transposed.index = clean_index

print(f"\nFirst 5 Cleaned IDs: {df_transposed.index[:5].tolist()}")

# --- STEP 4: MAP TARGETS TO DISEASE NAMES ---
# Create a dictionary mapping ID -> Disease Name using 'characteristics_ch1'
id_to_disease = dict(zip(df_meta['geo_accession'], df_meta['characteristics_ch1']))

# Map the dictionary to the dataframe index
df_transposed['Target_Label'] = df_transposed.index.map(id_to_disease)

# Check unique labels found (just to be safe)
print("\nUnique Disease Labels Found:")
print(df_transposed['Target_Label'].unique())

# --- STEP 5: ENCODE TO NUMBERS (binary: 0 = non-failing, 1 = heart failure) ---
# Ischemic and idiopathic dilated cardiomyopathy are pooled into one 'failing' class.
# (Output files keep the historical '_Multiclass' name.)
# We use the EXACT strings from the metadata file
label_mapping = {
    'normally functioning myocardium from unused donor heart': 0,             # Control
    'systolic heart failure due to ischemic cardiomyopathy': 1,               # Ischemic -> failing
    'systolic heart failure due to idiopathic dilated cardiomyopathy': 1     # Dilated  -> failing
}

# Apply the numeric mapping
df_transposed['Target'] = df_transposed['Target_Label'].map(label_mapping)

# --- STEP 6: VERIFY SUCCESS ---
# Check if any rows failed to map (NaN)
nan_count = df_transposed['Target'].isna().sum()

if nan_count > 0:
    print(f"\nCRITICAL ERROR: {nan_count} samples could not be mapped!")
    print(df_transposed[df_transposed['Target'].isna()]['Target_Label'])
else:
    print("\nSUCCESS: All samples mapped to 0 (non-failing) / 1 (failing).")
    print("Class Distribution:")
    print(df_transposed['Target'].value_counts().sort_index())

# --- STEP 7: SAVE THE FINAL FILE ---
# Drop the text label column (we only need the numbers for Machine Learning)
final_df = df_transposed.drop(columns=['Target_Label'])

# Move 'Target' to the first column for readability
cols = ['Target'] + [col for col in final_df.columns if col != 'Target']
final_df = final_df[cols]

output_filename = "data/GSE5406_Raw_Multiclass.csv"
final_df.to_csv(output_filename)

print(f"\nSAVED: '{output_filename}'")
print(f"Final Dataset Shape: {final_df.shape}")
print("You can now use this file for your Boxplots and Machine Learning.")