import pandas as pd
import numpy as np

# --- STEP 1: LOAD THE FILES ---
print("Loading data files...")

# Load the Raw Expression Matrix (Genes are rows, Samples are columns)
# This file was created by your R script
df_expression = pd.read_csv("data/GSE57338_Cleaned_Normalized.csv", index_col=0)

# Load the Metadata (Contains the labels)
# This file was downloaded by your R script
df_meta = pd.read_csv("data/GSE57338_Metadata.csv")

print(f"Expression Data Loaded: {df_expression.shape[0]} genes, {df_expression.shape[1]} samples")
print(f"Metadata Loaded: {df_meta.shape[0]} samples")

# --- STEP 2: TRANSPOSE THE DATA ---
# Machine Learning expects: Rows = Patients, Columns = Genes
df_transposed = df_expression.T

# --- STEP 3: CLEAN UP FILENAMES ---
# The index currently looks like: "GSM1379812_C02065_HuGene1.1ST.CEL"
# The metadata ID looks like: "GSM1379812"
# We must strip the extra junk so they match.

def clean_filename(filename):
    # Split by '_' and take the first part (GSMxxxxxx)
    # This turns "GSM1379812_C02065..." into just "GSM1379812"
    return str(filename).split('_')[0]

# Apply the cleaning function to the index
clean_index = [clean_filename(x) for x in df_transposed.index]
df_transposed.index = clean_index

print(f"\nFirst 5 Cleaned IDs: {df_transposed.index[:5].tolist()}")

# --- STEP 4: MAP TARGETS TO DISEASE NAMES ---
# Create a dictionary mapping ID -> Disease Name
# { 'GSM1379812': 'heart left ventricle, ischemic', ... }
id_to_disease = dict(zip(df_meta['geo_accession'], df_meta['source_name_ch1']))

# Map the dictionary to the dataframe index
df_transposed['Target_Label'] = df_transposed.index.map(id_to_disease)

# Check unique labels found (just to be safe)
print("\nUnique Disease Labels Found:", df_transposed['Target_Label'].unique())

# --- STEP 5: ENCODE TO NUMBERS (binary: 0 = non-failing, 1 = heart failure) ---
# Ischemic and idiopathic dilated cardiomyopathy are pooled into one 'failing' class.
# (Output files keep the historical '_Multiclass' name.)
# We use the EXACT strings from your metadata file
label_mapping = {
    'heart left ventricle, non-failing': 0,             # Control
    'heart left ventricle, ischemic': 1,                # Ischemic -> failing
    'heart left ventricle, idiopathic dilated CMP': 1  # Dilated  -> failing
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
# Drop the text label column (we only need the numbers)
final_df = df_transposed.drop(columns=['Target_Label'])

# Move 'Target' to the first column for readability
cols = ['Target'] + [col for col in final_df.columns if col != 'Target']
final_df = final_df[cols]

output_filename = "data/GSE57338_Raw_Multiclass.csv"
final_df.to_csv(output_filename)

print(f"\nSAVED: '{output_filename}'")
print(f"Final Dataset Shape: {final_df.shape}")
print("You can now use this file for your Boxplots and Machine Learning.")