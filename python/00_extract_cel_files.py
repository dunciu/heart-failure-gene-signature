import tarfile

# The file you downloaded
file_path = "data/GSE57338_RAW.tar"

# Open the file
with tarfile.open(file_path, 'r') as tar:
    # Extract all contents to a specific folder
    tar.extractall(path="data/extracted_data")
    
import os
import gzip
import shutil

source_folder = "data/extracted_data"
destination_folder = "data/unzipped_cel_files"

os.makedirs(destination_folder, exist_ok=True)

# Loop through all files
for filename in os.listdir(source_folder):
    if filename.endswith(".gz"):
        # Construct full file paths
        gz_path = os.path.join(source_folder, filename)
        new_filename = filename.replace(".gz", "") # Remove .gz extension
        extract_path = os.path.join(destination_folder, new_filename)

        print(f"Unzipping {filename}...")
        
        # Read the compressed file and write the uncompressed version
        with gzip.open(gz_path, 'rb') as f_in:
            with open(extract_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

print("Done! Files are in:", destination_folder)