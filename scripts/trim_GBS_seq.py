from pathlib import Path
import yaml
from numpy import random

# Get ROOT path
ROOT_DIR = Path(__file__).parents[1]

# Names based on config file
with open(ROOT_DIR.joinpath("config.yaml"), "r") as yaml_file:
    config = yaml.safe_load(yaml_file)
    
# Genetic data dirs
GENETIC_DIR = ROOT_DIR.joinpath("data") / config["genetic"]["dir"]["main"]
SNP_DIR = GENETIC_DIR / config["genetic"]["dir"]["gbs_full"]
SNP_TRIM_DIR = GENETIC_DIR / config["genetic"]["dir"]["gbs_trimmed"]

# Get all sequencing file names in data/GBS dir
gbs_filenames = [file.name for file in SNP_DIR.iterdir() if SNP_DIR.is_dir() and file.is_file()]

# Trimmed concat file
outfile_concat = SNP_TRIM_DIR.joinpath(config["genetic"]["files"]["trimmed_concat_fasta"])

# Number of bases to keep
trim_len = config["params"]["bp_to_keep"]

# Save to individual _trimmed tag
def trim_to_unique(in_dir, out_dir):
    # Create out_dir if needed
    if not out_dir.exists():
        out_dir.mkdir()
    # Sanity check
    if in_dir.is_dir():
        # Loop through all files and trim starting from a random sequence
        for full_fasta in in_dir.iterdir():
            if full_fasta.is_file():
                with open(full_fasta, "r") as fasta_in:
                    # Catch the sequence id
                    id = fasta_in.readline().strip("\n")
                    seq_full = ("".join([line.rstrip() for line in fasta_in]))
                    # Set random start from at least trimming length away from sequence end
                    random.seed(0)
                    trim_start = random.randint(0, (len(seq_full)-trim_len))
                    seq_trimmed = seq_full[trim_start:trim_start+trim_len]
                    print(seq_trimmed)
                    # Save to new dir/file with _trimmed tag
                    with open(out_dir / (full_fasta.stem + "_trimmed" + full_fasta.suffix), "w") as out_file:
                        out_file.write(f"{id}")
                        # Newline every 80th character
                        for i, nucleotide in enumerate(seq_trimmed):
                            if i % 80 == 0 :
                                out_file.write(f"\n")
                            else :
                                out_file.write(nucleotide)

# Concatenate to single fasta file
def trim_to_concat(trimmed_dir, outfile):
    with open(trimmed_dir / outfile, "w") as out_file:
        for trimmed_fasta in trimmed_dir.iterdir():
            if trimmed_fasta.is_file():
                with open(trimmed_fasta, "r") as in_file:
                    print(in_file.readline())   #TODO CONCAT TO single FILE

if __name__ == "__main__":
    # Get trimmed fasta from full fasta
    trim_to_unique(SNP_DIR, SNP_TRIM_DIR)

    # Get concat fasta from all trimmed fasta
    # trim_to_concat(SNP_TRIM_DIR, outfile_concat)  #TODO VERIFY OUTPUT WHEN WORKING FUNC