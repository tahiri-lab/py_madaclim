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
FASTA_DIR = GENETIC_DIR / config["genetic"]["dir"]["gbs_full"]

# Full fasta concat file name
outfile_concat = FASTA_DIR.joinpath(config["genetic"]["files"]["full_concat_fasta"])

# Concatenate to single fasta file
def trim_to_concat(in_dir, outfile):
    # Sorted list of all fasta files to concat
    trimmed_fastas = [fastafile for fastafile in sorted(in_dir.iterdir()) if in_dir.is_dir() and fastafile.is_file()]
    
    # Sanity check + Write to .fasta (output) 
    if in_dir.is_dir():
        with open(in_dir / outfile, "w") as concat_file:
            
            # Loop through all individual fasta
            for trimmed_fasta in trimmed_fastas:
                if trimmed_fasta.is_file(): # Sanity check
                    with open(trimmed_fasta, "r") as in_file:
                        data = in_file.read().split("\n")
                        id = data[0]   # Catch the sequence id and write to file
                        seq_full = ("".join(data[1:]))  # Catch the whole sequence
                        
                        # Write header and skip line
                        concat_file.write(f"{id}\n")
                        
                        # Newline every 60th character                        
                        for i, nucleotide in enumerate(seq_full):
                            if i % 60 == 0 and i > 0:
                                concat_file.write(f"\n{nucleotide}")
                            else:
                                concat_file.write(nucleotide)
                    concat_file.write(f"\n")                
                else:
                    raise FileNotFoundError
    else:
        raise NotADirectoryError
                        
if __name__ == "__main__":
    # Get concat fasta from all full fasta files
    trim_to_concat(FASTA_DIR, outfile_concat)
    