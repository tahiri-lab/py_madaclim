from pathlib import Path
import yaml
from numpy import random

def get_default_fasta_dir():
    # Get ROOT path
    ROOT_DIR = Path(__file__).parents[1]

    # Package src dir
    SRC_DIR = ROOT_DIR / "src"

    # Names based on config file
    with open(SRC_DIR.joinpath("config.yaml"), "r") as yaml_file:
        config = yaml.safe_load(yaml_file)
        
    # Genetic data dirs
    GENETIC_DIR = ROOT_DIR.joinpath("data") / config["genetic"]["dir"]["main"]
    SNP_DIR = GENETIC_DIR / config["genetic"]["dir"]["gbs_full"]
    SNP_TRIM_DIR = GENETIC_DIR / config["genetic"]["dir"]["gbs_trimmed"]

    return config, SNP_DIR, SNP_TRIM_DIR

def dir_path(in_path):
    path = Path(in_path)
    if path.is_dir():
        return path
    else:
        raise argparse.ArgumentTypeError(f"{path} is not a valid path")
    
def get_max_trim_length(startpos):
    pass

config, SNP_DIR, SNP_TRIM_DIR = get_default_fasta_dir()

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
    
    # Loop through all files and trim starting from a random sequence
    if in_dir.is_dir(): #Sanity check    
        for full_fasta in in_dir.iterdir():
            if full_fasta.is_file():
                with open(full_fasta, "r") as fasta_in:
                    data = fasta_in.read().split("\n")
                    id = data[0]   # Catch the sequence id
                    seq_full = ("".join(data[1:]))  # Catch the whole sequence
                    
                    # Set random start from at least trimming length away from sequence end
                    random.seed(0)
                    trim_start = random.randint(0, (len(seq_full) - trim_len))
                    seq_trimmed = seq_full[trim_start:trim_start + trim_len]
                    
                    # Save to new dir/file with _trimmed tag
                    out_trimmed = out_dir / (full_fasta.stem + "_trimmed" + full_fasta.suffix)
                    if not out_trimmed.is_file():
                        with open(out_trimmed, "w") as out_file:
                            out_file.write(f"{id}\n")
                            
                            # Newline every 60th character
                            for i, nucleotide in enumerate(seq_trimmed):
                                if i % 60 == 0 and i > 0:
                                    out_file.write(f"\n{nucleotide}")                                
                                else :
                                    out_file.write(nucleotide)
                    else:
                        raise FileExistsError
    else:
        raise NotADirectoryError
    
# Concatenate to single fasta file
def trim_to_concat(trimmed_dir, outfile):
    # Sorted list of all trimmed fasta files to concat
    trimmed_fastas = [fastafile for fastafile in sorted(trimmed_dir.iterdir())]
    
    # Sanity check + Write to trimmed.fasta (output) 
    if trimmed_dir.is_dir():
        with open(trimmed_dir / outfile, "w") as concat_file:
            
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
    # Get trimmed fasta from full fasta
    trim_to_unique(SNP_DIR, SNP_TRIM_DIR)

    # Get concat fasta from all trimmed fasta
    trim_to_concat(SNP_TRIM_DIR, outfile_concat)
    