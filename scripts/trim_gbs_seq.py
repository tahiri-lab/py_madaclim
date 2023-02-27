from pathlib import Path
import argparse
import yaml

from Bio.SeqIO.FastaIO import SimpleFastaParser
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

def indir_path(in_path):
    in_dir = Path(in_path)
    if in_dir.is_dir():
        return in_dir
    else:
        raise argparse.ArgumentTypeError(f"{in_dir} is not a valid path")
    
def outdir_path(out_path):
    out_dir = Path(out_path)
    # Create out_dir if needed
    if not out_dir.exists():
        out_dir.mkdir()
        return out_dir
    else:
        return out_dir
    
def get_seq_len(fasta):
    count_seq = 0
    length_seq = []
    # Parse single or multi_fasta
    with open(fasta, "r") as f:
        for title, seq in SimpleFastaParser(f):
            count_seq += 1
            length_seq.append(len(seq))
    min_seq_length = min(length_seq)
    
    return count_seq, min_seq_length

def get_start_pos(startpos_percent, seq_length):
    startpos = (startpos_percent / 100) * seq_length
    # Determine length of sequence from startpos to end of sequence
    start_to_end_length = seq_length - startpos
    return startpos, start_to_end_length
    
#!TODO REMOVE NOT USED
# def trim_len(fasta, trim_length):
#     num_seq, seq_length = get_seq_len(fasta)
#     if trim_length > seq_length:
#         raise ValueError(
#             f"Length to trim must be smaller than the {num_seq} sequence(s) of length {seq_length}"
#         )
#     else:
#         return trim_length

def _build_arg_parser():
    # For defaults
    config, SNP_DIR, SNP_TRIM_DIR = get_default_fasta_dir()

    parser = argparse.ArgumentParser(
    description = 
    """
    Takes multiple pre-aligned fasta files (or a single multifasta)
    and concatenates them in a single file.\n 
    Then trims all the sequences to a given length from either a random seed or a given start position.
    """,
    epilog = "Output can be used for tests in phylogenetic tree building with smaller sets for lower computing time"
    )

    parser.add_argument(
        "in_dir",
        nargs='?',
        type=indir_path,
        default=SNP_DIR,
        help="Directory name that contains the fasta files"
    )

    parser.add_argument(
        "out_dir",
        nargs="?",
        type=outdir_path,
        default=SNP_TRIM_DIR,
        help="Directory name that contains the trimmed fasta files"
    )

    parser.add_argument(
        "-l", "--length",
        nargs="?",
        type=int,
        default=config["params"]["bp_to_keep"],
        help="Length(int) of the sequence to trim. [Default = %(default)s]."
    )

    parser.add_argument(
        "-sp", "--startpos",
        nargs="?",
        metavar="N",
        type=int,
        choices=range(0,100),
        default=None,
        help="Start position from pth percentile of sequence instead of random. Choose from 0 to 99. [Default = %(default)s for random]"
    )

    parser.add_argument(
        "-c", "--concat",
        action="store_true",
        help="Add if input fasta is a single file that is already concatenated"
    )
    
    return parser
    
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
    args = _build_arg_parser()
    args.parse_args()
    # Get trimmed fasta from full fasta
    # trim_to_unique(SNP_DIR, SNP_TRIM_DIR)

    # # Get concat fasta from all trimmed fasta
    # trim_to_concat(SNP_TRIM_DIR, outfile_concat)
    