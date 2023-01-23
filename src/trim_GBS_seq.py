from definitions import config,SNP_DIR, SNP_TRIM_DIR  # Get global vars
from pathlib import Path

# Number of bases to keep
seq_len = config["genetic"]["bp_to_keep"]

# Get all sequencing file names in data/GBS dir
gbs_filenames = [file.name for file in SNP_DIR.iterdir() if SNP_DIR.is_dir() and file.is_file()]

# Save to individual _trimmed tag
def trim_to_unique(in_dir, out_dir):
    for full_fasta in in_dir.iterdir():
        if full_fasta.is_file():
            with open(full_fasta, "r") as fasta_in:
                id = fasta_in.readline().strip("\n") # Catch the sequence id
                seq_trimmed = ("".join([line.rstrip() for line in fasta_in])[:seq_len+1])  # Get first bp up to seq_len
                # Save to new dir/file with _trimmed tag
                with open(out_dir / (full_fasta.stem + "_trimmed" + full_fasta.suffix), "w") as out_file:
                    out_file.write(f"{id}")
                    # Newline every 80th character
                    for i, nucleotide in enumerate(seq_trimmed):
                        if i % 80 == 0 :
                            out_file.write(f"\n")
                        else :
                            out_file.write(nucleotide)

def trim_to_concat(trimmed_dir, outfile):
    with open(trimmed_dir / outfile, "r"):
        for trimmed_fasta in trimmed_dir.iterdir():
            if trimmed_fasta.is_file():
                print(trimmed_fasta.read())

if __name__ == "__main__":
    # Get trimmed fasta from full fasta
    # trim_to_unique(SNP_DIR, SNP_TRIM_DIR)

    # Get concat fasta from all trimmed fasta
            