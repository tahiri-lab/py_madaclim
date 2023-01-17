import os

# Full fasta data
gbs_dir = "./data/GBS"
full_fasta_filenames = os.listdir(gbs_dir)

# Trimmed fasta data
gbs_trimmed = "./data/GBS-trimmed"


# Looping through all the files
for filename in full_fasta_filenames:
    file = os.path.join(gbs_dir, filename)
    if os.path.isfile(file):    #filecheck
        
        # Read each files, save id(header) and first 1000 bp
        with open(file) as f:
            fasta_id = f.readline().strip('\n')
            fasta_seq1000 = (''.join([line.rstrip() for line in f])[:1000])
            #TODO write to file in a single concatenated fasta
        