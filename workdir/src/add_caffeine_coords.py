import pandas as pd
import argparse

#input_file = '../data/gbif_coffea_ex3_formatted.csv'
#caffeine_file = '../data/no_caffeine_nodes_w_specimen.csv'

def add_caffeine(input_file, caffeine_file):
    
    gbif_df = pd.read_csv(input_file)
    node_names_df = pd.read_csv(caffeine_file)

    # Merge the two DataFrames based on the matching specimen_id columns
    # Assuming 'specimen_id' in gbif_df and 'Specimen_id' in node_names_df are the corresponding columns
    merged_df = pd.merge(gbif_df, node_names_df[['Species_name', 'caffeine_percent']], 
                        left_on='Node Name', right_on='Species_name', how='left')

    # Drop the 'Specimen_id' column if you no longer need it after the merge
    merged_df = merged_df.drop(columns=['Species_name', 'specimen_id'])
    merged_df = merged_df.rename(columns={'Node Name': 'specimen_id'})
    merged_df = merged_df[['specimen_id', 'longitude', 'latitude', 'caffeine_percent']]
    merged_df = merged_df.dropna(subset=['caffeine_percent'])


    output_file = input_file.replace('formatted', 'w_caffeine')

    # Save the merged DataFrame to a new CSV file
    merged_df.to_csv(output_file, index=False)

    # Print the output file name
    print(f"Data saved to: {output_file}")
    
    
def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(description='Process and merge data files.')
    parser.add_argument('input_file', type=str, help='Path to the input GBIF data CSV file.')
    parser.add_argument('caffeine_file', type=str, help='Path to the caffeine content CSV file.')
    args = parser.parse_args()

    # Process the data
    add_caffeine(args.input_file, args.caffeine_file)

if __name__ == "__main__":
    main()