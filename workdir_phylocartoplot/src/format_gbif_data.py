import sys
import pandas as pd
import os
import re

#csv_file = '../input/gbif_coffea_ex3.csv'
#nodes_file = '../input/node_names.csv'
    
def format_gbif(csv_file, nodes_file):
    df_gbif = pd.read_csv(csv_file, usecols=['specimen_id', 'longitude', 'latitude'])
    df_node = pd.read_csv(nodes_file)

    #print(df_node)

    # Function to extract part of the specimen_id for matching
    def extract_name(specimen_id):
        # Use regex to extract pattern (e.g., everything before the underscore)
        return re.sub(r'_\d+', '', specimen_id)

    # Apply extraction function to both DataFrames
    df_gbif['key'] = df_gbif['specimen_id'].apply(extract_name)
    df_node['key'] = df_node['node_name'].apply(lambda x: re.sub(r'^C_|_[\dA-Za-z]+$', '', x))

    # Create a dictionary for mapping key to Node Name
    mapping = df_node.set_index('key')['node_name'].to_dict()

    # Map the Node Name into a new column in df_gbif
    df_gbif['node_name'] = df_gbif['key'].map(mapping)



    # Drop the key column (optional)
    df_gbif.drop(columns='key', inplace=True)
    df_gbif = df_gbif.dropna(subset=['node_name'])

    #print(df_gbif)


    df_new = df_gbif[['node_name','longitude', 'latitude']]

    # Renaming columns
    df_new = df_new.rename(columns={'node_name': 'specimen_id'})

    #print(df_new)
        
    base_name, extension = os.path.splitext(csv_file)

    formatted_csv_file = base_name + '_formatted' + extension

    df_gbif.to_csv(formatted_csv_file, index=False)

    print(f"Data saved to {formatted_csv_file}")



if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python format_gbif_data.py <gbif_file> <nodes_file> \n Output in data folder")
        sys.exit(1)
    
    csv_file = sys.argv[1]
    nodes_file = sys.argv[2]
    
    format_gbif(csv_file, nodes_file)


#csv_file = '../data/gbif_coffea_5years.csv'
#nodes_file = '../data/node_names.csv'