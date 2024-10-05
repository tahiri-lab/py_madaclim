import matplotlib.pyplot as plt
from matplotlib.patches import ConnectionPatch
import matplotlib.colors as mcolors  
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import pandas as pd
from Bio import Phylo
import argparse

def custom_label(clade):
  if clade.is_terminal():
    return clade.name
  else:
    return None

def calc_node_positions(tree, x_start, x_end, y_start, y_step):
  if tree.is_terminal():
    x_pos = x_start
    y_pos = y_start
    y_start += y_step
  else:
    x_pos = (x_start + x_end) / 2
    y_pos = y_start

    child_y_start = y_start
    for child in tree.clades:
      child_x_pos, child_y_pos, y_start = calc_node_positions(child, x_start, x_end, y_start, y_step)
      x_start = child_x_pos

    y_pos = (y_start + child_y_start) / 2

  tree.position = (x_pos, y_pos)
  return x_pos, y_pos, y_start

def get_x_offset(node_name, offsets_dict):
  return offsets_dict.get(node_name, 0)  # default offset is 0 if not found

def plot_adjusted_node(ax, node, y_offset, offsets_dict, gps):
    x, y = node.position
    x_offset = get_x_offset(node.name, offsets_dict)
    x += x_offset
    y += y_offset

    # Default color if node name is not found
    color = "grey"

    # Check if node name is in the gps DataFrame and set color accordingly
    if node.name in gps['specimen_id'].values:
        color = gps[gps['specimen_id'] == node.name]['color'].values[0]
        
    ax.plot(x, y, 'o', markersize=8, markerfacecolor=color, markeredgewidth=2, markeredgecolor="black")
    return x, y

def value_to_color(val):
    if val == 0.00:
        return 'red'
    elif 0.00 < val < 0.7:
        # Gradient from orange to yellow
        cmap = plt.get_cmap('viridis')  # Yellow to red colormap
        norm = mcolors.Normalize(vmin=0.02, vmax=0.06)
        return mcolors.to_hex(cmap(norm(val)))
    elif val == 0.7:
        return 'green'
    else:
        return 'grey'

def main():
  parser = argparse.ArgumentParser(description='PhyloCartoPlot Main script')
  parser.add_argument('--nwk', type=str, required=True, help='Path to the file for tree in nwk format')
  parser.add_argument('--gps', type=str, required=True, help='Path to the file for gps coordinates')
  parser.add_argument('--offset', type=str, required=True, help='Path to the file for ajustments')

  args = parser.parse_args()
  
    # Create a new map with PlateCarree projection
  fig = plt.figure(figsize=(26, 11))

  # --------------------------------------
  # ------------  Phylogenetic MAP -------
  # --------------------------------------

  # Load the tree
  tree = Phylo.read(args.nwk, "newick")

  # Load offsets from CSV
  offsets_df = pd.read_csv(args.offset)
  offsets_dict = pd.Series(offsets_df.XOffset.values, index=offsets_df.NodeName).to_dict()

  # Calculate positions for all nodes
  y_step = 1
  calc_node_positions(tree.root, 0, 1, 0, y_step)

  # Create a figure for the subplot
  ax_tree = fig.add_subplot(121)

  gps = pd.read_csv(args.gps)
  gps['color'] = gps['caffeine_percent'].apply(value_to_color)
  text_colors = dict(zip(gps['specimen_id'], gps['color']))

  # Plot the tree
  Phylo.draw(tree, do_show=False, axes=ax_tree, label_func=custom_label, label_colors=text_colors)
  ax_tree.set_title("Coffea species with their geolocation per caffeine content", fontsize=18)

  # Set axes limits to verify the data range
  ax_tree.set_xlim(0, 1)
  ax_tree.set_ylim(0, max(node.position[1] for node in tree.get_terminals()) + 2)

  node_positions = {clade.name: clade.position for clade in tree.find_clades()}

  # Generate DataFrame with node coordinates (commented out as unnecessary here)
  rows = []
  for clade in tree.find_clades():
    if clade.is_terminal():
      label = clade.name
      x, y = plot_adjusted_node(ax_tree, clade, y_step, offsets_dict, gps)  # Adjust offsets if necessary
      rows.append([label, (x, y)])

  # Create DataFrame with node coordinates (commented out as unnecessary here)
  df = pd.DataFrame(rows, columns=["ID", "Coordinates"])

  # --------------------------------------
  # ------------  GRAPH MAP --------------
  # --------------------------------------


  # Create subplot 2 with the map plot
  ax2 = fig.add_subplot(122, projection=ccrs.PlateCarree())
  extent = [43, 51, -27, -11]
  ax2.set_extent(extent)

  # Plot points from GPS dataframe on the map
  for index, row in gps.iterrows():
      ax2.plot(
          row["longitude"], row["latitude"], 
          'o',  # Circle marker
          markersize=8,  # Same size as in ax.plot
          markerfacecolor=row['color'],  # Use the color from the 'color' column
          markeredgewidth=2,  # Same edge width
          markeredgecolor='black'  # Same edge color
      )

  # Add coastlines and country borders for context
  ax2.coastlines(resolution='10m')
  ax2.add_feature(cfeature.LAND)
  ax2.add_feature(cfeature.OCEAN)
  ax2.add_feature(cfeature.COASTLINE)
  ax2.add_feature(cfeature.BORDERS)

  ax2.set_xlabel("Longitude")
  ax2.set_ylabel("Latitude")
  ax2.set_title("Species Coordinates", fontsize=18)
  ax2.legend(["0 %dmb caffeine"], loc='upper right')

  # --------------------------------------
  # ------------  Line mapping -----------
  # --------------------------------------

  # Group gps DataFrame by ID and create a dictionary of lists of coordinates
  gps_grouped = gps.groupby('specimen_id')[['longitude', 'latitude', 'color']].apply(
      lambda x: list(zip(x['longitude'], x['latitude'], x['color']))).to_dict()
  #print(gps_grouped)
  for index, row in df.iterrows():
      # Get corresponding list of coordinates and color from gps DataFrame
      if row['ID'] in gps_grouped:
          species_coords_list = gps_grouped[row['ID']]
          # Create connection patches for each coordinate in the list
          for species_coords in species_coords_list:
              longitude, latitude, color = species_coords  # Unpack coordinates and color
              con = ConnectionPatch(
                  xyA=row['Coordinates'], coordsA="data",
                  xyB=(longitude, latitude), coordsB="data",
                  axesA=ax_tree, axesB=ax2,
                  color=color,  # Use the corresponding color for each specimen_id
                  linewidth=1, linestyle="--", alpha=0.5,
                  zorder=2
              )
              fig.add_artist(con)

  plt.tight_layout()
  output_file = r'..\images\figure2.svg'  
  plt.savefig(output_file, format='svg')
  output_file = r'..\images\figure2.png'  
  plt.savefig(output_file, format='png')
  
  plt.close(fig)


if __name__ == "__main__":
  main()
  

  #nwk_file = r"..\input\aligned_caffeine_tree.nwk"
  #gps_coords = r'..\tmp\file_w_caffeine.csv'
  #offsets_file = r"..\input\offsets_caff.csv"
# example use
  # python .\tree_to_map_caffeine_content.py --nwk nwk_file --gps gps_coords --offset offsets_file

