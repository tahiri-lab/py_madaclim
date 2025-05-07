import matplotlib.pyplot as plt
from matplotlib.patches import ConnectionPatch
import matplotlib.colors as mcolors
import matplotlib.lines as mlines
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
            child_x_pos, child_y_pos, y_start = calc_node_positions(
                child, x_start, x_end, y_start, y_step
            )
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
    if node.name in gps["specimen_id"].values:
        color = gps[gps["specimen_id"] == node.name]["color"].values[0]

    ax.plot(
        x,
        y,
        "o",
        markersize=8,
        markerfacecolor=color,
        markeredgewidth=2,
        markeredgecolor="black",
    )
    return x, y


def value_to_color(val):
    if val == 0.00:
        return "grey"
    elif 0.00 < val < 0.7:
        # Gradient from orange to yellow
        cmap = plt.get_cmap("viridis")
        norm = mcolors.Normalize(vmin=0.02, vmax=0.06)
        return mcolors.to_hex(cmap(norm(val)))
    elif val == 0.7:
        return "red"
    else:
        return "black"


def value_to_color_v03(val):
    color_map = {
        1: "red",
        2: "blue",
        0.00: "green",
        4: "purple",
        5: "orange",
        6: "cyan",
        0.7: "magenta",
        0.8: "grey",
        0.9: "brown",
    }
    return color_map.get(
        val, "black"
    )  # Default to 'grey' if the value is not in the dictionary


def main():
    parser = argparse.ArgumentParser(description="PhyloCartoPlot Main script")
    parser.add_argument(
        "--nwk", type=str, required=True, help="Path to the file for tree in nwk format"
    )
    parser.add_argument(
        "--gps", type=str, required=True, help="Path to the file for gps coordinates"
    )
    parser.add_argument(
        "--offset", type=str, required=True, help="Path to the file for ajustments"
    )

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
    offsets_dict = pd.Series(
        offsets_df.XOffset.values, index=offsets_df.NodeName
    ).to_dict()

    # Calculate positions for all nodes
    y_step = 1
    calc_node_positions(tree.root, 0, 1, 0, y_step)

    # Create a figure for the subplot
    ax_tree = fig.add_subplot(121)

    gps = pd.read_csv(args.gps)
    gps["color"] = gps["caffeine_percent"].apply(value_to_color)
    text_colors = dict(zip(gps["specimen_id"], gps["color"]))

    # Plot the tree
    Phylo.draw(
        tree,
        do_show=False,
        axes=ax_tree,
        label_func=custom_label,
        label_colors=text_colors,
    )

    for txt in ax_tree.texts:
        txt.set_fontsize(16)

    # ax_tree.set_title("Coffea species with their geolocation per caffeine content", fontsize=18)
    ax_tree.set_frame_on(False)  # Remove the border/frame
    ax_tree.axis("off")  # Turn off the axes, including ticks and labels
    # Set axes limits to verify the data range
    ax_tree.set_xlim(-0.05, 1)
    ax_tree.set_ylim(0, max(node.position[1] for node in tree.get_terminals()) + 2)

    node_positions = {clade.name: clade.position for clade in tree.find_clades()}

    # Generate DataFrame with node coordinates (commented out as unnecessary here)
    rows = []
    for clade in tree.find_clades():
        if clade.is_terminal():
            label = clade.name
            x, y = plot_adjusted_node(
                ax_tree, clade, y_step, offsets_dict, gps
            )  # Adjust offsets if necessary
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
    for _, row in gps.iterrows():
        is_zero = row["caffeine_percent"] == 0.0
        ax2.plot(
            row["longitude"],
            row["latitude"],
            "o",
            markersize=8,
            markerfacecolor=row["color"],
            markeredgewidth=2,
            markeredgecolor="black",
            alpha=0.2 if is_zero else 1.0,  # ← 50% transparent grey dots
        )
    # for index, row in gps.iterrows():
    #     ax2.plot(
    #         row["longitude"],
    #         row["latitude"],
    #         "o",  # Circle marker
    #         markersize=8,  # Same size as in ax.plot
    #         markerfacecolor=row["color"],  # Use the color from the 'color' column
    #         markeredgewidth=2,  # Same edge width
    #         markeredgecolor="black",  # Same edge color
    #     )

    # Add coastlines and country borders for context
    ax2.coastlines(resolution="10m")
    ax2.add_feature(cfeature.LAND)
    ax2.add_feature(cfeature.OCEAN)
    ax2.add_feature(cfeature.COASTLINE)
    ax2.add_feature(cfeature.BORDERS)

    ax2.set_xlabel("Longitude")
    ax2.set_ylabel("Latitude")
    # ax2.set_title("Species Coordinates", fontsize=18)
    # ax2.legend(["0 %dmb caffeine"], loc="upper right")

    legend_handles = [
        mlines.Line2D(
            [], [], color="grey", marker="o", linestyle="None", label="0% caffeine"
        ),
        mlines.Line2D(
            [], [], color="red", marker="o", linestyle="None", label="0.7% caffeine"
        ),
        mlines.Line2D(
            [],
            [],
            color="purple",
            marker="o",
            linestyle="None",
            label="0.02–0.06% caffeine",
        ),
    ]
    ax2.legend(
        handles=legend_handles,
        loc="lower right",  # put it in the bottom-right
        frameon=True,
        fontsize=12,  # increase to 12pt
    )

    import matplotlib as mpl

    # build a ScalarMappable for your viridis gradient
    norm = mpl.colors.Normalize(vmin=0.02, vmax=0.06)
    cmap = plt.get_cmap("viridis")
    sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    # sm.set_array([])              # only needed for older Matplotlib versions

    # after you’ve drawn ax2 …
    cbar = fig.colorbar(
        sm,
        ax=ax2,
        orientation="vertical",  # or 'horizontal'
        pad=0.02,  # spacing between plot and bar
    )
    cbar.set_label("Caffeine %", fontsize=12)
    cbar.ax.tick_params(labelsize=10)

    # --------------------------------------
    # ------------  Line mapping -----------
    # --------------------------------------

    # Group gps DataFrame by ID and create a dictionary of lists of coordinates
    gps_grouped = (
        gps.groupby("specimen_id")[["longitude", "latitude", "color"]]
        .apply(lambda x: list(zip(x["longitude"], x["latitude"], x["color"])))
        .to_dict()
    )
    caff_dict = dict(zip(gps["specimen_id"], gps["caffeine_percent"]))

    for index, row in df.iterrows():
        specimen = row["ID"]
        if specimen in gps_grouped:
            for longitude, latitude, color in gps_grouped[specimen]:
                caf = caff_dict.get(specimen, 0.0)
                lw = 2.0 if caf > 0 else 0.8
                # if the line is grey (0% caffeine) make it 50% opaque, else keep your 0.2
                line_alpha = 0.2 if color == "grey" else 0.2

                con = ConnectionPatch(
                    xyA=row["Coordinates"],
                    coordsA="data",
                    xyB=(longitude, latitude),
                    coordsB="data",
                    axesA=ax_tree,
                    axesB=ax2,
                    color=color,
                    linewidth=lw,  # whatever logic you already have for line thickness
                    linestyle="--",
                    alpha=line_alpha,  # ← use our conditional transparency
                    zorder=2,
                )
                fig.add_artist(con)

    # for index, row in df.iterrows():
    #     specimen = row["ID"]
    #     if specimen in gps_grouped:
    #         # fetch its caffeine %
    #         caf = caff_dict.get(specimen, 0.0)
    #         # pick a thicker line for non-zero caffeine
    #         lw = 2.0 if caf > 0 else 0.8

    #         for longitude, latitude, color in gps_grouped[specimen]:
    #             con = ConnectionPatch(
    #                 xyA=row["Coordinates"],
    #                 coordsA="data",
    #                 xyB=(longitude, latitude),
    #                 coordsB="data",
    #                 axesA=ax_tree,
    #                 axesB=ax2,
    #                 color=color,
    #                 linewidth=lw,         # ← use our dynamic line width
    #                 linestyle="--",
    #                 alpha=0.2,
    #                 zorder=2,
    #             )
    #             fig.add_artist(con)

    # for index, row in df.iterrows():
    #     # Get corresponding list of coordinates and color from gps DataFrame
    #     if row["ID"] in gps_grouped:
    #         species_coords_list = gps_grouped[row["ID"]]
    #         # Create connection patches for each coordinate in the list
    #         for species_coords in species_coords_list:
    #             longitude, latitude, color = (
    #                 species_coords  # Unpack coordinates and color
    #             )
    #             con = ConnectionPatch(
    #                 xyA=row["Coordinates"],
    #                 coordsA="data",
    #                 xyB=(longitude, latitude),
    #                 coordsB="data",
    #                 axesA=ax_tree,
    #                 axesB=ax2,
    #                 color=color,  # Use the corresponding color for each specimen_id
    #                 linewidth=0.8,
    #                 linestyle="--",
    #                 alpha=0.2,
    #                 zorder=2,
    #             )
    #             fig.add_artist(con)

    plt.tight_layout()
    output_file = r"..\images\figure3.svg"
    plt.savefig(output_file, format="svg")
    output_file = r"..\images\figure3.png"
    plt.savefig(output_file, format="png")

    plt.close(fig)


if __name__ == "__main__":
    main()

    # nwk_file = r"..\input\aligned_caffeine_tree.nwk"
    # gps_coords = r'..\tmp\file_w_caffeine.csv'
    # offsets_file = r"..\input\offsets_caff.csv"
# example use
# python .\tree_to_map_caffeine_content.py --nwk nwk_file --gps gps_coords --offset offsets_file
