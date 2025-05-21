import matplotlib.pyplot as plt
from matplotlib.patches import ConnectionPatch
import matplotlib.colors as mcolors
import matplotlib.lines as mlines
import matplotlib.cm as cm
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import pandas as pd
import numpy as np
from Bio import Phylo
import argparse

label_map = {
    "C_abbayesii_A601": "C. abbayesii (A601)",
    "C_ambodirianensis_A572": "C. ambodirianensis (A572)",
    "C_ambongensis": "C. ambongensis",
    "C_andrambovatensis_A310": "C. andrambovatensis (A310)",
    "C_ankaranensis_A525": "C. ankaranensis (A525)",
    "C_arenesiana_A403": "C. arenesiana (A403)",
    "C_augagneuri_A966": "C. augagneuri (A966)",
    "C_bernardiniana_MDC": "C. bernardiniana (MDC)",
    "C_bertrandii_A5": "C. bertrandii (A5)",
    "C_bissetiae": "C. bissetiae",
    "C_boinensis": "C. boinensis",
    "C_boiviniana_A980": "C. boiviniana (A980)",
    "C_bonnieri_A535": "C. bonnieri (A535)",
    "C_costei_A956": "C. costei (A956)",
    "C_dolichophylla_A206": "C. dolichophylla (A206)",
    "C_dubardii_A969": "C. dubardii (A969)",
    "C_farafanganensis_A208": "C. farafanganensis (A208)",
    "C_heimii_A516": "C. heimii (A516)",
    "C_homollei_A945": "C. homollei (A945)",
    "C_humbertii_RNF785": "C. humbertii (RNF785)",
    "C_humblotiana_BM19_20": "C. humblotiana (BM19, 20)",
    "C_jumellei_A974": "C. jumellei (A974)",
    "C_kianjavatensis_A602": "C. kianjavatensis (A602)",
    "C_kihansiensis_APD2922": "C. kihansiensis (APD2922)",
    "C_labatii_APD3096": "C. labatii (APD3096)",
    "C_lancifolia_A320": "C. lancifolia (A320)",
    "C_leroyi_A315": "C. leroyi (A315)",
    "C_liaudii_A1013": "C. liaudii (A1013)",
    "C_macrocarpa-PET": "C. macrocarpa (PET)",
    "C_mauritiana_BM17_25": "C. mauritiana (BM17, 25)",
    "C_mauritiana_Makes4": "C. mauritiana (Makes4)",
    "C_mcphersonii_A977": "C. mcphersonii (A977)",
    "C_millotii_A222": "C. millotii (A222)",
    "C_mogeneti_A975": "C. mogeneti (A975)",
    "C_montis-sacri_A321": "C. montis-sacri (A321)",
    "C_myrtifolia_MBR_A9": "C. myrtifolia (MBR, A9)",
    "C_perrieri_A12": "C. perrieri (A12)",
    "C_pervilleana_A957": "C. pervilleana (A957)",
    "C_ratsimamangae_A528": "C. ratsimamangae (A528)",
    "C_resinosa_A8": "C. resinosa (A8)",
    "C_richardii_A575": "C. richardii (A575)",
    "C_sahafaryensis_A978": "C. sahafaryensis (A978)",
    "C_sakarahae_A304": "C. sakarahae (A304)",
    "C_tetragona_A252": "C. tetragona (A252)",
    "C_tsirananae_A515": "C. tsirananae (A515)",
    "C_vatovavyensis_A830": "C. vatovavyensis (A830)",
    "C_vianneyi_A946": "C. vianneyi (A946)",
    "Tricalysia": "C. Tricalysia",
}


def custom_label(clade):
    if clade.is_terminal():
        # return clade.name
        return label_map.get(clade.name, clade.name)
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
        base_cmap = plt.get_cmap("viridis")
        newcolors = base_cmap(np.linspace(0, 0.8, 256))
        viridis_no_yellow = mcolors.ListedColormap(newcolors, name="viridis_no_yellow")
        cmap = viridis_no_yellow
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
    # text_colors = dict(zip(gps["specimen_id"], gps["color"]))
    raw_color_map = dict(zip(gps["specimen_id"], gps["color"]))

    # build a map from *displayed* label → color
    text_colors = {label_map.get(raw, raw): col for raw, col in raw_color_map.items()}

    # Plot the tree
    Phylo.draw(
        tree,
        do_show=False,
        axes=ax_tree,
        label_func=custom_label,
        label_colors=text_colors,
    )

    for txt in ax_tree.texts:
        txt.set_fontsize(24)
        txt.set_fontstyle("italic")

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

    ax2.coastlines(resolution="10m")
    ax2.add_feature(cfeature.LAND)
    ax2.add_feature(cfeature.OCEAN)
    ax2.add_feature(cfeature.COASTLINE)
    ax2.add_feature(cfeature.BORDERS)

    ax2.set_xlabel("Longitude")
    ax2.set_ylabel("Latitude")

    legend_handles = [
        mlines.Line2D(
            [], [], color="grey", marker="o", linestyle="None", label="0% caffeine"
        ),
        mlines.Line2D(
            [],
            [],
            color="purple",
            marker="o",
            linestyle="None",
            label="0.02–0.06% caffeine",
        ),
        mlines.Line2D(
            [], [], color="red", marker="o", linestyle="None", label="0.7% caffeine"
        ),
    ]
    ax2.legend(
        handles=legend_handles,
        loc="lower right",
        frameon=True,
        fontsize=12,
    )

    import matplotlib as mpl

    # build a ScalarMappable for your viridis gradient
    base_cmap = plt.get_cmap("viridis")
    newcolors = base_cmap(np.linspace(0, 0.8, 256))
    viridis_no_yellow = mcolors.ListedColormap(newcolors, name="viridis_no_yellow")
    cmap = viridis_no_yellow

    norm = mpl.colors.Normalize(vmin=0.02, vmax=0.06)
    sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)

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

    plt.tight_layout()
    output_file = r"..\images\tree2map.svg"
    plt.savefig(output_file, format="svg")
    output_file = r"..\images\tree2map.png"
    plt.savefig(output_file, format="png")

    plt.close(fig)


if __name__ == "__main__":
    main()

    # nwk_file = r"..\input\aligned_caffeine_tree.nwk"
    # gps_coords = r'..\tmp\file_w_caffeine.csv'
    # offsets_file = r"..\input\offsets_caff.csv"
# example use
# python .\tree_to_map_caffeine_content.py --nwk nwk_file --gps gps_coords --offset offsets_file
