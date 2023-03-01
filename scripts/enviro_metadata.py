import pandas as pd
from pathlib import Path

root_dir = Path(__file__).parent
enviro_meta = pd.read_html("https://madaclim.cirad.fr/environ-data-format/", index_col="Layer")[0]
# print(enviro_meta)
enviro_meta.to_csv(root_dir / "enviro_metadata.csv")