<h1 align="center">py-madaclim: a simple Python API to interact with the Madaclim database</h1>
    <!-- badges should work when repo will go public -->
    <p align="center">
        <img alt="GitHub top language" src="https://img.shields.io/github/languages/top/tahiri-lab/py-madaclim?logoColor=blue">
        <img src="https://img.shields.io/github/contributors/tahiri-lab/py-madaclim?color=orange&logo=github"></img>
        <img src="https://img.shields.io/github/last-commit/tahiri-lab/py-madaclim?color=purple&logo=github"></img>
        <img src="https://img.shields.io/website/https/tahirinadia.github.io.svg">
    </p>
<!-- table of contents -->
<details open>
    <summary>Table of Contents</summary>
        <ol style>
            <li>
                <a href=#package-description>Package Description</a>
            </li>
            <li><a href=#installation>Installation</a></li>
                <details><summary>OS-specific steps</summary>
                <ul>
                    <li ><a href=#install-linux>Linux/UNIX-based systems</a></li>
                    <li ><a href=#install-win>Windows 10/11</a></li>
                    <li><a href=#install-mac>macOS</a></li>
                </ul>
                </details>
            <li><a href=#workflow>General workflow</a></li>
            <li><a href=#example>Getting started (quick example)</a></li>
            <li><a href=#refs>References</a></li>
            <li><a href=#contact>Contact us</a></li>
        </ol>
</details>
<!-- package description -->
<section>
    <h2 id="package-description">py-madaclim Package Description</h2>
        <img src="example/example_coll_plot.png" alt="Example of a MadaclimCollection plot with layer 79">
        <p>
            <code>py-madaclim</code> is a Python 3+ package that allows to interact with the <a href="https://madaclim.cirad.fr/">Madaclim db</a>, an open-source climate and environmental database for Madagascar.
        </p>
        <p>
            Fetch and explore the raster-based data with metadata information support, create new datasets from existent spreadsheets/csv/dataframes from any Coordinate Reference System (CRS) and explore/manipulate your data with visualization and transformation tools.
        </p>
        <p>
            Official documentation can be found here: <a href="">py-madaclim-docs</a>
        </p>
</section>

<!-- Installation -->
<section>
<h2 id="installation">Installation</h2>
<h3>Environment setup and requirements</h3>
<p>
<code>py-madaclim</code> is working with Python 3.10 and 3.11. We have provided two ways to setup a working environment for both versions:</p>
<ul>
    <li>Using <code>pip</code> and <code>venv</code> for Python=3.10</li>
    <li>Using <a href="https://conda.io">Conda</a> for Python=3.11</li>
</ul>
<p>
The requirements for the conda setup can be found in <a href="https://github.com/tahiri-lab/coffeaPhyloGeo/blob/main/conda_requirements.txt">conda_requirements.txt</a> and for the pip setup in <a href="https://github.com/tahiri-lab/coffeaPhyloGeo/blob/main/venv_requirements.txt">venv_requirements.txt</a>. OS-specific installation steps are listed below:
</p>
<!-- linux -->
<h3 id="install-linux">Debian/Linux systems</h3>
<h4>Steps for <code>pip</code> installation</h4>
<ol>
<li>Clone the repo and create a new venv</li>

```bash
git clone git@github.com:tahiri-lab/py_madaclim.git
cd py_madaclim
python -m venv ~/.pyenv/py_mada_env    #python=3.10
source ~/.pyenv/py_mada_env/bin/activate
```
<li>Activate the environment and install the requirements</li>

```bash
pip install -r venv_requirements.txt    # reqs before py_madaclim
pip install .    # to install py-madaclim
```
</ol>
<h4>Steps for <code>conda</code> installation</h4>
<ol start=0>
<li>First follow these <a href="https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html">instructions</a> to install conda on your machine</li>
<li>Clone the repo and create a conda environment</li>

```bash
git clone git@github.com:tahiri-lab/py_madaclim.git
cd py_madaclim
conda create -n py_mada_env --file conda_requirements.txt
```
<li>Activate the environment and install the requirements</li>

```bash
conda activate py_mada_env
pip install .    # using pip inside conda env
```
</ol>


<h3 id="install-win">Windows</h3>
<ol>
<li>#TODO!</li>
</ol>
<h3 id="install-mac">macOS</h3>
<ol>
<li>#TODO!</li>
</ol>
</section>

<!-- Example -->
<section>
<h2 id="example">Getting Started: Quick Example</h2>
<h3>Madaclim db metadata with the <code>info</code> module</h3>

```python
# Get available methods and properties for MadaclimLayers
>>> from py_madaclim.info import MadaclimLayers
>>> mada_info = MadaclimLayers()
>>> print(mada_info)
MadaclimLayers(
	all_layers = DataFrame(79 rows x 6 columns)
	categorical_layers = DataFrame(Layers 75, 76, 77, 78 with a total of 79 categories
	public methods -> download_data, fetch_specific_layers, get_categorical_combinations
			 get_layers_labels, select_geoclim_type_layers
)

# To access all layers as a dataframe
>>> mada_info.all_layers
geoclim_type  layer_number layer_name                       layer_description  is_categorical    units
0         clim             1      tmin1   Monthly minimum temperature - January           False  °C x 10
...

# Built-in method to download the Madaclim raster files
from pathlib import Path
cwd = Path.cwd()
mada_info.download_data(save_dir=cwd)
```
<h3>Explore the Madaclim rasters with the <code>raster_manipulation</code> module</h3>
<p>MadaclimRasters basic functionalities</p>

```python
#TODO!
```

</section>

<!-- References -->
<section>
    <h2 id="refs">References</h2>
        <ul>
            <li><a href="https://madaclim.cirad.fr">Madaclim</a> @ CIRAD</li>
            <li><a href="https://madaclim.cirad.fr">Tahiri lab</a> @ Université de Sherbrooke</li>
        </ul>
</section>

<!-- Contact -->
<section>
    <h2 id="contact">Contact Us</h2>
        <p>For any questions, feedback or to get in touch with us : <a href = "mailto: Nadia.Tahiri@USherbrooke.ca">Nadia.Tahiri@USherbrooke.ca</a></p>
        <p>For our lab's other research projects, visit our <a href="https://tahirinadia.github.io/">website</a></p>
</section>
