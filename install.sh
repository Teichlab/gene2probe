NEW_ENV="gene2probe_env_tr"
conda create -y -n $NEW_ENV python=3.11

# Packages to install
conda_forge_bioconda_packages=(
    "ipython"
    "blast"
    "bedtools"
    "nb_conda_kernels"
)

pip_packages=(
    "biopython"
)

# Install packages from conda-forge or bioconda
mamba install -y -n $NEW_ENV -c conda-forge -c bioconda "${conda_forge_bioconda_packages[@]}"

# Install packages from pypi
mamba run -n $NEW_ENV pip install "${pip_packages[@]}"

#other installations
conda run -n $NEW_ENV pip install git+https://github.com/Teichlab/gene2probe.git 

conda activate $NEW_ENV

