############
Installation
############

**********************
Software depenedencies
**********************

Our pipeline uses a number of software packages. To facilitate their installation you can use conda package manager.

Automated installation with conda
=================================

After downloading the proper conda version from `https://docs.conda.io/en/latest/miniconda.html`, you can install it by typing in the directory with the installation script:

.. prompt:: bash $

  bash ./Miniconda3-latest-Linux-x86_64.sh

Then, you need to install git and Snakemake:

.. prompt:: bash $

  conda install git snakemake-minimal

With git, you can download the GitHub repository with our pipeline (to the directory of your choice):

.. prompt:: bash $

  mkdir {taqlore_dir}
  cd {taqlore_dir}
  git pull https://github.com/twrzes/TAQLoRe.git

Then, you can run the Snakemake file which should install all the dependencies. Alternatively, especially if you use HPC to run our pipeline (which is preferred method of running it) and you do not have access to Internet on all nodes, you can pre-install the environment by using `conda env create`:

.. prompt:: bash $

  conda env create --file {taqlore_dir}/TAQLoRe/envs/taqlore.yaml

Manual installation
===================

If you prefer not to use conda, you can install all the dependencies manually. Software versions in brackets denote ones that we used to develop our pipeline so they should work without any problems. Software dependencies include:

- Python (3.5.1)
- Perl (5.22.1)
- R (3.5.1)
- Bedtools (2.26.0)
- LAST (979)
- GMAP (20190315)
- Git (2.21)
- Snakemake (5.4.5)

Also, the workflow requires the following Python modules:

- pandas (0.24.2)
- numpy (1.16.3)
- tqdm (4.32.1)
- pysam (0.15.2)
- pybedtools (0.8.0)
- scipy (1.2.1)
- natsort (6.0.0)

It also requires following R libraries:

- heatmap3
- ggfortify
- ggplot2
- FactoMineR
- factoextra
- sva

After installing all the dependencies simply clone the GitHub repository (after creating a directory when you want to store the instance of the pipeline):

.. prompt:: bash $

  mkdir {taqlore_dir}
  cd {taqlore_dir}
  git pull https://github.com/twrzes/TAQLoRe.git
