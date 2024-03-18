Cite this code using 10.5281/zenodo.10835976

## Clone the repo

```bash
git clone https://gitlab.com/mitchellmcm27/eclogite-tcg.git
cd eclogite-tcg
```

This will also pull in the "TCG_SLB" code, which provides convenient Python classes and scripts for working with the Stixrude & Lithgow-Bertelloni (2011, 2021) databases through TCG, as a git subtree.

## Start a Docker container

Build an image using the provided Dockerfile, using `-t` to give the image a human-readable tag.
Start an interactive container, binding the local **eclogite-tcg** directory to **workspaces/eclogite-tcg** inside the container.

```bash
docker build -t eclogite-tcg .
docker run -it --rm -v $PWD:/workspaces/eclogite-tcg eclogite-tcg
```
Alternatively, open the **eclogite-tcg** directory in VSCode and use the Dev Containers extension to automatically build and re-open the repository inside a Docker container.

## Thermodynamic database and reactions

A custom thermodynamic database based on Stixrude and Lithgow-Bertelloni (2021) is included as **tcg_slb/database/tcg_stx21_database.tar.gz**.
The scripts, source code, and data for generating this database are also provided, but doing so is not necessary.

Reactions are included as **\*.rxml** files.
Because generating the C++ code for these reactions can take several hours, pre-built binaries and Python bindings are included, which are compatible with the Docker image.

If reactions are edited and need to be re-built, you can do so as follows:

```bash
cd tcg_slb
scripts/generate_reactions_eclogite -v 21
scripts/build-reactions database/reactions/eclogitization_2024_stx21_rx.rxml
```

## Model calculations

Model calculation scripts are in the **models** directory.

```bash
cd models
```

Four python scripts are given as follows:

- `python3 parallel_pd.py` generates a (_T_,_P_) pseudosection for comparing density with Perple_X results.
- `python3 parallel_profile.py` generates a 1-d profile through (_T_,_P_)-space for comparing phase mode with Perple_X results.
- `python3 parallel_experiment2.py` runs the geodynamic model of crustal thickening at the Moho.
- `damkohler-fit.ipynb` shows how Damkohler number is fit to empirical data.

In most cases, you can pass the name of any pre-defined composition that exists in the **models/perple_x/compositions.json** file. 
For example, `python3 parallel_pd.py -c hacker_2015_md_xenolith`.
By default the **parallel_\*** scripts use all available CPU cores.

Arguments can be passed as follows to customize the model runs:

| CLI argument    |  Purpose                           | Default |
|-----------------|------------------------------------|---------|
|   `-n [int]`    | Number of CPU processes to use     | `mp.cpu_count()` |
|   `-e [float]`  | Scaled ending time                 |  1               |
|   `-c [string]` | Bulk composition to use, by name   | an array of 4 compositions |
|   `-r [string]` | Reaction to use, by name           | eclogitization_2024_stx21_rx |
|   `-q`          | Run in "quick mode" (_Da_ ≤ 1e4)   | False |
|   `-f`          | Force model to recalculate results | False |
|   `-p [string]`  | String to prefix output directory  | None |

## Model outputs

All outputs are saved to the **models/output** directory.
Outputs will be automatically grouped into subdirectories based on the name of the reaction and composition.

## Perple_X: Adding new compositions

Perple_X outputs for equilibrium thermodynamic properties are included for several compositions defined in the **models/perple_x** directory.
The python scripts read the text files output by Perple_X to initialize the reactive thermodynamic models.

If a composition needs to be added that does not already exist, use the following recipe:

- Add the oxide percentages to **models/perple_x/compositions.json** file, making sure to use the same template as the existing compositions.
- Within the **perple_x** directory, run `julia solve_composition.jl [name]` where `[name]` is the unique identifier, the key for the composition object in **compositions.json** (e.g., `sammon_2021_lower_crust`).
- The composition can now be used in python scripts by passing the `-c [name]` argument as described above.

Note that the Docker image already has a installation of Perple_X.

## Additional analyses in R

Some analyses are more convenient to perform in R after the models have been run.
For this purpose, the model saves a summary of key outputs to a **_summary.csv** file.

An R script for reading and working with this file is provided in the **./r** directory.
