**This is the code associated with "Reactive thermodynamics of crustal eclogitization and foundering" by McMillan, M., Sim, S.J., and and Wilson, C.R.**

Cite the paper: TBD

Cite the code:
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10835976.svg)](https://doi.org/10.5281/zenodo.10835976)

## Installation

The following downloads a prebuilt Docker image, starts an interactive container, and binds a local directory **eclogite-output** to the output folder inside the container.

```bash
docker run -it --rm -v ./eclogite-output:/home/tcg/shared/eclogite-tcg/models/output registry.gitlab.com/mitchellmcm27/eclogite-tcg
```

This Docker image includes prebuilt binaries for the thermodynamic database (endmembers and phases) and reactions.
It also includes dependencies such as
- an installation of [ThermoCodegen (TCg)](https://gitlab.com/ENKI-portal/ThermoCodegen), which is required for running the models,
- the present repository at **~/shared/eclogite-tcg/**,
- *TCg_SLB*, which provides convenient Python classes and scripts for working with the Stixrude & Lithgow-Bertelloni (2011, 2021) databases,
- Python (including numpy, matplotlib, and scipy),
- Julia, and
- the equilibrium thermodynamics software [Perple_X](https://github.com/jadconnolly/Perple_X) (v7.0.10).

As a test, run the following command within the container:

```bash
cd models && python3 parallel_profile.py
```

The paper's resluts can be replicated by running the following command and inspecting the outputs:

```bash
python3 parallel_experiment2.py -q
```

Model outputs should appear on your local machine in a newly created **eclogite-output** directory.

## Thermodynamic database

ThermoCodegen was used to generate a custom thermodynamic database using the data from Stixrude and Lithgow-Bertelloni (2021).
The compiled database is included as **tcg_slb/database/tcg_stx21_database.tar.gz**.
Although the scripts, source code, and data for generating this database are provided, doing so is not necessary as long as the **.tar.gz** file is in place.

## Reactions 

Descriptions of the eclogitization reactions are included as **\*.rxml** files.
Because generating the C++ code for these reactions can take some time, the provided Docker image includes pre-built binaries.
If reactions are edited and need to be re-built, do so as follows:

```bash
cd tcg_slb
scripts/generate_reactions_eclogite -v 21
scripts/build_reactions database/reactions/[name].rxml
```
It is recommended, but not necessary, to pass the path to the exact **.rxml** file that needs to be built, as shown above.

The `scripts/generate_reactions` and `scripts/generate_reactions_eclogite` files provide examples of how to generate a set of reaction descriptions (**.rxml** files) between endmembers defined in a specific thermodynamic database.

## Reactive eclogitization model

Model calculations are provided as Python scripts in the **models** directory as follows:

- `parallel_experiment2.py` runs the reactive geodynamic model of crustal thickening.
- `parallel_pd.py` generates a (_T_,_P_) pseudosection for comparing reactive density with equilibrium (Perple_X).
- `parallel_profile.py` generates a 1-d profile through (_T_,_P_)-space for comparing reactive phases equilibrium.
- `damkohler-fit.ipynb` shows how Damkohler number is fit to empirical data.

By default the **parallel_\*** scripts use all available CPU cores.

Arguments can be passed as follows to customize the model runs:

| CLI argument    |  Purpose                           | Default                      |
|-----------------|------------------------------------|------------------------------|
|   `-n [int]`    | Number of CPU processes to use     | `mp.cpu_count()`             |
|   `-e [float]`  | Dimensionless ending time          | 1                            |
|   `-c [string]` | Bulk composition to use, by name   | hacker_2015_md_xenolith, or an array of 4 compositions   |
|   `-r [string]` | Reaction to use, by name           | eclogitization_2024_stx21_rx |
|   `-q`          | "Quick" run (_Da_ ≤ 1e4 only)      | False                        |
|   `-f`          | Force model to re-calculate              | False                        |
|   `-p [string]` | String to prefix output directory  | None                         |

All of the arguments are optional.
In most cases, you will use the `-c` argument to specifiy a bulk composition, provided that Perple_X output data exist for it in the **perple_x/output** folder.
The only requirement to generate such data is an oxide bulk composition in mol% or wt% (see Perple_x section below).
The `-q` flag useful for reducing computational time.

### Model output

All outputs are saved to the **models/output** directory.
Outputs will be automatically grouped into subdirectories based on the name of the script, reaction, and composition.

### Model initialization using Perple_X

The models depend on equilibrium thermodynamic data for their initial compositions.
Perple_X determines the equilibrium thermodynamic properties of a composition using a Gibbs free energy minimization routine.
The model scripts described above read the text files generated by Perple_X to initialize themselves.

### Adding new compositions

Perple_X equilibrium data are already included for the bulk compositions discussed in the accompanying manuscript.
Additional compositions are available in the **perple_x/compositions.json** file, which can be edited to add more bulk compositions.

A Julia interface for programatically interacting with Perple_X (based on [StatGeochem.jl](https://osf.io/tjhmw/)) is included as **perplex_api.jl**.
The script **solve_composition.jl** reads **compositions.json** file and passes the data to Perple_X.

To generate equilibrium data for a new composition:

- Add the oxide percentages to **models/perple_x/compositions.json** file, making sure to use the same format as the existing entries.
- The dictionary key for each composition object must be unique, as it will be used to refer to the composition throughout the codebase.
- Within the **perple_x** directory, run `julia solve_composition.jl [name]`, where `[name]` is the composition's dictionary key in **compositions.json** (e.g., `sammon_2021_lower_crust`).
- Because a pseudosection calculation in Perple_X (the **vertex** program) can take some time, **vertex** only runs if its output data don't already exist. If **vertex** data need to be updated, use the `-f` argument with `solve_composition.jl` to force Perple_X to re-calculate everything from scratch.
- After Perple_X completes successfully, the equilibrium data can be used in the reactive models using the `-c [name]` argument, as described above.

### Additional plotting in R

Some plots are more convenient to make in R after the models have been run.
For this purpose, the model saves a summary of key outputs to a **_summary.csv** file, and an R script for reading and working with this file is provided in the **./r** directory.
