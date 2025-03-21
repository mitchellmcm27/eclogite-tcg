**This is the code associated with "Reactive thermodynamics of crustal eclogitization and foundering" by McMillan, M., Sim, S.J., and and Wilson, C.R.**

Citation: https://doi.org/10.1016/j.epsl.2025.119302


## Installation

Docker is the fastest way to start running the models.

First, clone this repository:

```bash
git clone https://gitlab.com/mitchellmcm27/eclogite-tcg.git
```

Within the repository, build an image from the provided Dockerfile, giving it a convenient tag **eclogite**:

```bash
cd eclogite-tcg
docker build -t eclogite -f ./docker/Dockerfile.dev
```

Then run the **eclogite** image, making sure to bind the **eclogite-tcg** directory so that any changes will be reflected on your local machine:

```bash
docker run -it --rm -v $PWD:/home/tcg/shared/eclogite-tcg eclogite
```

This will start an interactive shell in the container from which you can run the models.

Users of [VS Code](https://code.visualstudio.com/) can clone this repository and open it in Docker using the provided **.devcontainers.json** (requires the Dev Containers extension). This will also take care of mounting folders and building the reactions.

The Docker image includes prebuilt binaries for the thermodynamic database (endmembers and phases) and the reactions.
It also includes the following dependencies:
- [ThermoCodegen (TCg)](https://gitlab.com/ENKI-portal/ThermoCodegen),
- a copy of this repository at **~/shared/eclogite-tcg/**,
- *tcg_slb_database*, a convenient Python interface by Cian Wilson for the Stixrude & Lithgow-Bertelloni (2011, 2021) databases,
- Python, Julia, and R languages,
- [Perple_X](https://github.com/jadconnolly/Perple_X) (v7.0.10)

As a test, run the following command within the container:

```bash
cd models && ./profile_1d
```

The paper's results can be replicated by running the following (from the **models** directory):

```bash
./thickening_model -q
```

## Thermodynamic database

ThermoCodegen generates a custom thermodynamic database using the Stixrude and Lithgow-Bertelloni (2021) data.
The compiled database is included as **tcg_slb_database/database/tcg_slb21_database.tar.gz**.
Although the scripts, source code, and data for generating this database are provided, doing so is not necessary as long as the **.tar.gz** file is present.

## Eclogitization reactions 

Descriptions of the eclogitization reactions are included as **\*.rxml** files.
Generating and compiling the C++ code for these reactions can take some time.
If reactions are edited and need to be re-built, do so as follows:

```bash
cd tcg_slb_database
scripts/generate_reactions_eclogite -v slb21
scripts/build_reactions database/reactions/[name].rxml
```
It is highly recommended, but not necessary, to pass the path to the exact **.rxml** file that needs to be built, as shown above.

The `scripts/generate_reactions` and `scripts/generate_reactions_eclogite` files provide examples of how to generate a set of reaction descriptions (**.rxml** files) between endmembers defined in a specific thermodynamic database.

## Reactive eclogitization model

Model calculations are provided as Python scripts in the **models** directory as follows:

- **./thickening_model** runs the reactive geodynamic model of crustal thickening.
- **./pseudosection** generates a (_T_,_P_) pseudosection and plots reactive density compared with equilibrium (Perple_X) and mantle.
- **./profile_1d** generates a 1-d profile through (_T_,_P_)-space for comparing reactive phases equilibrium.
- **damkohler-fit.ipynb** shows how Damkohler number is fit to empirical data.
- **viscosity-B.ipynb** shows how viscosity coefficient _B_ is calculated for the Rayleigh-Taylor instability analysis.

By default the scripts use all available CPU cores.

Arguments can be passed as follows to customize the model runs:

| Argument    |  Purpose                           | Default                      |
|-----------------|------------------------------------|------------------------------|
|   `-n [int]`    | Number of CPU processes to use     | `mp.cpu_count()`             |
|   `-c [string]` | Bulk composition to use, by name   | hacker_2015_md_xenolith, or an array of 4 compositions   |
|   `-r [string]` | Reaction to use, by name           | eclogitization_2024_slb21_rx |
|   `-q`          | "Quick" mode (_Da_ ≤ 1e4)      | False                        |
|   `-f`          | Force model to re-calculate              | False                        |
|   `-p`          | Specific PTt path (`a`,`b`, etc.)  | all                            |
|   `-s`          | Steepness of T (`steep`,`shallow`) | None (normal)       |
|   `-t`          | Thickening curve (`tanh`,`exponential`,`decay`) | None (linear)   |

All arguments are optional.

In most cases, you will use the `-c` argument to specifiy a bulk composition, provided that Perple_X data exist for it in the **perple_x/output/** folder (see below).
The `-q` flag is useful for reducing computational time by avoiding unnecessarily large _Da_ values.

The **thickening_model** script will attempt to load the previous output from a **.pickle** file, if it exists. This is convenient for adding or modifying plots and other post-processing. To override this behavior, use the `-f` argument to force the model to re-calculate from scratch.

### Model output

All outputs are saved to the **models/output/** directory.
Outputs will be automatically grouped into subdirectories based on the name of the script, reaction, and composition.

### Model initialization using Perple_X

The models depend on equilibrium thermodynamic data for their initial compositions.
Perple_X determines the equilibrium thermodynamic properties of a composition using a Gibbs free energy minimization routine.
The model scripts described above read the text files generated by Perple_X to initialize themselves.

### Equilibrium data for compositions

Perple_X equilibrium data are already included for the bulk compositions discussed in the accompanying manuscript.
Additional compositions are available in the **perple_x/compositions.json** file, which can be edited to add more bulk compositions.

A Julia interface for programatically interacting with Perple_X (based on [StatGeochem.jl](https://osf.io/tjhmw/)) is included as **perplex_api.jl**.
The script **solve_composition.jl** reads **compositions.json** file, passes the data to Perple_X, and produces the necessary files.

To generate equilibrium data for a new composition:

- Add major element oxide percentages to **models/perple_x/compositions.json** file, making sure to use the same format as the existing entries.
- Make sure that the dictionary key for each composition object is unique, as it will be used to identify the composition throughout the code.
- Within the **perple_x** directory, run `julia solve_composition.jl [name]`, where `[name]` is the composition's dictionary key in **compositions.json** (e.g., `sammon_2021_lower_crust`).
- Because a pseudosection calculation in Perple_X (the **vertex** program) can take some time, **vertex** only runs if its output data don't already exist. If **vertex** data need to be updated, use the `-f` argument with `solve_composition.jl` to force Perple_X to re-calculate everything from scratch.
- After Perple_X completes successfully, the equilibrium data can be used in the reactive models using the `-c [name]` argument, as described above.

### Summary CSV file and plotting in R

Some plots may be more convenient to make in R after the models have been run.
For this purpose, the model saves a summary of key outputs to the **_summary.csv** file.
An R script for reading and working with this file is provided in the **./r/** directory.
