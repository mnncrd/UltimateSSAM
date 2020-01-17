# UltimateSSAM

UtimateSSAM (Ultimate Secondary Structure Assignment Method) is a secondary structure assignment tool aiming at combining and centralising the preexisting tools to implement a more generic approach.

Four modes are available: `dssp`, `dsspcompare`, `ssam`, and `ssamcompare`.

The `dssp` mode performs a DSSP-like method, while the `dsspcompare` mode performs a DSSP-like method **AND** compares the result to DSSP.

As for now, `ssam` is equivalent to `dssp`, and `ssamcompare` is equivalent to `dsspcompare`, but they will be updated to include more secondary structure assignment methods, truly making them the Ultimate Secondary Structure Assignment Method.

## Requirements

This program is coded in Python and only uses standard Python libraries. You will need:

* Pyhton &ge; 3.7
* DSSP 3.0 

A Conda environment is available for a quick installation of all requirements.

## Installation

### Conda environment

To install the conda environment, run:

```
conda env create --file environment.yml
```
### 

To install UltimateSSAM
```
git clone https://github.com/mnncrd/UltimateSSAM.git
```

## Examples

First, load the Conda environment

```
conda activate UltimateSSAM
```

### DSSP mode

To run the `dssp` mode, run the following command:

```
python3 ssam.py dssp -i input -o output
```

With:
* `input`: a `.pdb` or `.cif` file
* `output`: a `.dssp` file

For hydrogen atoms *à la* DSSP, run: 

```
python3 ssam.py dssp -i input -o output -hy
```

### DSSPCOMPARE mode

To run the `dsspcompare` mode, run the following command:

```
python3 ssam.py dsspcompare -i input -o output -oc output-compare
```

With:
* `input`: a `.pdb` or `.cif` file
* `output`: a `.dssp` file
* `output-compare`: a `.csv` file

## Results

### DSSP mode

A `output.dssp` will be produced.

### DSSPCOMPARE mode

A `output.dssp` and a `output-compare.csv` will be produced.

## Author

Manon Curaudeau
