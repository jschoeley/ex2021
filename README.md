# Bounce backs amid continued losses: Life expectancy changes since COVID-19

v0.92
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6653291.svg)](https://zenodo.org/record/6653291)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6653179.svg)](https://zenodo.org/record/6653179)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6653120.svg)](https://zenodo.org/record/6653120)

## Introduction

This is a repository to accompany 'Bounce backs amid continued losses: Life expectancy changes since COVID-19', forthcoming in Nature Human Behavior. The replication files for this paper include customised functionality written in the [**R**](https://www.r-project.org/) statistical programming language.

## Prerequisites

As a pre-requisite to running this locally, you will need a working installation of [**R**](https://www.r-project.org/) with all of the necessary dependencies installed.

## Running the Code

To run this code, do something like:

```console
git clone https://github.com/jschoeley/ex2021.git
```

and then execute each of the scripts in ascending numeric order which will undertake sequential tasks like downloading and harmonizing data, running life table simulations, and outputting tables and visualizations. Output files are tagged with the same numeric code as the source file generating the file.

## Structure

- _cfg_ relates to: configuration files
- _dat_ relates to: input source data
- _out_ relates to: output data, figures, and sensitivity analysis
- _src_ relates to: code to replicate the wrangling, analysis and visualisation
- _tmp_ relates to: a subdir to store temporary files
- _ass_ relates to: a place to store repo assets

### License

This work is free. You can redistribute it and/or modify it under the terms of the GNU Public license and subject to all prior terms and licenses imposed by the free, public data sources provided by the HMD-STMF, CoverAge-DB, UK-ONS, and US-CDC (i.e. the 'data originators'). The code comes without any warranty, to the extent permitted by applicable law.
