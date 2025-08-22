# ew-workflows repository

This repository includes functions for running the geochemical reactive transport code [SCEPTER](https://github.com/cdr-laboratory/SCEPTER) and processing its results. It has been tested with the version of SCEPTER that we cloned, here [SCEPTER](https://github.com/carbonplan/SCEPTER).

The latest release comes alongside our analysis for "Swapping carbonate for silicate in agricultural enhanced rock weathering" by Tyler Kukla, Yoshiki Kanzaki, Freya Chay, Noah J. Planavsky, and Christopher T. Reinhard.

**Learn more about SCEPTER and its strengths and limitations in the following model description papers:** 
Kanzaki, Y., Zhang, S., Planavsky, N. J., & Reinhard, C. T. (2022). Soil Cycles of Elements simulator for Predicting TERrestrial regulation of greenhouse gases: SCEPTER v0. 9. Geoscientific Model Development, 15(12), 4959-4990.

Kanzaki, Y., Chiaravalloti, I., Zhang, S., Planavsky, N. J., & Reinhard, C. T. (2024). In silico calculation of soil pH by SCEPTER v1. 0. Geoscientific Model Development, 17(10), 4515-4532


## overview
The repository relies on Argo-workflows to run batches of SCEPTER simulations. Quickstart instructions for running SCEPTER with or without Argo can be found in `scepter-quickstart.md`. 

## contents
```
├── inputs/scepter/       # input files for SCEPTER simulations
|    ├── clim/            # climate inputs 
|    ├── dust/            # dust flux inputs for variable dust application runs
|    └── params/          # parameter input .yaml for argo workflows
|
├── run_scepter/                      # functions for SCEPTER simulations and postprocessing
|    ├── argo_helper_fxns.py          # functions to support argo simulations
|    ├── batch_helperFxns.py          # functions to support batch simulations
|    ├── build_composite_multiyear.py # functions to build output dir from runs chained together
|    ├── cdr_fxns_flx.py              # functions for calculating CDR terms fluxes
|    ├── cdr_fxns_postproc.py         # functions for calculating CDR relative to baseline
|    ├── cflx_proc.py                 # functions for calculating CDR-relevant fluxes
|    ├── sa_postproc_fxns.py          # functions to process sensitivity analysis results
|    └── scepter_helperFxns.py        # functions to support SCEPTER simulations
|
├── scripts/
|    ├── figures   # scripts to make figures
|    ├── process   # scripts to process data
|    ├── run       # scripts to run the model
|    └── setup     # scripts to setup batch and spinup runs
|
├── environment.yaml
├── LICENSE
├── pyproject.toml
├── README.md              # this file
└── scepter-quickstart.md  # instructions for running SCEPTER using this repository
```

## license
This work is licensed under the MIT license. 

## about us
CarbonPlan is a nonprofit organization that uses data and science for climate action. We aim to improve the transparency and scientific integrity of climate solutions with open data and tools. Find out more at [carbonplan.org](https://carbonplan.org/).
