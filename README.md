This repository contains the code and data used in [Williamson *et al.* (2025)]():

"On real-time calibrated prediction for complex model-based decision support in pandemics: Part I"

Daniel B. Williamson, Trevelyan J. McKinley, Xiaoyu Xiong, James M. Salter, Robert Challen, Leon Danon, Ben Youngman and Doug McNeall

All simulator model runs were conducted on a MacBook Pro with an M1 Pro chip and 16GB of RAM.

## File descriptions

The repository contains the following folders:

* `perfect_model_experiments/DeathOnly/`: folder containing code for **Section 4.2.1 - Assimilating death counts only**.
* `perfect_model_experiments/Death_HospitalCases/`: folder containing code for **Section 4.2.1 - Assimilating death counts and hospital cases**.
* `Application/`: folder containing code for **Section 5 - Application: UK COVID-19 lockdown, March 2020**.

A brief description of the files in each folder is provided below, with further details in later sections.

* `data/`: folder containing real data or simulated data.
* `inputs/`: folder containing additional data and scripts necessary for running the 
models.
* `BPF.R`: R functions to run the particle filters. Arguments to the main `BPF` function
are defined at the top of the file.
* `BPF_foreCasts.R`:  R functions to run the particle filters. This variant differs from `BPF.R` by preserving all particle trajectories (`saveAll = TRUE`) and incorporating sampled observation error into these trajectories to facilitate comparison with observed data.
* `stochModel.R`: Code to simulate data. This file is only used for the perfect model experiments.
* `checkWavexDesign.R`: code to check the design points at each new wave and convert the 
 design into the correct format for use in the model.
* `discreteStochModel.cpp`: C++ implementation file containing a discrete-time stochastic model for COVID-19 simulation that can be called from R.
* `TruncSkellams.cpp`: C++ implementation file containing functions for computing log probabilities and sampling from truncated Skellam distributions.
* `tnorm.cpp`: C++ implementation for fast truncated Gaussian sampling that can be called from R.
* `bessel.h`: C header file that defines constants and parameters for calculating Bessel functions in the R programming language.
* `wave1Design.R`: code to generate a Wave 1 design across the input space.
* `ensembleForecasts.R`: code to generate ensemble forecasts for any specified wave.
* `wavexRuns.R`: code to run the model across multiple design points for a specified wave.
* `emulateDGP_wave1.R`: implementation of history matching (HM) with the Deep Gaussian Process (DGP) emulator configured for Wave 1 design points.
* `emulateDGP_wavex.R`: implementation of HM with the DGP emulator configured for arbitrary Wave x design points, enabling multi-wave analysis.
* `emulator_fns.R`: custom functions supporting the HM framework with a DGP emulator.


## Workflow

The basic workflow is described below, but please read all the corresponding sections for details on each of these steps.

1. Run the `dataProcess.R` script in the `data/` folder to clean
 up the raw data and extract it into a usable form for the model.
2. Generate an initial Wave 1 design using `wave1Design.R`.
3. If conducting a simulation study, then use `stochModel.R` to generate simulated data.
4. Run the ensemble on a laptop machine. This uses `wavexRuns.R` and other associated files
 to run the model for each design point. This give the model outputs along with some summary plots of the ensemble trajectories.
5. Perform history matching by running `emulateDGP_wave1.R` for Wave 1 and `emulateDGP_wavex.R` for subsequent waves. 
6. Check the new design and set up the files in the correct format for the next wave
 using `checkWavexDesign.R`.
7. Return to step 4 and repeat until sufficient waves have been run.

## Death and hospitalisation data

The `data/` and `inputs/` folders contain the raw data files, and some scripts for tidying them up and extracting relevant time periods over which to fit the models. We have put the original links where the data were downloaded below, but note that some of these are not accessible anymore. Nevertheless, the resulting downloaded data are included.

The raw data consist of:

* Deaths reported within 28 days of a positive test by age and region: `data/nation_2021-05-10.csv`, downloaded from
[https://api.coronavirus.data.gov.uk/v2/data?areaType=nation&areaCode=E92000001&metric=newDeaths28DaysByDeathDateAgeDemographics&format=csv&release=2021-05-10](https://api.coronavirus.data.gov.uk/v2/data?areaType=nation&areaCode=E92000001&metric=newDeaths28DaysByDeathDateAgeDemographics&format=csv&release=2021-05-10) [link defunct]. As the record begins on 2nd March 2020, we assume that deaths from COVID-19 are zero between 15th January–2nd March.

Change the working directory of R to `data` and run:

```
source("dataProcess.R")
```
This will generate the file 'deathTotalByAge.rds'.

**Note**: daily hospitalisation data from January 15th to March 23rd, 2020 is not publicly available and was originally obtained from the COVID-19 Hospitalisation in England Surveillance System (CHESS). For demonstration purposes in this repository, we have generated random numbers to represent this period, saved as `'data/hospitalCasesTotal_Eng_forDemo.rds'`. Both this file and the `'deathTotalByAge.rds' file` are used in the `'wavexRuns.R'` and `'ensembleForecasts.R'` scripts in the `'Application'` folder to assimilate death counts and hospital cases data up to the first lockdown period.


## Other data files, scripts and objects necessary for sampling from the input spaces

The `inputs/` folder also contains additional files required to run the model. These are:

* `inputs/age_seeds.csv`: this contains the proportion of the UK population in each of the age-classes of the model, derived from the Office for National Statistics. (2020). 2011 Census: Aggregate Data. [data collection]. UK Data Service. SN: 7427, DOI: [http://doi.org/10.5257/census/aggregate-2011-2](http://doi.org/10.5257/census/aggregate-2011-2).
* `inputs/coMix_matrix.csv`: this is a contact matrix between different age-classes used *after* the first lockdown.
The original data (`20200327_comix_social_contacts.xlsx`) are from [Jarvis *et al.* (2020)](https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-020-01597-8#availability-of-data-and-materials) and
were downloaded from here 

[https://cmmid.github.io/topics/covid19/comix-impact-of-physical-distance-measures-on-transmission-in-the-UK.html](https://cmmid.github.io/topics/covid19/comix-impact-of-physical-distance-measures-on-transmission-in-the-UK.html)

and we used the `All_contacts_imputed` sheet with column `A` and row `1` removed.
* The `inputs/dataTools.R` file contains helper function required to run the model.
* `inputs/hospStays.rds`: this is a fitted finite mixture model (FMM) used to sample from the input
space for the parameters guiding the length of hospital stays.
* `inputs/hospThresh.rds`: is a threshold on the p.d.f. of the FMM in `inputs/hospStays.rds`, which
is used to define a boundary for the input space.
* `inputs/parRanges.rds`: contains the input ranges for the parameters with a hypercube input space.
* `inputs/pathways.rds`: this is a fitted FMM used to sample from the input space for the parameters 
guiding the probabilities of transitioning down different epidemiological pathways.
* `inputs/pathThresh.rds`: is a threshold on the p.d.f. of the FMM in `inputs/pathways.rds`, which
is used to define a boundary for the input space.
* `inputs/POLYMOD_matrix.csv`: this is a contact matrix between different age-classes used *before* 
the first lockdown. This was extracted from the [`socialmixr`](https://cran.r-project.org/web/packages/socialmixr/vignettes/socialmixr.html) package, which uses data in [Mossong *et al.* (2008)](https://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.0050074).

## Wave 1 design

The `wave1Design.R` file contains code to generate a Wave 1 space-filling design across the original input space. This uses a Latin Hypercube design for some parameters, and then a space-filling design for the parameters with non-rectangular input spaces, as described in [Williamson *et al.* (2025)](). This script takes a single argument which corresponds to a name for the outputs (which is appended to the word "wave" to store the model runs in a standardised way across the waves). For example, to produce a Wave 1 design in a folder called `wave1`, you can run the following from a terminal window:

```
R CMD BATCH --no-save --no-restore --slave '--args 1' wave1Design.R
```

Alternatively, the script can be run from within R in the usual way, noting that you will have to change the default wave name manually. 

The contents of the `wave1` folder are:

* `design.pdf`: a plot of the Wave 1 design.
* `disease.rds`: a `tibble` object containing the design points but converted to the necessary
form for use in the model.
* `inputs.rds`: a `tibble` object containing the design points.

## Simulation study

The code used to simulate a synthetic data set from the underlying model can be found in the `stochModel.R` file. This also requires a set of input parameters to use for the simulation. We did this by selecting a set from an initial set of design points generated using the `wave1Design.R` file. 

To generate simulated data, first run the `wave1Design.R`  file to create a Wave 1 Latin Hypercube Sampling (LHS) design. Then execute `stochModel.R`  to simulate data using one of these design points. On line 31, you can specify which parameter set to use for simulation. In the example, we use the tenth set (`10`) as the 'true' parameters and simulate a 90-day period. The script automatically creates an output folder named `data`, containing the simulated results. 

Change the working directory of R to `perfect_model_experiments/DeathOnly` and run:

```
source("stochModel.R")
```
This runs a small number of replicates, and then picks the run that is the closest in L2 distance to the median across the replicates as a representative model run. 

The contents of the output folder "data" contains files:

* `disSims.rds`: a `tibble` object storing the simulated hidden states of the system
* `pars.rds`: a `tibble` object containing the parameters used to simulate the data.
* `sims.pdf`: a plot of the simulated data and  the hidden states.

Copy the generated 'data' folder to `perfect_model_experiments/Death_HospitalCases` to assimilate death counts and hospital cases for the perfect model experiment.

## Running the model

The `wavexRuns.R` file contains code for running an ensemble of design points for a given wave. It takes one argument: the "wave" folder where both design points and data are stored. For example, to run the first design point of the wave 1 design (stored in the `wave1` folder), using the data (log-likelihood estimates stored in files named `runs_md.rds`) also located in the `wave1` folder, we can execute:


```
R CMD BATCH --no-save --no-restore --slave '--args 1 ' wavexRuns.R
```

## History matching

The files `emulateDGP_wave1.R`,`emulateDGP_wavex.R``emulator_fns.R` contain code for building and validating emulators, and for generating new designs. 

**Note**: the Wave 1 code works slightly differently from that of the subsequent waves. 

### Wave 1

For the Wave 1 emulators, check the prior ranges on L16--30.

To select the active variables for the deep Gaussian Process (DGP) emulator, we first perform a sensitivity analysis and then choose the first 10 variables with the strongest main effects as the active variables (L85-117).

Once this has been done, and the active variables selected, then we can fit the DGP to the active variables using MCMC, which is done on L121--142. As with any MCMC method it is useful to assess convergence, and so if the trace plots do not seem to have converged then the chain can be run for longer. If you need to do this then the code on L154--165 can be uncommented and repeated as required.

**Note**: for many of the runs we conducted, the mixing of the chains aren't as good as we would ideally like. However, the MCMC here is used mainly to find the area of high posterior mass, and the predictions from the DGP actually only a point prediction corresponding to the median value from the last half of the chain. To improve mixing the chain could be run for longer and then thinned, but this would take a long time for large-scale problems, and so instead here we take a more pragmatic approach and run it until it looks like it has converged at least in the last half of the chain. We then produce validation plots, and if the validation looks OK we proceed with the rest of the process. We found in tests that running the chains long enough to get good mixing did not noticably improve the validation, as long as the median value was roughly correct.

The code then stores any "doubt points" on L220--254 (as described in the paper), and then builds a custom emulator object on L256--284. The `hmer` package allows us to build custom emulators, using a function called `Proto_emulator()`, but since the whole history matching approach is different to the standard approach, we wrap all of the custom implausibility measures, target log-likelihoods, and necessary files for sampling from the non-hypercube input spaces into this object. The `emulator_fns.R` file contains custom functions to do all of these steps, which are used in the construction of the `emulator` object.

Then we run a final check to ensure that the current best point is retained in the NROY space (L286--294). A new design is then generated by initially sampling a large number of points from the original input space (1,000 in this case; L311--367), and then extracing the subset of these that are retained in the NROY space. The space removed at this wave is estimated on L370--373 (noting that the "target" log-likelihood is stored in the `emulator` object, and as such we have to set a dummy `targets` object with `val = 1` in order for the code to work. This is because, as described above, we have many context-specific things used here which are different to the usual HM approach (such as non-hypercube input spaces, use of Voronoi regions to retain space around "doubt points", custom implausibility measures etc.), and as such we have to employ a few tricks to allow `hmer` to handle all of these aspects correctly. L381--382 then extract the subset of the baseline points that remain in the NROY space, which we use as a basis for building a new set of 1,000 points in the NROY space by utilising the slice sampling method implemented in `hmer` (L381--382). The Wave 2 training and validation points are then sub-sampled from this set of baseline points using a maximin design (L385--394). If you want to change the number of training and validation points, then you can do this on L385 and L388 respectively (noting that if you want more than 1,000 design points, then you will have to generate a new set of baseline points accordingly).

Finally, the new design is plotted (L397--421), and then to aid subsequent waves we then augment the set of baseline NROY points back up to 1,000 to replace the points removed for the new design (L428--429). The `inputsWave2.csv` then contains the training and validation points for Wave 2, and so this can be copied into a new folder with the same naming convention as used in Wave 1 (so if the Wave 1 runs are stored in `wave1`, then create a new folder called `wave2` and copy the `inputsWave2.csv` file into this folder before running the new wave).

### Wave 2 onwards

For the Wave 2+ emulators, `emulateDGP_wavex.R` contains the same basic structure as the Wave 1 code described above, but with a few key differences. As in Wave 1, the code requires some manual checking as it is running, and so it is best to run this code
interactively in R.

The first thing to note is that since a consistent naming convention is used for the folders containing the outputs of each wave, then the `emulateDGP_wavex.R` code can automatically link to the outputs of the previous waves where required (and can copy across seeds, prior ranges etc.). It takes one argument: the current "wave" folder where outputs of previous waves (copied across) and current wave are stored. For example, for Wave 2, we can execute:

```
R CMD BATCH --no-save --no-restore --slave '--args 2' emulateDGP_wavex.R
```

Then the code re-builds the emulators for the previous waves, which are required to sample new points and check the NROY space etc. Note that in standard `hmer` code, we could save the previous wave emulators and simply load them back in again as a single object. The complexity here is that the emulators are built using the `dgpsi` package, which uses Python, rather than R, under-the-hood. As such we can save the Python objects as `.pkl` files, but we cannot include the Python objects embedded in an R object saved as an `.rds` file. We overcome this by saving all of the relevant components at each wave, and then looping over the previous waves to rebuild the emulators. This is done on L47--84.

Then, to ensure maximising the information in the training data, we load in all the previous wave design points, and then extract any of these that remain in the previous wave NROY space. We can then augment the training data with these additional points. This whole process is done on L87--149. As before, check the number of training and validation points on L120 and L125 respectively (these may have to be amended if you
removed any failed runs when running the ensemble).

Then the code runs similarly to Wave 1, where we define active variables (L170--201), build the DGP emulator
using the active variables (L210--289), check and record "doubt points" (L301--334), set up the new custom emulator object for use in `hmer` (L337--371) and check that the best point is retained in the previous NROY space (L374--387).

Then we load in the set of baseline NROY points from the previous wave, and use these to estimate the space removed at the current wave (L390--396). Then, as before we extract the subset of the previous wave baseline points that remain in the updated NROY space, and then use these to generate 1,000 baseline points in the updated NROY space (L400--405). These are then sub-sampled using a maximin design to generate a new set of training and validation points (L408--417). Again, the number of training and validation points can be amended on L408 and L411 respectively if required. The new design points are then plotted (L423--443), before the updated set of baseline NROY points are augmented back up to 1,000 to replace the points removed for the new design (L447--448).

The resulting `inputsWaveX.csv` file then contains the training and validation points for the new Wave `X`,
and so this can be copied into a new folder with the same naming convention as used in the previous waves (so if we are at the Wave 3, and the runs are stored in `wave1`, `wave2` and `wave3`, then create a new folder called `wave4` and copy the `inputsWave4.csv` file into this folder before running the new wave).

## Forecasts

The `ensembleForecasts.R` file provides code for running forecasts for an ensemble of a particular wave and plotting the trajectories. You can set the wave number by changing line 33. Note that `ensembleForecasts.R` uses `BPF_forecasts.R`, which differs from `BPF.R` in that it saves all particle trajectories by setting `saveAll = TRUE` and adds sampled observation error to the trajectories in order to plot them against the observations. You can change line 35 to specify when you would like the forecasting to stop.


