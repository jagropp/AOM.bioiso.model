# AOM.bioiso.model
Source code to calculate carbon and hydrogen isotope evolution in anaerobic oxidation of methane (AOM), depending on availbility of sulfate. This code is part of a manuscript under review (Wegener G., Gropp J., Taubner H., Halevy I. and Elvert M.), with the title "Sulfate-dependent reversibility of intracellular reactions explains the opposing isotope effects in the anaerobic oxidation of methane". In this manuscript, we describe the evolution of carbon and hydrogen isotopic composition of methane and the reversibility of the metabolic pathway as a function of sulfate concentrations. 

## Using the repository
There are two version of the model here:
1. A two-box model based on the model presented in Yoshinaga et al., 2014 (*Nat. Geo. Sci.*, doi.org/10.1038/ngeo2069). This model includes methane and DIC.
2. A four-box model that is presented in the manuscript, which includes the metabolites methane, DIC, CH3-H4MPT, CHO-MFR and HS-CoB.

There are two main functions to run the models:

### Two-box model
Run the function `run_AOM_model_twobox()`.

### Four-box model
Consists of three functions, which are initiallized by running the function `run_AOM_model()`, which will use the best fit model parameters to generate a figure for the change of isptope compositions of methane and DIC with time. 
In the function `run_bioiso_model` the user can manually change the initial model conditions (concentrations, kinetic fractionation factors from the literature, equilibrium fractionation factors).
The function `AOM_ODEs_solver` is the ODE solver for the model.

## Additional information
Please contact me for questions/help with using this code at jonathan[dot]gropp[at]weizmann[dot]ac[dot]il.
