# sturgeoninterval

Data and code to run Bayesian multi-state capture-mark-recapture models for the manuscript: 

Jacobs, G.R., D. Gorsky, Z. Biesinger, M.A.H. Webb, S.J. Wenger, C.W. Osenberg. Warmer springs may lead to more frequent spawning in lake sturgeon. 

## Folder contents

Data

* `int.data_240826.RData` - File containing data objects conditioned for model input via R

Scripts

* `IMSCJS_m9.0.3.txt` - JAGS program including interactive *G\*X* effects on transition probabilities
* `IMSCJS_m9.0.3_additive.txt` - JAGS program including additive *G+X* effects on transition probabilities
* `IMSCJS_m9.0.3_null.txt` - JAGS program omitting *X* effects from transition probabilities
* `IMSCJS_m9.0.3_RUNALL.r` - R script using conditioned data and JAGS programs to fit 11 candidate models compared in the manuscript

## Workflow

1) run `IMSCJS_m9.0.3_RUNALL.r` (note: run time was 33 hours on a machine with an i9-13900 processor and 64 GB RAM)
 
	- Imports conditioned data and runs JAGS programs to fit 11 candidate models compared in the manuscript
	- Produces 11 folders containing JAGS model outputs, with the naming convention "runjagsfiles9.0.3_VARIABLE_TYPE"
		- "VARIABLE" is one of the X temperature variables described in the manuscript. 
		- "TYPE" is either "interaction" for models run with the `IMSCJS_m9.0.3.txt` script, "additive" for models run with the `IMSCJS_m9.0.3_additive.txt` script, or "null" for the model run with the `IMSCJS_m9.0.3_null.txt` script. 
		- Example: the folder `runjagsfiles9.0.3_Dgt5lag1_additive` contains the results of the JAGS model evaluating the additive effect of *Sex* (*G*) and the temperature (*X*) variable *Dgt5lag1*. 

2) Summarize, visualize, and compare JAGS models using your favorite statistical program(s). 
