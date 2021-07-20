README for code associated with paper "Modeling multi-scale occupancy for monitoring rare and highly mobile species"


FOLDERS:

ct_vs_dt_sims: runs analysis for Methods subsection "Simulation"
	Contains:
	runallscens_idds.R: wrapper script with functions and JAGS code
        runscen_idds.R: simulates data and runs models for a specific simulation scenario

data_analysis: see "OTHER FILES"

power_analysis: runs analysis for Methods subsection "Power analysis"
	Contains:
	pop_sim_and_detect_functions.R: functions to simulate wolverine populations, movement, and detection by camera grids
	runallscens_idds.R: wrapper script with simulation setup and JAGS code
        runscen_idds.R: simulates data and runs models for a specific simulation scenario

OTHER FILES:

The data used in this paper are not owned by any of the coauthors on this manuscript. I have provided code
that I used to run the data analysis and analyze simulation results in highlymobile_occ.Rmd.
The data analysis code is provided in the code chunk titled "wolverine_data_analysis".
JAGS model files for the data analysis are provided in wolverinepoiswa.txt and wolverinedtdowa.txt
Please contact R.L. Emmet for questions about the data analysis code.


INSTRUCTIONS TO RUN CODE IN FOLDERS:

1. Begin a new R session/ensure global environment is empty.
2. Install any packages required to run the code (including external software JAGS).
3. Change line 4 in runallscens_idds.R and line 3 in runallscens_pwr.R.
4. Run code.