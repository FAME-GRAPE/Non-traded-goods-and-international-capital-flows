===========================================================================================

Title:   	"Non-Traded Goods, Factor Markets Frictions, and International Capital Flows"
Authors: 	Jacek Rothert; Jacob Short
Ms. No.: 	RED-D-20-00201
Editor:  	Jonathan Heathcote

===========================================================================================


Software and OS: MATLAB 2021b, Windows 

Run in the following order:
 
1) JRJS_NT_main.m FOUR times
1a) use current, default numbers
1b) line 71, set deep_params.teta = 0.5 
1c) line 71, set deep_params.teta = 2
1d) line 73, set deep_params.omc = 0.5; line 71, set deep_params.teta = 0.1

2) JRJS_NT_counterfactuals.m

3) dD_over_dtaus.m

Expected computation time (see notes below)






DATA INPUTS
===========================================================================================
JRJS_data.xlsx
% Excel file with cross-country statistics - input into JRJS_NT_main.m 
===========================================================================================



MODEL CALIBRATION
===========================================================================================
RUN: JRJS_NT_main.m - main MATLAB script to calibrate the model for all countries; 

DEPENDENCIES:
a) initial_calibration.m
	(a-i) calls eqm_stst.m
b) eqm_path.m
c) main_calibration.m
	(c-i) calls eqm_path.m

INPUTS: 
- JRJS_updateddata.xlsx

OUTPUTS (only one output, for each run of the file):
- Calibration_Results_teta01_ppsi1000.mat (set line 47, 48, and 49 to 0.1, 1000, and 0.2, respectively)
- Calibration_Results_teta05_ppsi1000.mat (set line 47, 48, and 49 to 0.5, 1000, and 0.2, respectively)
- Calibration_Results_teta2_ppsi1000.mat (set line 47, 48, and 49 to 2, 1000, and 0.2, respectively)
- Calibration_Results_teta01_ppsi1000_omc05.mat (set line 47, 48, and 49 to 0.1, 1000, and 0.5, respectively)


EXPECTED RUN TIME: 13-15 hours each (parallel computing, with 32 nodes)
===========================================================================================




COUNTERFACTUALS
===========================================================================================
RUN: JRJS_NT_counterfactuals.m - main MATLAB script to run counterfactuals 
(after model calibration)

DEPENDENCIES
a) calibration_1sector.m
	(a-i) calls eqm_path_1sector.m
b) counterfactual.m
	(b-i) calls eqm_path.m
	(b-ii) calls eqm_path_flex_labor.m

INPUTS: 
- Calibration_Results_teta01_ppsi1000.mat 
- Calibration_Results_teta05_ppsi1000.mat 
- Calibration_Results_teta2_ppsi1000.mat 
- Calibration_Results_teta01_ppsi1000_omc05.mat 

OUTPUTS: 
- OneSec.mat
- results_counterfactuals.mat
- prints parts of Table 4 and Table 5 in the MATLAB Command Window

EXPECTED RUN TIME: 4-5 hours (parallel computing, with 32 nodes)
===========================================================================================




SENSITIVITY
===========================================================================================
RUN: dD_over_dtaus

INPUTS: 
- OneSec.mat (output of JRJS_NT_counterfactuals.m)
- results_counterfactuals.mat (output of JRJS_NT_counterfactuals.m)


OUTPUTS: 
- prints Table 6 in the MATLAB Command Window

EXPECTED RUN TIME: 1-3 seconds
===========================================================================================


