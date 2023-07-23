Title:   	"Non-Traded Goods, Factor Markets Frictions, and International Capital Flows"

Authors: 	Jacek Rothert; Jacob Short

Ms. No.: 	RED-D-20-00201

Editor:  	Jonathan Heathcote

Research funded by National Science Center (Narodowe Centrum Nauki), grant # 2019/35/B/HS4/00769

----


Software and OS: MATLAB 2021b, Stata, Windows 

Run in the following order:
 
1) Matlab files in folder "model_codes"
	see README file contained in folder

2) Matlab and Stata files in folder "paper_figures_table"
	see README file contained in folder

3) There is an Excel file JRJS_data.xlsx with cross-country statistics - input into JRJS_NT_main.m (containted in "input_data"). The description of data construction is included in "input_data" folder.


MODEL CALIBRATION
----
RUN: JRJS_NT_main.m - main MATLAB script to calibrate the model for all countries; 

DEPENDENCIES:
a) initial_calibration.m
	(a-i) calls eqm_stst.m
b) eqm_path.m
c) main_calibration.m
	(c-i) calls eqm_path.m

INPUTS: JRJS_updateddata.xlsx

OUTPUTS (only one output, for each run of the file):
- Calibration_Results_teta01_ppsi1000.mat (set line 47, 48, and 49 to 0.1, 1000, and 0.2, respectively)
- Calibration_Results_teta05_ppsi1000.mat (set line 47, 48, and 49 to 0.5, 1000, and 0.2, respectively)
- Calibration_Results_teta2_ppsi1000.mat (set line 47, 48, and 49 to 2, 1000, and 0.2, respectively)
- Calibration_Results_teta01_ppsi1000_omc05.mat (set line 47, 48, and 49 to 0.1, 1000, and 0.5, respectively)


EXPECTED RUN TIME: 13-15 hours each (parallel computing, with 32 nodes)



COUNTERFACTUALS
----
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


FIGURES AND TABLES
----
Associated files are run after model calibration and counterfactuals are completed.

The following files will generate the figures and tables in the paper.

Run in the following order:
 
1) output4paper_figures_tables.m
2) figures_tables_4paper.do

DATA INPUTS: data_laborflows_timeseries.dta is the data file with time series of labor and value added statistics for each country 

EXPECTED RUN TIME: Less than a minute



