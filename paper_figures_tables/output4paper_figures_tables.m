clear 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Matlab file for "Non-traded goods, factor market frictions, and
%international capital flows", Jacek Rothert and Jacob Short.
%
%
%THIS FILE WILL COMPILE THE CROSS COUNTRY DATA FROM
% DATA, CALIBRATION, AND COUNTERFACTUALS NEEDED TO
% GENERATE THE FIGURES AND TABLES FOR THE PAPER
%
% Input files:
% 1. results_counterfactuals_v2.mat - contains counterfactual results
% 2. Calibration_Results_teta01_ppsi1000_v2.mat - contains benchmark
% calibration
%
% Output files:
% 1. modeldata_4stata.xlsx - output of observed, calibrated, and
% counterfactual data to be read into Stata
% 2. LNtimeseries_bench.xlsx, LNtimeseries_psiL.xlsx,
% LNtimeseries_tauL.xlsx - output of the time series for labor from
% benchmark model and counterfactuals to be read into Stata.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%main directory 
mainpath='{INSERT MAIN PATH CONTAINING MODEL CODE';

%Calibration and Results Path
resultsfolder=strcat(mainpath,'{INSERT PATH CONTAINING RESULTS FROM CALIBRATION AND COUNTERFACTUALS}');

%Output Path for excel files to be used to generate figures and tables
exceloutput=strcat(resultsfolder,'{INSERT PATH TO STORE OUTPUT OF EXCEL FILES}');


%create matrices which will be outputed to generate model fit scatter plots
cd(resultsfolder)

%counterfactual files contains all results
%NOTE BEFORE LOADING RESULTS COUNTERFACTUALS FOR OUTPUTING DATA
% WE NEED TO LOAD results_counterfactuals_v2.mat AND THEN
% LOAD Calibration_Results_teta01_ppsi1000_v2.mat IN ORDER TO 
% REPLACE THE Calibrated_Table AND OTHER CALIBRATED DATA SO THAT IT
% IS THE BENCHMARK CALIBRATION. 
% results_counterfactuals_v2.mat ENDED WITH THE CALIBRATION FROM THETA=2.
load('results_counterfactuals_v2.mat')
load('Calibration_Results_teta01_ppsi1000_v2.mat')
save('Calibration_Counterfactuals_Results_DDY0.mat')



%SET UP THE CROSS-COUNTRY DATA THAT WILL BE OUTPUTED
%FIGURES in paper
%   1. DDYT v. dydata - need DDYT_data, dydata
%   2. Private flows v dydata - DDYTpriv
%   3. debt price wedges v dydata - tauScost for 1 sector and 2 sector
%   4/5. Employment and VA shares v RGDP per capita and dydata
%           -LN0, LNT data and VAnontraded0 
%           -RGDP per capita and VAnontradedT needed from data files
%   6. model fit
%       -dydata and dy_cal
%       -dpdata and drer_cal
%       -IoverY and IoverY_cal
%       -DDYT and DDYT_cal
%   7. gT/gN v. dydata - gT_cal and gN_cal
%   8. Counterfactual inflows
%       -DDYT data
%       -Tau_S CF - one sector, two sector  
%       -psiL CF - one sector, two sector
%       -should store a couple other CF flows which may want to be used
%   9. debt price wedge v dydata
%       -tauScost from one sector, two sector
%   10. labor allocation wedge v dy
%       -avgtauL for benchmark and other CFs
%   11. counterfactual labor
%       -LN0,LN1,LN2,LNT, maxdLN for benchmark and counterfactuals
%       -yN0,yN1,yNT and maxdyN for benchmark and CFs
%
%TABLES in paper
%   1. debt price wedge and actual and predicted flows
%       -tauScost for one sector, two sector bench and CFs
%   2. RGDP growth and labor flows
%       -LN0, LNT data and dydata
%   3. Parameterization (nothing here)
%   4. MAIN RESULTS TABLE
%       -DDYT for data, one sector, two sector bench and CFs
%       -DDYT for other specifications and CFs
%       -correl(DDYT,dy)
%       -SSE DDYT 
%       -min/max/avg - tauScost for 1 sector, 2 sector bench and CFs and
%       other specifications and CFs
%       -correl(tauS,dy)
%       -correl(S,I) - SoverY and IoverY for all models and CFs
%
%
%



colnames={'iso','year1','yearT','TT',...  %data
    'dy_data','drer_data','DDY0_data','IoverY_data','SoverY_data',...  %data
    'DDY0priv_data','LN0_data','LNT_data','yN0_data', ... %data
    'DDY0_1sec','dy_1sec','IoverY_1sec','SoverY_1sec','tauScost_1sec','tauS_1sec',...  %one sector
    'DDY0tauS_1sec','dytauS_1sec','IoverYtauS_1sec','SoverYtauS_1sec',...  %one sector CF
    'DDY0tauSK_1sec','dytauSK_1sec','IoverYtauSK_1sec','SoverYtauSK_1sec',...  %one sector tauSKCF
    'dy_bench','drer_bench','DDY0_bench','IoverY_bench','SoverY_bench',... %benchmark
    'SIcorr_bench','gT_bench','gN_bench','tauS_bench','tauK_bench',... %benchmark model
    'tauScost_bench','avgtauL_bench','maxdLN_bench',... %benchmark model
    'yN0_bench','yN1_bench','yNT_bench','maxdyN_bench',... %benchmark model
    'DDY0tauS_bench','IoverYtauS_bench','SoverYtauS_bench','dytauS_bench',... %bench tauS CF
    'yN0tauS_bench','yN1tauS_bench','yNTtauS_bench','maxdyNtauS_bench',... %bench tauS CF
    'DDY0tauSK_bench','IoverYtauSK_bench','SoverYtauSK_bench',... %bench tauSK CF
    'yN0tauSK_bench','yN1tauSK_bench','yNTtauSK_bench','dytauSK_bench',... %bench tauSK CF
    'maxdyNtauSK_bench',...
    'DDY0tauL_bench','LN0tauL_bench','LN1tauL_bench','LN2tauL_bench',... %bench tauL CF
    'LNTtauL_bench','maxdLNtauL_bench',... %bench tauL CF 
    'yN0tauL_bench','yN1tauL_bench','yNTtauL_bench','maxdyNtauL_bench',... %bench tauL CF
    'DDY0psi_bench','LN0psi_bench','LN1psi_bench','LN2psi_bench',... %bench psi CF
    'LNTpsi_bench','maxdLNpsi_bench',... %bench psi CF 
    'yN0psi_bench','yN1psi_bench','yNTpsi_bench','maxdyNpsi_bench',... %bench psi CF
    'DDY0psiL_bench','IoverYpsiL_bench','SoverYpsiL_bench',... %bench psiL CF
    'dypsiL_bench','tauScostpsiL_bench','avgtauLpsiL_bench',... %bench psiL CF
    'LN0psiL_bench','LN1psiL_bench','LN2psiL_bench','LNTpsiL_bench',... %bench psiL CF
    'maxdLNpsiL_bench',... %bench psiL CF 
    'yN0psiL_bench','yN1psiL_bench','yNTpsiL_bench','maxdyNpsiL_bench',... %bench psiL CF 
    'DDY0psiSKL_bench','IoverYpsiSKL_bench','SoverYpsiSKL_bench','dypsiSKL_bench',... %bench psiSKL CF
    'tauScostpsiSKL_bench','avgtauLpsiSKL_bench',... %bench psiSKL CF
    'LN0psiSKL_bench','LN1psiSKL_bench','LN2psiSKL_bench','LNTpsiSKL_bench',... %bench psiSKL CF
    'maxdLNpsiSKL_bench',... %bench psiSKL CF
    'yN0psiSKL_bench','yN1psiSKL_bench','yNTpsiSKL_bench','maxdyNpsiSKL_bench',... %bench psiSKL CF
    'DDY0tauS_theta05','IoverYtauS_theta05','SoverYtauS_theta05','dytauS_theta05',... %Theta05 tauS CF
    'yN0tauS_theta05','yN1tauS_theta05','yNTtauS_theta05','maxdyNtauS_theta05',... %Theta05 tauS CF
    'DDY0psiL_theta05','IoverYpsiL_theta05','SoverYpsiL_theta05','dypsiL_theta05',... %Theta05 psiL CF
    'tauScostpsiL_theta05','avgtauLpsiL_theta05',... %Theta05 psiL CF
    'LN0psiL_theta05','LN1psiL_theta05','LN2psiL_theta05','LNTpsiL_theta05',... %Theta05 psiL CF
    'maxdLNpsiL_theta05',... %Theta05 psiL CF
    'yN0psiL_theta05','yN1psiL_theta05','yNTpsiL_theta05','maxdyNpsiL_theta05',... %Theta05 psiL CF     
    'DDY0tauS_theta2','IoverYtauS_theta2','SoverYtauS_theta2','dytauS_theta2',... %Theta2 tauS CF
    'yN0tauS_theta2','yN1tauS_theta2','yNTtauS_theta2','maxdyNtauS_theta2',... %Theta2 tauS CF
    'DDY0psiL_theta2','IoverYpsiL_theta2','SoverYpsiL_theta2','dypsiL_theta2',... %Theta2 psiL CF
    'tauScostpsiL_theta2','avgtauLpsiL_theta2',... %Theta2 psiL CF
    'LN0psiL_theta2','LN1psiL_theta2','LN2psiL_theta2','LNTpsiL_theta2',... %Theta2 psiL CF
    'maxdLNpsiL_theta2',... %Theta2 psiL CF
    'yN0psiL_theta2','yN1psiL_theta2','yNTpsiL_theta2','maxdyNpsiL_theta2',... %Theta2 psiL CF 
    'DDY0tauS_omc05','IoverYtauS_omc05','SoverYtauS_omc05','dytauS_omc05',... %omc05 tauS CF
    'yN0tauS_omc05','yN1tauS_omc05','yNTtauS_omc05','maxdyNtauS_omc05',... %omc05 tauS CF
    'DDY0psiL_omc05','IoverYpsiL_omc05','SoverYpsiL_omc05','dypsiL_omc05',... %omc05 psiL CF
    'tauScostpsiL_omc05','avgtauLpsiL_omc05',... %omc05 psiL CF
    'LN0psiL_omc05','LN1psiL_omc05','LN2psiL_omc05','LNTpsiL_omc05',... %omc05 psiL CF
    'maxdLNpsiL_omc05',... %omc05 psiL CF
    'yN0psiL_omc05','yN1psiL_omc05','yNTpsiL_omc05','maxdyNpsiL_omc05',... %omc05 psiL CF     
    };

modeldata_4paper=array2table(zeros(size(sum_stats,1),length(colnames)),'VariableNames',colnames);




%% DATA 
%'iso','year1','yearT','TT',...  %data
%    'dy_data','drer_data','DDY0_data','IoverY_data','SoverY_data',...  %data
%    'DDY0priv_data','LN0_data','LNT_data','yN0_data', ... %data
    

modeldata_4paper.iso=sum_stats.iso;
modeldata_4paper.year1=sum_stats.year1;
modeldata_4paper.yearT=sum_stats.yearT;
modeldata_4paper.TT=sum_stats.yearT-sum_stats.year1+1;
modeldata_4paper.dy_data=sum_stats.dydata;
modeldata_4paper.drer_data=sum_stats.dpdata;
modeldata_4paper.DDY0_data=sum_stats.DD_Y0;
modeldata_4paper.IoverY_data=sum_stats.IoverY;
modeldata_4paper.SoverY_data=sum_stats.SoverY;
modeldata_4paper.DDY0priv_data=sum_stats.DDpriv_Y0;
modeldata_4paper.LN0_data=sum_stats.LN0;
modeldata_4paper.LNT_data=sum_stats.LNT;
modeldata_4paper.yN0_data=sum_stats.pN0yN0_over_gdp0;

%RGDP per capita and yNT need to be pulled from data files.

%%




%% One Sector
%    'DDY0_1sec','dy_1sec','IoverY_1sec','SoverY_1sec','tauScost_1sec','tauS_1sec',...  %one sector
%    'DDY0tauS_1sec','dytauS_1sec','IoverYtauS_1sec','SoverYtauS_1sec',...  %one sector CF

modeldata_4paper.DDY0_1sec=OneSec_BENCH_DD_Y0;
modeldata_4paper.dy_1sec=OneSec_BENCH_dy;
modeldata_4paper.IoverY_1sec=OneSec_BENCH_IoverY;
modeldata_4paper.SoverY_1sec=OneSec_BENCH_SoverY;
modeldata_4paper.tauScost_1sec=OneSec_BENCH_tauScost;
modeldata_4paper.tauS_1sec=OneSec_BENCH_tauS;
modeldata_4paper.DDY0tauS_1sec=OneSec_COUNTER_tauS_DD_Y0;
modeldata_4paper.dytauS_1sec=OneSec_COUNTER_tauS_dy;
modeldata_4paper.IoverYtauS_1sec=OneSec_COUNTER_tauS_IoverY;
modeldata_4paper.SoverYtauS_1sec=OneSec_COUNTER_tauS_SoverY;
modeldata_4paper.DDY0tauSK_1sec=OneSec_COUNTER_tauSK_DD_Y0;
modeldata_4paper.dytauSK_1sec=OneSec_COUNTER_tauSK_dy;
modeldata_4paper.IoverYtauSK_1sec=OneSec_COUNTER_tauSK_IoverY;
modeldata_4paper.SoverYtauSK_1sec=OneSec_COUNTER_tauSK_SoverY;


%%


%% Benchmark
%    'dy_bench','drer_bench','DDY0_bench','IoverY_bench','SoverY_bench',... %benchmark
%    'SIcorr_bench','gT_bench','gN_bench','tauS_bench','tauK_bench',... %benchmark model
%    'tauScost_bench','avgtauL_bench','maxdLN_bench',... %benchmark model
%    'yN0_bench','yN1_bench','yNT_bench','maxdyN_bench',

modeldata_4paper.dy_bench=crosscountrydata.Bench.calib.dy;
modeldata_4paper.drer_bench=crosscountrydata.Bench.calib.drer;
modeldata_4paper.DDY0_bench=crosscountrydata.Bench.calib.DDY0;
modeldata_4paper.IoverY_bench=crosscountrydata.Bench.calib.IoverY;
modeldata_4paper.SoverY_bench=crosscountrydata.Bench.calib.SoverY;
modeldata_4paper.SIcorr_bench=Calibrated_Table.SIcorr_cal;
modeldata_4paper.gT_bench=Calibrated_Table.gT_cal;
modeldata_4paper.gN_bench=Calibrated_Table.gN_cal;
modeldata_4paper.tauS_bench=Calibrated_Table.tauS_cal;
modeldata_4paper.tauK_bench=Calibrated_Table.tauK_cal;
modeldata_4paper.tauScost_bench=crosscountrydata.Bench.calib.tauScost;
modeldata_4paper.avgtauL_bench=crosscountrydata.Bench.calib.avgtauL;
%replace with the following once the counterfactual files are updated with
%new output
modeldata_4paper.maxdLN_bench=crosscountrydata.Bench.calib.maxdLN;
modeldata_4paper.yN0_bench=crosscountrydata.Bench.calib.yN0;
modeldata_4paper.yN1_bench=crosscountrydata.Bench.calib.yN1;
modeldata_4paper.yNT_bench=crosscountrydata.Bench.calib.yNT;
modeldata_4paper.maxdyN_bench=crosscountrydata.Bench.calib.maxdyN;


%%


%% Benchmark tauS CF
%    'DDY0tauS_bench','IoverYtauS_bench','SoverYtauS_bench',... %bench tauS CF
%    'yN0tauS_bench','yN1tauS_bench','yNTtauS_bench','dytauS_bench',... %bench tauS CF
%    'maxdyNtauS_bench',...
    
modeldata_4paper.DDY0tauS_bench=crosscountrydata.Bench.tau_S.DDY0;
modeldata_4paper.IoverYtauS_bench=crosscountrydata.Bench.tau_S.IoverY;
modeldata_4paper.SoverYtauS_bench=crosscountrydata.Bench.tau_S.SoverY;
modeldata_4paper.dytauS_bench=crosscountrydata.Bench.tau_S.dy;
modeldata_4paper.yN0tauS_bench=crosscountrydata.Bench.tau_S.yN0;
modeldata_4paper.yN1tauS_bench=crosscountrydata.Bench.tau_S.yN1;
modeldata_4paper.yNTtauS_bench=crosscountrydata.Bench.tau_S.yNT;
modeldata_4paper.maxdyNtauS_bench=crosscountrydata.Bench.tau_S.maxdyN;



%%

%% Benchmark tauSK CF
%    'DDY0tauSK_bench','IoverYtauSK_bench','SoverYtauSK_bench',... %bench tauSK CF
%    'yN0tauSK_bench','yN1tauSK_bench','yNTtauSK_bench','dytauSK_bench',... %bench tauSK CF
%    'maxdyNtauSK_bench',...
    
modeldata_4paper.DDY0tauSK_bench=crosscountrydata.Bench.tau_SK.DDY0;
modeldata_4paper.IoverYtauSK_bench=crosscountrydata.Bench.tau_SK.IoverY;
modeldata_4paper.SoverYtauSK_bench=crosscountrydata.Bench.tau_SK.SoverY;
modeldata_4paper.dytauSK_bench=crosscountrydata.Bench.tau_SK.dy;
modeldata_4paper.yN0tauSK_bench=crosscountrydata.Bench.tau_SK.yN0;
modeldata_4paper.yN1tauSK_bench=crosscountrydata.Bench.tau_SK.yN1;
modeldata_4paper.yNTtauSK_bench=crosscountrydata.Bench.tau_SK.yNT;
modeldata_4paper.maxdyNtauSK_bench=crosscountrydata.Bench.tau_SK.maxdyN;



%%



%% Benchmark tauL CF
%    'DDY0tauL_bench','LN0tauL_bench','LN1tauL_bench','LN2tauL_bench',... %bench tauL CF
%    'LNTtauL_bench','maxdLNtauL_bench',... %bench tauL CF 
%    'yN0tauL_bench','yN1tauL_bench','yNTtauL_bench','maxdyNtauL_bench',... %bench tauL CF
    
modeldata_4paper.DDY0tauL_bench=crosscountrydata.Bench.tau_L.DDY0;
modeldata_4paper.LN0tauL_bench=crosscountrydata.Bench.tau_L.LN0;
modeldata_4paper.LN1tauL_bench=crosscountrydata.Bench.tau_L.LN1;
modeldata_4paper.LN2tauL_bench=crosscountrydata.Bench.tau_L.LN2;
modeldata_4paper.LNTtauL_bench=crosscountrydata.Bench.tau_L.LNT;
modeldata_4paper.maxdLNtauL_bench=crosscountrydata.Bench.tau_L.maxdLN;
modeldata_4paper.yN0tauL_bench=crosscountrydata.Bench.tau_L.yN0;
modeldata_4paper.yN1tauL_bench=crosscountrydata.Bench.tau_L.yN1;
modeldata_4paper.yNTtauL_bench=crosscountrydata.Bench.tau_L.yNT;
modeldata_4paper.maxdyNtauL_bench=crosscountrydata.Bench.tau_L.maxdyN;



%%



%% Benchmark psi CF
%    'DDY0psi_bench','LN0psi_bench','LN1psi_bench','LN2psi_bench',... %bench tauL CF
%    'LNTpsi_bench','maxdLNpsi_bench',... %bench psi CF 
%    'yN0psi_bench','yN1psi_bench','yNTpsi_bench','maxdyNpsi_bench',... %bench tauL CF
    
modeldata_4paper.DDY0psi_bench=crosscountrydata.Bench.psi.DDY0;
modeldata_4paper.LN0psi_bench=crosscountrydata.Bench.psi.LN0;
modeldata_4paper.LN1psi_bench=crosscountrydata.Bench.psi.LN1;
modeldata_4paper.LN2psi_bench=crosscountrydata.Bench.psi.LN2;
modeldata_4paper.LNTpsi_bench=crosscountrydata.Bench.psi.LNT;
modeldata_4paper.maxdLNpsi_bench=crosscountrydata.Bench.psi.maxdLN;
modeldata_4paper.yN0psi_bench=crosscountrydata.Bench.psi.yN0;
modeldata_4paper.yN1psi_bench=crosscountrydata.Bench.psi.yN1;
modeldata_4paper.yNTpsi_bench=crosscountrydata.Bench.psi.yNT;
modeldata_4paper.maxdyNpsi_bench=crosscountrydata.Bench.psi.maxdyN;



%%




%% Benchmark psiL CF
%    'DDY0psiL_bench','IoverYpsiL_bench','SoverYpsiL_bench','dypsiL_bench',... %bench psiL CF
%    'tauScostpsiL_bench','avgtauLpsiL_bench',... %bench psiL CF
%    'LN0psiL_bench','LN1psiL_bench','LN2psiL_bench','LNTpsiL_bench',... %bench psiL CF
%    'maxdLNpsiL_bench',... %bench psiL CF 
%    'yN0psiL_bench','yN1psiL_bench','yNTpsiL_bench','maxdyNpsiL_bench',... %bench psiL CF 

modeldata_4paper.DDY0psiL_bench=crosscountrydata.Bench.psiL.DDY0;
modeldata_4paper.IoverYpsiL_bench=crosscountrydata.Bench.psiL.IoverY;
modeldata_4paper.SoverYpsiL_bench=crosscountrydata.Bench.psiL.SoverY;
modeldata_4paper.dypsiL_bench=crosscountrydata.Bench.psiL.dy;
modeldata_4paper.tauScostpsiL_bench=crosscountrydata.Bench.psiL.tauScost;
modeldata_4paper.avgtauLpsiL_bench=crosscountrydata.Bench.psiL.avgtauL;
modeldata_4paper.LN0psiL_bench=crosscountrydata.Bench.psiL.LN0;
modeldata_4paper.LN1psiL_bench=crosscountrydata.Bench.psiL.LN1;
modeldata_4paper.LN2psiL_bench=crosscountrydata.Bench.psiL.LN2;
modeldata_4paper.LNTpsiL_bench=crosscountrydata.Bench.psiL.LNT;
modeldata_4paper.maxdLNpsiL_bench=crosscountrydata.Bench.psiL.maxdLN;
modeldata_4paper.yN0psiL_bench=crosscountrydata.Bench.psiL.yN0;
modeldata_4paper.yN1psiL_bench=crosscountrydata.Bench.psiL.yN1;
modeldata_4paper.yNTpsiL_bench=crosscountrydata.Bench.psiL.yNT;
modeldata_4paper.maxdyNpsiL_bench=crosscountrydata.Bench.psiL.maxdyN;


%%

%% Benchmark psiSKL CF
%    'DDY0psiSKL_bench','IoverYpsiSKL_bench','SoverYpsiSKL_bench',... %bench psiSKL CF
%    'tauScostpsiSKL_bench','avgtauLpsiSKL_bench',... %bench psiSKL CF
%    'LN0psiSKL_bench','LN1psiSKL_bench','LN2psiSKL_bench','LNTpsiSKL_bench',... %bench psiSKL CF
%    'maxdLNpsiSKL_bench',... %bench psiL CF
%    'yN0psiSKL_bench','yN1psiSKL_bench','yNTpsiSKL_bench','maxdyNpsiSKL_bench',... %bench psiSKL CF

modeldata_4paper.DDY0psiSKL_bench=crosscountrydata.Bench.psiSKL.DDY0;
modeldata_4paper.IoverYpsiSKL_bench=crosscountrydata.Bench.psiSKL.IoverY;
modeldata_4paper.SoverYpsiSKL_bench=crosscountrydata.Bench.psiSKL.SoverY;
modeldata_4paper.dypsiSKL_bench=crosscountrydata.Bench.psiSKL.dy;
modeldata_4paper.tauScostpsiSKL_bench=crosscountrydata.Bench.psiSKL.tauScost;
modeldata_4paper.avgtauLpsiSKL_bench=crosscountrydata.Bench.psiSKL.avgtauL;
modeldata_4paper.LN0psiSKL_bench=crosscountrydata.Bench.psiSKL.LN0;
modeldata_4paper.LN1psiSKL_bench=crosscountrydata.Bench.psiSKL.LN1;
modeldata_4paper.LN2psiSKL_bench=crosscountrydata.Bench.psiSKL.LN2;
modeldata_4paper.LNTpsiSKL_bench=crosscountrydata.Bench.psiSKL.LNT;
modeldata_4paper.maxdLNpsiSKL_bench=crosscountrydata.Bench.psiSKL.maxdLN;
modeldata_4paper.yN0psiSKL_bench=crosscountrydata.Bench.psiSKL.yN0;
modeldata_4paper.yN1psiSKL_bench=crosscountrydata.Bench.psiSKL.yN1;
modeldata_4paper.yNTpsiSKL_bench=crosscountrydata.Bench.psiSKL.yNT;
modeldata_4paper.maxdyNpsiSKL_bench=crosscountrydata.Bench.psiSKL.maxdyN;


%%


%% Theta=0.5 tauS CF
%    'DDY0tauS_theta05','IoverYtauS_theta05','SoverYtauS_theta05',... %Theta05 tauS CF
%    'yN0tauS_theta05','yN1tauS_theta05','yNTtauS_theta05','maxdyNtauS_theta05',... %Theta05 tauS CF
    
modeldata_4paper.DDY0tauS_theta05=crosscountrydata.Theta05.tau_S.DDY0;
modeldata_4paper.IoverYtauS_theta05=crosscountrydata.Theta05.tau_S.IoverY;
modeldata_4paper.SoverYtauS_theta05=crosscountrydata.Theta05.tau_S.SoverY;
modeldata_4paper.dytauS_theta05=crosscountrydata.Theta05.tau_S.dy;
modeldata_4paper.yN0tauS_theta05=crosscountrydata.Theta05.tau_S.yN0;
modeldata_4paper.yN1tauS_theta05=crosscountrydata.Theta05.tau_S.yN1;
modeldata_4paper.yNTtauS_theta05=crosscountrydata.Theta05.tau_S.yNT;
modeldata_4paper.maxdyNtauS_theta05=crosscountrydata.Theta05.tau_S.maxdyN;


%%


%% Theta=0.5 psiL CF
%    'DDY0psiL_theta05','IoverYpsiL_theta05','SoverYpsiL_theta05',... %Theta05 psiL CF
%    'tauScostpsiL_theta05','avgtauLpsiL_theta05',... %Theta05 psiL CF
%    'LN0psiL_theta05','LN1psiL_theta05','LN2psiL_theta05','LNTpsiL_theta05',... %Theta05 psiL CF
%    'maxdLNpsiL_theta05',... %Theta05 psiL CF
%    'yN0psiL_theta05','yN1psiL_theta05','yNTpsiL_theta05','maxdyNpsiL_theta05',... %Theta05 psiL CF     

modeldata_4paper.DDY0psiL_theta05=crosscountrydata.Theta05.psiL.DDY0;
modeldata_4paper.IoverYpsiL_theta05=crosscountrydata.Theta05.psiL.IoverY;
modeldata_4paper.SoverYpsiL_theta05=crosscountrydata.Theta05.psiL.SoverY;
modeldata_4paper.dypsiL_theta05=crosscountrydata.Theta05.psiL.dy;
modeldata_4paper.tauScostpsiL_theta05=crosscountrydata.Theta05.psiL.tauScost;
modeldata_4paper.avgtauLpsiL_theta05=crosscountrydata.Theta05.psiL.avgtauL;
modeldata_4paper.LN0psiL_theta05=crosscountrydata.Theta05.psiL.LN0;
modeldata_4paper.LN1psiL_theta05=crosscountrydata.Theta05.psiL.LN1;
modeldata_4paper.LN2psiL_theta05=crosscountrydata.Theta05.psiL.LN2;
modeldata_4paper.LNTpsiL_theta05=crosscountrydata.Theta05.psiL.LNT;
modeldata_4paper.maxdLNpsiL_theta05=crosscountrydata.Theta05.psiL.maxdLN;
modeldata_4paper.yN0psiL_theta05=crosscountrydata.Theta05.psiL.yN0;
modeldata_4paper.yN1psiL_theta05=crosscountrydata.Theta05.psiL.yN1;
modeldata_4paper.yNTpsiL_theta05=crosscountrydata.Theta05.psiL.yNT;
modeldata_4paper.maxdyNpsiL_theta05=crosscountrydata.Theta05.psiL.maxdyN;



%%


%% Theta=2 tauS CF
%    'DDY0tauS_theta2','IoverYtauS_theta2','SoverYtauS_theta2',... %Theta2 tauS CF
%    'yN0tauS_theta2','yN1tauS_theta2','yNTtauS_theta2','maxdyNtauS_theta2',... %Theta2 tauS CF   

modeldata_4paper.DDY0tauS_theta2=crosscountrydata.Theta2.tau_S.DDY0;
modeldata_4paper.IoverYtauS_theta2=crosscountrydata.Theta2.tau_S.IoverY;
modeldata_4paper.SoverYtauS_theta2=crosscountrydata.Theta2.tau_S.SoverY;
modeldata_4paper.dytauS_theta2=crosscountrydata.Theta2.tau_S.dy;
modeldata_4paper.yN0tauS_theta2=crosscountrydata.Theta2.tau_S.yN0;
modeldata_4paper.yN1tauS_theta2=crosscountrydata.Theta2.tau_S.yN1;
modeldata_4paper.yNTtauS_theta2=crosscountrydata.Theta2.tau_S.yNT;
modeldata_4paper.maxdyNtauS_theta2=crosscountrydata.Theta2.tau_S.maxdyN;




%%


%% Theta=2 psiL CF
%    'DDY0psiL_theta2','IoverYpsiL_theta2','SoverYpsiL_theta2',... %Theta2 psiL CF
%    'tauScostpsiL_theta2','avgtauLpsiL_theta2',... %Theta2 psiL CF
%    'LN0psiL_theta2','LN1psiL_theta2','LN2psiL_theta2','LNTpsiL_theta2',... %Theta2 psiL CF
%    'maxdLNpsiSKL_theta2',... %Theta2 psiL CF
%    'yN0psiL_theta2','yN1psiL_theta2','yNTpsiL_theta2','maxdyNpsiL_theta2',... %Theta2 psiL CF 

modeldata_4paper.DDY0psiL_theta2=crosscountrydata.Theta2.psiL.DDY0;
modeldata_4paper.IoverYpsiL_theta2=crosscountrydata.Theta2.psiL.IoverY;
modeldata_4paper.SoverYpsiL_theta2=crosscountrydata.Theta2.psiL.SoverY;
modeldata_4paper.dypsiL_theta2=crosscountrydata.Theta2.psiL.dy;
modeldata_4paper.tauScostpsiL_theta2=crosscountrydata.Theta2.psiL.tauScost;
modeldata_4paper.avgtauLpsiL_theta2=crosscountrydata.Theta2.psiL.avgtauL;
modeldata_4paper.LN0psiL_theta2=crosscountrydata.Theta2.psiL.LN0;
modeldata_4paper.LN1psiL_theta2=crosscountrydata.Theta2.psiL.LN1;
modeldata_4paper.LN2psiL_theta2=crosscountrydata.Theta2.psiL.LN2;
modeldata_4paper.LNTpsiL_theta2=crosscountrydata.Theta2.psiL.LNT;
modeldata_4paper.maxdLNpsiL_theta2=crosscountrydata.Theta2.psiL.maxdLN;
modeldata_4paper.yN0psiL_theta2=crosscountrydata.Theta2.psiL.yN0;
modeldata_4paper.yN1psiL_theta2=crosscountrydata.Theta2.psiL.yN1;
modeldata_4paper.yNTpsiL_theta2=crosscountrydata.Theta2.psiL.yNT;
modeldata_4paper.maxdyNpsiL_theta2=crosscountrydata.Theta2.psiL.maxdyN;


%%




%% omc=0.5 tauS CF
%    'DDY0tauS_omc05','IoverYtauS_omc05','SoverYtauS_omc05',... %omc05 tauS CF
%    'yN0tauS_omc05','yN1tauS_omc05','yNTtauS_omc05','maxdyNtauS_omc05',... %omc05 tauS CF
    
modeldata_4paper.DDY0tauS_omc05=crosscountrydata.omc05.tau_S.DDY0;
modeldata_4paper.IoverYtauS_omc05=crosscountrydata.omc05.tau_S.IoverY;
modeldata_4paper.SoverYtauS_omc05=crosscountrydata.omc05.tau_S.SoverY;
modeldata_4paper.dytauS_omc05=crosscountrydata.omc05.tau_S.dy;
modeldata_4paper.yN0tauS_omc05=crosscountrydata.omc05.tau_S.yN0;
modeldata_4paper.yN1tauS_omc05=crosscountrydata.omc05.tau_S.yN1;
modeldata_4paper.yNTtauS_omc05=crosscountrydata.omc05.tau_S.yNT;
modeldata_4paper.maxdyNtauS_omc05=crosscountrydata.omc05.tau_S.maxdyN;


%%


%% omc=0.5 psiL CF
%    'DDY0psiL_omc05','IoverYpsiL_omc05','SoverYpsiL_omc05',... %omc05 psiL CF
%    'tauScostpsiL_omc05','avgtauLpsiL_omc05',... %omc05 psiL CF
%    'LN0psiL_omc05','LN1psiL_omc05','LN2psiL_omc05','LNTpsiL_omc05',... %omc05 psiL CF
%    'maxdLNpsiL_omc05',... %omc05 psiL CF
%    'yN0psiL_omc05','yN1psiL_omc05','yNTpsiL_omc05','maxdyNpsiL_omc05',... %omc05 psiL CF     

modeldata_4paper.DDY0psiL_omc05=crosscountrydata.omc05.psiL.DDY0;
modeldata_4paper.IoverYpsiL_omc05=crosscountrydata.omc05.psiL.IoverY;
modeldata_4paper.SoverYpsiL_omc05=crosscountrydata.omc05.psiL.SoverY;
modeldata_4paper.dypsiL_omc05=crosscountrydata.omc05.psiL.dy;
modeldata_4paper.tauScostpsiL_omc05=crosscountrydata.omc05.psiL.tauScost;
modeldata_4paper.avgtauLpsiL_omc05=crosscountrydata.omc05.psiL.avgtauL;
modeldata_4paper.LN0psiL_omc05=crosscountrydata.omc05.psiL.LN0;
modeldata_4paper.LN1psiL_omc05=crosscountrydata.omc05.psiL.LN1;
modeldata_4paper.LN2psiL_omc05=crosscountrydata.omc05.psiL.LN2;
modeldata_4paper.LNTpsiL_omc05=crosscountrydata.omc05.psiL.LNT;
modeldata_4paper.maxdLNpsiL_omc05=crosscountrydata.omc05.psiL.maxdLN;
modeldata_4paper.yN0psiL_omc05=crosscountrydata.omc05.psiL.yN0;
modeldata_4paper.yN1psiL_omc05=crosscountrydata.omc05.psiL.yN1;
modeldata_4paper.yNTpsiL_omc05=crosscountrydata.omc05.psiL.yNT;
modeldata_4paper.maxdyNpsiL_omc05=crosscountrydata.omc05.psiL.maxdyN;



%%



%benchmark
filename='modeldata_4stata.xlsx';

cd(exceloutput)
writetable(modeldata_4paper,filename)





LNtimeseries_bench=array2table(zeros(size(sum_stats,1),71));

LNtimeseries_bench.Var1=sum_stats.iso;

for j=1:70

    varnum=strcat('Var',num2str(j+1));
    
   LNtimeseries_bench.(varnum)=transpose(crosscountry_laborflows.Bench.calib(j,:)); 
    
end



LNtimeseries_psiL=array2table(zeros(size(sum_stats,1),71));

LNtimeseries_psiL.Var1=sum_stats.iso;

for j=1:70

    varnum=strcat('Var',num2str(j+1));
    
   LNtimeseries_psiL.(varnum)=transpose(crosscountry_laborflows.Bench.psiL(j,:)); 
    
end


LNtimeseries_tauL=array2table(zeros(size(sum_stats,1),71));

LNtimeseries_tauL.Var1=sum_stats.iso;

for j=1:70

    varnum=strcat('Var',num2str(j+1));
    
   LNtimeseries_tauL.(varnum)=transpose(crosscountry_laborflows.Bench.tau_L(j,:)); 
    
end



cd(exceloutput)
writetable(LNtimeseries_bench,'LNtimeseries_bench.xlsx')
writetable(LNtimeseries_psiL,'LNtimeseries_psiL.xlsx')
writetable(LNtimeseries_tauL,'LNtimeseries_tauL.xlsx')



