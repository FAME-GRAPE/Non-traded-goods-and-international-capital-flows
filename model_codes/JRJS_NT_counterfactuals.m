%{

Replication file - counterfactuals
"Non-Traded Goods, Factor Markets Frictions, and International Capital Flows"
by Jacek Rothert (USNA), and Jacob Short (BoC)
Review of Economic Dynamics, 2022 (accepted)
(c) JR&JS, September 2021


INPUTS: 
- Calibration_Results_teta01_ppsi1000.mat           - ( output of JRJS_NT_main.m )
- Calibration_Results_teta05_ppsi1000.mat           - ( output of JRJS_NT_main.m ) 
- Calibration_Results_teta2_ppsi1000.mat            - ( output of JRJS_NT_main.m ) 
- Calibration_Results_teta01_ppsi1000_omc05.mat     - ( output of JRJS_NT_main.m ) 

OUTPUTS: 
- OneSec.mat
- results_counterfactuals.mat
- prints Table4 and Table5 in the MATLAB Command Window

DEPENDENCIES
a) calibration_1sector.m
	(a-i) calls eqm_path_1sector.m
b) counterfactual.m
	(b-i) calls eqm_path.m
	(b-ii) calls eqm_path_flex_labor.m

EXPECTED RUN TIME: XYZ hours, (parallel computing, using XYZ nodes)

%}




clear all

%% PATHS
% ==================================================================================================

%calibrate 1sector or not
calibrate_1sector = 1;
%number of cases for counterfactuals: benchmark, theta=0.5, theta=0.2, omc=0.5
num_cases = 4;
mainpath                    = 'C:\Users\R2D2\Dropbox (UW)\research\JRJS_CapitalFlows\submissions\RED\Data_Computing_Codes\model_codes';
% option to have separate folders for Excel, m-files, mat.files with results, etc.
ExcelFolder                 = mainpath;
MatlabFolder                = mainpath;
ResultsFolder               = mainpath;
CountryFolder               = mainpath;
folders_destination.excel   = ExcelFolder;
folders_destination.matlab  = MatlabFolder;
folders_destination.results = ResultsFolder;
folders_destination.country = CountryFolder;

mainpath_source0 = mainpath;
mainpath_source  = mainpath_source0;
filename1sector  = 'OneSec.mat';


cd(mainpath_source)
%save path names to overwrite whenever loading a calibration results file
pathsfilename='CFpathnames.mat';
save(pathsfilename)


cd(mainpath_source)
filename = 'Calibration_Results_teta01_ppsi1000.mat'; % load benchmark
load(filename)
load(pathsfilename)
% ==================================================================================================

cd(folders_destination.matlab)
globals.opcje_iter= optimset('Display','iter','MaxIter',500);


%% ONE-SECTOR MODEL RESULTS
% ====================================================================================
OneSec_BENCH_dy                 = zeros(size(sum_stats,1),1);
OneSec_BENCH_DD_Y0              = zeros(size(sum_stats,1),1);
OneSec_BENCH_IoverY             = zeros(size(sum_stats,1),1);
OneSec_BENCH_SoverY             = zeros(size(sum_stats,1),1);
OneSec_BENCH_tauScost           = zeros(size(sum_stats,1),1);
OneSec_BENCH_flag_calibration   = zeros(size(sum_stats,1),1);
OneSec_BENCH_flag_eqpath        = zeros(size(sum_stats,1),1);
OneSec_BENCH_g                  = zeros(size(sum_stats,1),1);
OneSec_BENCH_tauS               = zeros(size(sum_stats,1),1);
OneSec_BENCH_tauK               = zeros(size(sum_stats,1),1);

OneSec_COUNTER_tauS_dy       = zeros(size(sum_stats,1),1);
OneSec_COUNTER_tauS_DD_Y0    = zeros(size(sum_stats,1),1);
OneSec_COUNTER_tauS_IoverY   = zeros(size(sum_stats,1),1);
OneSec_COUNTER_tauS_SoverY   = zeros(size(sum_stats,1),1);

OneSec_COUNTER_tauK_dy       = zeros(size(sum_stats,1),1);
OneSec_COUNTER_tauK_DD_Y0    = zeros(size(sum_stats,1),1);
OneSec_COUNTER_tauK_IoverY   = zeros(size(sum_stats,1),1);
OneSec_COUNTER_tauK_SoverY   = zeros(size(sum_stats,1),1);

OneSec_COUNTER_tauSK_dy       = zeros(size(sum_stats,1),1);
OneSec_COUNTER_tauSK_DD_Y0    = zeros(size(sum_stats,1),1);
OneSec_COUNTER_tauSK_IoverY   = zeros(size(sum_stats,1),1);
OneSec_COUNTER_tauSK_SoverY   = zeros(size(sum_stats,1),1);

yguesses = zeros(3,size(sum_stats,1));
for i = 1:size(sum_stats,1)
    yguesses(1,i) = sum_stats.dydata(i) - 0.0256;
end

N_lbdagrid = 101;
globals.N_lbdagrid = N_lbdagrid;

globals_fmin            = globals;
globals_fsolve          = globals;
globals_fmin.use_fmin   = 1;
globals_fsolve.use_fmin = 0;

temp_table_1sec = array2table(...
    zeros(size(sum_stats,1),length(colnames)),'VariableNames',colnames,'RowNames',sum_stats.iso'...
    );


OneSec_elast_DDY0 = ones(N_lbdagrid,size(sum_stats,1));
OneSec_elast_tauS = ones(N_lbdagrid,size(sum_stats,1));



if calibrate_1sector == 1
    
    %% set options for parallel computing at the Bank of Canada cluster
    %setenv('MATLAB_WORKER_ARGS','-p long --time=0-96:00:00 --mem=30GB');
    %ppddy0cf1=parpool('edith92',12);
    
    parfor i_country = 1:size(sum_stats,1)
        
        disp(sum_stats.iso(i_country))
        
        %% calibration for the one-sector model
        yguess = yguesses(:,i_country);
        
        f_calibrate_1sector = @(y) calibration_1sector(y,sum_stats(i_country,:),deep_params,globals_fsolve);
        [ youtput, ~, flag_1sec_calibration ] = fsolve(f_calibrate_1sector,yguess,globals_fsolve.opcje_off);
        
        if flag_1sec_calibration ~= 1
            
            f_calibrate_1sector = @(y) calibration_1sector(y,sum_stats(i_country,:),deep_params,globals_fmin);
            [ youtput, ~, flag_1sec_calibration ] = fminsearch(f_calibrate_1sector,yguess,globals_fmin.opcje_off);
            
        end
        
        f_calibrate_1sector = @(y) calibration_1sector(y,sum_stats(i_country,:),deep_params,globals_fsolve);
        [residual,one_sec_results] = f_calibrate_1sector(youtput);
        
        
        %one_sec_results = calibrate_one_sector(sum_stats(i_country,:),deep_params,globals);
        OneSec_BENCH_dy(i_country)                  = one_sec_results.dy;
        OneSec_BENCH_DD_Y0(i_country)               = one_sec_results.DD_Y0;
        OneSec_BENCH_IoverY(i_country)              = one_sec_results.IoverY;
        OneSec_BENCH_SoverY(i_country)              = one_sec_results.SoverY;
        OneSec_BENCH_tauScost(i_country)            = one_sec_results.tauScost;
        OneSec_BENCH_flag_calibration(i_country)    = flag_1sec_calibration;
        OneSec_BENCH_flag_eqpath(i_country)         = one_sec_results.flag_eqpath;
        OneSec_BENCH_g(i_country)                   = one_sec_results.g;
        OneSec_BENCH_tauS(i_country)                = one_sec_results.tauS;
        OneSec_BENCH_tauK(i_country)                = one_sec_results.tauK;
        
        youtput0 = youtput;
        
        %% counter-factual with tauS = 0
        youtput = youtput0;
        youtput(2) = 0;
        [~,one_sec_results] = f_calibrate_1sector(youtput);
        OneSec_COUNTER_tauS_dy(i_country)       = one_sec_results.dy;
        OneSec_COUNTER_tauS_DD_Y0(i_country)    = one_sec_results.DD_Y0;
        OneSec_COUNTER_tauS_IoverY(i_country)   = one_sec_results.IoverY;
        OneSec_COUNTER_tauS_SoverY(i_country)   = one_sec_results.SoverY;
        
        
        %% counter-factual with tauK = 0
        youtput = youtput0;
        youtput(3) = 0;
        [~,one_sec_results] = f_calibrate_1sector(youtput);
        OneSec_COUNTER_tauK_dy(i_country) = one_sec_results.dy;
        OneSec_COUNTER_tauK_DD_Y0(i_country) = one_sec_results.DD_Y0;
        OneSec_COUNTER_tauK_IoverY(i_country) = one_sec_results.IoverY;
        OneSec_COUNTER_tauK_SoverY(i_country) = one_sec_results.SoverY;
        
        %% counter-factual with tauS = tauK = 0
        youtput = youtput0;
        youtput(2) = 0;  youtput(3) = 0;
        [~,one_sec_results] = f_calibrate_1sector(youtput);
        
        OneSec_COUNTER_tauSK_dy(i_country) = one_sec_results.dy;
        OneSec_COUNTER_tauSK_DD_Y0(i_country) = one_sec_results.DD_Y0;
        OneSec_COUNTER_tauSK_IoverY(i_country) = one_sec_results.IoverY;
        OneSec_COUNTER_tauSK_SoverY(i_country) = one_sec_results.SoverY;
        
        taugrid = linspace(...
            min(0,youtput0(2)),...
            max(0,youtput0(2)), N_lbdagrid);
        youtput = youtput0;
        for i_tau = 1:N_lbdagrid
            %disp([i_country,i_tau])
            youtput(2) = taugrid(i_tau);
            [~,one_sec_results] = f_calibrate_1sector(youtput);
            % record DDY0 vs tauS
            OneSec_elast_DDY0(i_tau,i_country) = one_sec_results.DD_Y0;
            OneSec_elast_tauS(i_tau,i_country) = youtput(2);
        end
        
        
    end
    
    % clear the local variable from parallel computing at the Bank of Canada cluster
    %delete(ppddy0cf1);

    cd(mainpath_source)
    save(filename1sector)
    
else
    cd(mainpath_source)
    load(filename1sector)
    load(pathsfilename)
end
%}

mainpath_source = mainpath_source0;

    
    

%% table with one-sector model results
% ==================================================================================================
colnames = {'iso','flag_calibration','flag_eqpath',...
    'dy_err','IoverY_err','DD_Y0_err',...
    'dy_data','IoverY_data','DD_Y0_data',...
    'g','tauS','tauK',...
    'tauScost',...
    'SoverY','IoverY'};
rownames = sum_stats.iso;
Calibration_Table_1sector               = ...
    array2table(NaN(length(rownames),length(colnames)),...
    'VariableNames',colnames,...
    'RowNames',rownames);
Calibration_Table_1sector.iso           = sum_stats.iso;
Calibration_Table_1sector.flag_calibration   = OneSec_BENCH_flag_calibration;
Calibration_Table_1sector.flag_eqpath   = OneSec_BENCH_flag_eqpath;
Calibration_Table_1sector.g             = OneSec_BENCH_g;
Calibration_Table_1sector.tauS          = OneSec_BENCH_tauS;
Calibration_Table_1sector.tauK          = OneSec_BENCH_tauK;
Calibration_Table_1sector.dy_err        = OneSec_BENCH_dy - sum_stats.dydata;
Calibration_Table_1sector.IoverY_err    = OneSec_BENCH_IoverY - sum_stats.IoverY;
Calibration_Table_1sector.DD_Y0_err     = OneSec_BENCH_DD_Y0 - sum_stats.DD_Y0;
Calibration_Table_1sector.dy_data       = sum_stats.dydata;
Calibration_Table_1sector.IoverY_data   = sum_stats.IoverY;
Calibration_Table_1sector.DD_Y0_data    = sum_stats.DD_Y0;
Calibration_Table_1sector.tauScost      = OneSec_BENCH_tauScost;
Calibration_Table_1sector.SoverY        = OneSec_BENCH_SoverY;
Calibration_Table_1sector.IoverY        = OneSec_BENCH_IoverY;

% ====================================================================================
%
%
disp(Calibration_Table_1sector)

if(num_cases==4)
%ALL COUNTERFACTUALS ARE RUN
cases = {...
    'teta01_ppsi1000','teta05_ppsi1000','teta2_ppsi1000','teta01_ppsi1000_omc05'};
elseif(num_cases==3)
%ALL COUNTERFACTUALS ARE RUN
cases = {...
    'teta01_ppsi1000','teta05_ppsi1000','teta2_ppsi1000'};
elseif(num_cases==2)
%2 COUNTERFACTUALS ARE RUN
cases = {...
    'teta01_ppsi1000','teta05_ppsi1000'};
elseif(num_cases==1)
cases = {...
    'teta01_ppsi1000'};
end

counter_types = {'tau_S','tau_K','tau_SK','psi','psiL','tau_L','tau_KL','tau_SKL','psiSKL'};


%% pre-define Table 4 (all results with counterfactuals)
% ==================================================================================================
colnames = {'min_DDY0','max_DDY0','cor_DDY0_dy','sse_DDY0',...
    'min_tauScost','max_tauScost','avg_tauScost',...
    'corr_tauScost_dy','corr_tauScost_1sec','corr_SI',...
    'max_dellN','max_dyN'};


rownames = {'Data',...
    'OneSec','OneSec_counter_tauS','OneSec_counter_tauSK',...
    'Bench',...
    'Bench_counter_tau_S','Bench_counter_tau_K','Bench_counter_tau_SK',...
    'Bench_counter_psi','Bench_counter_tau_L','Bench_counter_psiL',...
    'Bench_counter_tau_KL','Bench_counter_tau_SKL','Bench_counter_psiSKL',...
    'Theta05',...
    'Theta05_counter_tau_S','Theta05_counter_tau_K','Theta05_counter_tau_SK',...
    'Theta05_counter_psi','Theta05_counter_tau_L','Theta05_counter_psiL'...
    'Theta05_counter_tau_KL','Theta05_counter_tau_SKL','Theta05_counter_psiSKL',...
    'Theta2',...
    'Theta2_counter_tau_S','Theta2_counter_tau_K','Theta2_counter_tau_SK',...
    'Theta2_counter_psi','Theta2_counter_tau_L','Theta2_counter_psiL',...
    'Theta2_counter_tau_KL','Theta2_counter_tau_SKL','Theta2_counter_psiSKL',...
    'omc05',...
    'omc05_counter_tau_S','omc05_counter_tau_K','omc05_counter_tau_SK',...
    'omc05_counter_psi','omc05_counter_tau_L','omc05_counter_psiL',...
    'omc05_counter_tau_KL','omc05_counter_tau_SKL','omc05_counter_psiSKL'};
TABLE4 = array2table(NaN(length(rownames),length(colnames)),...
    'VariableNames',colnames,'RowNames',rownames);
% ==================================================================================================



% temporary tables that make it easier to record counterfactuals
% ==================================================================================================

temp_colnames = {'DD_Y0','tauScost','SoverY','IoverY'};

temp_table_S_DD_Y0          = zeros(size(sum_stats.iso,1),1);
temp_table_K_DD_Y0          = temp_table_S_DD_Y0;
temp_table_SK_DD_Y0         = temp_table_S_DD_Y0;
temp_table_psi_DD_Y0       = temp_table_S_DD_Y0;
temp_table_psiL_DD_Y0      = temp_table_S_DD_Y0;
temp_table_L_DD_Y0          = temp_table_S_DD_Y0;
temp_table_KL_DD_Y0         = temp_table_S_DD_Y0;
temp_table_SKL_DD_Y0        = temp_table_S_DD_Y0;
temp_table_psiSKL_DD_Y0        = temp_table_S_DD_Y0;
temp_table_ORIG_DD_Y0       = temp_table_S_DD_Y0;

temp_table_S_tauScost       = temp_table_S_DD_Y0;
temp_table_K_tauScost       = temp_table_S_DD_Y0;
temp_table_SK_tauScost      = temp_table_S_DD_Y0;
temp_table_psi_tauScost    = temp_table_S_DD_Y0;
temp_table_psiL_tauScost   = temp_table_S_DD_Y0;
temp_table_L_tauScost       = temp_table_S_DD_Y0;
temp_table_KL_tauScost      = temp_table_S_DD_Y0;
temp_table_SKL_tauScost     = temp_table_S_DD_Y0;
temp_table_psiSKL_tauScost     = temp_table_S_DD_Y0;
temp_table_ORIG_tauScost    = temp_table_S_DD_Y0;

%% add tables with S/Y and I/Y
temp_table_S_SoverY       = temp_table_S_DD_Y0;
temp_table_K_SoverY       = temp_table_S_DD_Y0;
temp_table_SK_SoverY      = temp_table_S_DD_Y0;
temp_table_psi_SoverY    = temp_table_S_DD_Y0;
temp_table_psiL_SoverY   = temp_table_S_DD_Y0;
temp_table_L_SoverY       = temp_table_S_DD_Y0;
temp_table_KL_SoverY      = temp_table_S_DD_Y0;
temp_table_SKL_SoverY     = temp_table_S_DD_Y0;
temp_table_psiSKL_SoverY     = temp_table_S_DD_Y0;
temp_table_ORIG_SoverY    = temp_table_S_DD_Y0;

temp_table_S_IoverY       = temp_table_S_DD_Y0;
temp_table_K_IoverY       = temp_table_S_DD_Y0;
temp_table_SK_IoverY      = temp_table_S_DD_Y0;
temp_table_psi_IoverY     = temp_table_S_DD_Y0;
temp_table_psiL_IoverY    = temp_table_S_DD_Y0;
temp_table_L_IoverY       = temp_table_S_DD_Y0;
temp_table_KL_IoverY      = temp_table_S_DD_Y0;
temp_table_SKL_IoverY     = temp_table_S_DD_Y0;
temp_table_psiSKL_IoverY  = temp_table_S_DD_Y0;
temp_table_ORIG_IoverY    = temp_table_S_DD_Y0;

%% dy and drer
temp_table_S_dy       = temp_table_S_DD_Y0;
temp_table_K_dy       = temp_table_S_DD_Y0;
temp_table_SK_dy      = temp_table_S_DD_Y0;
temp_table_psi_dy     = temp_table_S_DD_Y0;
temp_table_psiL_dy    = temp_table_S_DD_Y0;
temp_table_L_dy       = temp_table_S_DD_Y0;
temp_table_KL_dy      = temp_table_S_DD_Y0;
temp_table_SKL_dy     = temp_table_S_DD_Y0;
temp_table_psiSKL_dy  = temp_table_S_DD_Y0;
temp_table_ORIG_dy    = temp_table_S_DD_Y0;

temp_table_S_drer       = temp_table_S_DD_Y0;
temp_table_K_drer       = temp_table_S_DD_Y0;
temp_table_SK_drer      = temp_table_S_DD_Y0;
temp_table_psi_drer     = temp_table_S_DD_Y0;
temp_table_psiL_drer    = temp_table_S_DD_Y0;
temp_table_L_drer       = temp_table_S_DD_Y0;
temp_table_KL_drer      = temp_table_S_DD_Y0;
temp_table_SKL_drer     = temp_table_S_DD_Y0;
temp_table_psiSKL_drer  = temp_table_S_DD_Y0;
temp_table_ORIG_drer    = temp_table_S_DD_Y0;


%% ellN and yNshare
temp_table_S_ellNpath      = zeros(deep_params.T,size(sum_stats.iso,1));
temp_table_K_ellNpath      = temp_table_S_ellNpath;
temp_table_SK_ellNpath     = temp_table_S_ellNpath;
temp_table_psi_ellNpath    = temp_table_S_ellNpath;
temp_table_psiL_ellNpath   = temp_table_S_ellNpath;
temp_table_L_ellNpath      = temp_table_S_ellNpath;
temp_table_KL_ellNpath     = temp_table_S_ellNpath;
temp_table_SKL_ellNpath    = temp_table_S_ellNpath;
temp_table_psiSKL_ellNpath = temp_table_S_ellNpath;
temp_table_ORIG_ellNpath   = temp_table_S_ellNpath;


temp_table_S_yNshare_path      = zeros(deep_params.T,size(sum_stats.iso,1));
temp_table_K_yNshare_path      = temp_table_S_yNshare_path;
temp_table_SK_yNshare_path     = temp_table_S_yNshare_path;
temp_table_psi_yNshare_path    = temp_table_S_yNshare_path;
temp_table_psiL_yNshare_path   = temp_table_S_yNshare_path;
temp_table_L_yNshare_path      = temp_table_S_yNshare_path;
temp_table_KL_yNshare_path     = temp_table_S_yNshare_path;
temp_table_SKL_yNshare_path    = temp_table_S_yNshare_path;
temp_table_psiSKL_yNshare_path = temp_table_S_yNshare_path;
temp_table_ORIG_yNshare_path   = temp_table_S_yNshare_path;


temp_table_S_DDelast      = zeros(size(sum_stats.iso,1),1);
temp_table_SK_DDelast     = temp_table_S_DDelast;
temp_table_SKL_DDelast    = temp_table_S_DDelast;
temp_table_psiSKL_DDelast = temp_table_S_DDelast;


temp_table_DD_Y0    = array2table(zeros(size(sum_stats,1),length(counter_types)),'VariableNames',counter_types);
temp_table_tauScost = temp_table_DD_Y0;
temp_table_SoverY = temp_table_DD_Y0;
temp_table_IoverY = temp_table_DD_Y0;
mean_tauL_vs_dy   = zeros(size(sum_stats.iso,1),1+length(cases));
mean_tauL_vs_dy(:,1+length(cases)) = sum_stats.dydata;
% ==================================================================================================
%% fill out Table 4 with data and one-sector model statistics
% ==================================================================================================

% data
TABLE4.min_DDY0('Data') = min(sum_stats.DD_Y0);
TABLE4.max_DDY0('Data') = max(sum_stats.DD_Y0);
TABLE4.cor_DDY0_dy('Data') = corr(sum_stats.DD_Y0,sum_stats.dydata);
TABLE4.corr_SI('Data')      = corr(sum_stats.IoverY,sum_stats.SoverY);

% calibrated 1-sector
TABLE4.min_DDY0('OneSec')           = min(OneSec_BENCH_DD_Y0);
TABLE4.max_DDY0('OneSec')           = max(OneSec_BENCH_DD_Y0);
TABLE4.cor_DDY0_dy('OneSec')        = corr(OneSec_BENCH_DD_Y0,OneSec_BENCH_dy);
TABLE4.sse_DDY0('OneSec')           = sum((sum_stats.DD_Y0 - OneSec_BENCH_DD_Y0).^2)/size(sum_stats,1);
TABLE4.min_tauScost('OneSec')       = min(OneSec_BENCH_tauScost);
TABLE4.max_tauScost('OneSec')       = max(OneSec_BENCH_tauScost);
TABLE4.avg_tauScost('OneSec')       = mean(OneSec_BENCH_tauScost);
TABLE4.corr_tauScost_dy('OneSec')   = corr(OneSec_BENCH_tauScost,sum_stats.dydata);
TABLE4.corr_SI('OneSec')            = corr(OneSec_BENCH_SoverY,OneSec_BENCH_IoverY);

% counterfactual 1-sector tauS = 0
TABLE4.min_DDY0('OneSec_counter_tauS')           = min(OneSec_COUNTER_tauS_DD_Y0);
TABLE4.max_DDY0('OneSec_counter_tauS')           = max(OneSec_COUNTER_tauS_DD_Y0);
TABLE4.cor_DDY0_dy('OneSec_counter_tauS')        = corr(OneSec_COUNTER_tauS_DD_Y0,OneSec_BENCH_dy);
TABLE4.sse_DDY0('OneSec_counter_tauS')           = sum((sum_stats.DD_Y0 - OneSec_COUNTER_tauS_DD_Y0).^2)/size(sum_stats,1);
TABLE4.min_tauScost('OneSec_counter_tauS')       = 0;
TABLE4.max_tauScost('OneSec_counter_tauS')       = 0;
TABLE4.avg_tauScost('OneSec_counter_tauS')       = 0;
TABLE4.corr_tauScost_dy('OneSec_counter_tauS')   = 0;
TABLE4.corr_SI('OneSec_counter_tauS')            = corr(OneSec_COUNTER_tauS_SoverY,OneSec_COUNTER_tauS_IoverY);


% counterfactual 1-sector tauS = tauK = 0
TABLE4.min_DDY0('OneSec_counter_tauSK')           = min(OneSec_COUNTER_tauSK_DD_Y0);
TABLE4.max_DDY0('OneSec_counter_tauSK')           = max(OneSec_COUNTER_tauSK_DD_Y0);
TABLE4.cor_DDY0_dy('OneSec_counter_tauSK')        = corr(OneSec_COUNTER_tauSK_DD_Y0,OneSec_BENCH_dy);
TABLE4.sse_DDY0('OneSec_counter_tauSK')           = sum((sum_stats.DD_Y0 - OneSec_COUNTER_tauSK_DD_Y0).^2)/size(sum_stats,1);
TABLE4.min_tauScost('OneSec_counter_tauSK')       = 0;
TABLE4.max_tauScost('OneSec_counter_tauSK')       = 0;
TABLE4.avg_tauScost('OneSec_counter_tauSK')       = 0;
TABLE4.corr_tauScost_dy('OneSec_counter_tauSK')   = 0;
TABLE4.corr_SI('OneSec_counter_tauSK')            = corr(OneSec_COUNTER_tauSK_SoverY,OneSec_COUNTER_tauSK_IoverY);


% ==================================================================================================

DDY0_vs_tauS_debt = zeros(globals.N_lbdagrid,size(sum_stats,1));
DDY0_vs_tauS_taus = zeros(globals.N_lbdagrid,size(sum_stats,1));

if(num_cases==4)
cases_nms = {'OneSec','Bench','Theta05','Theta2','omc05'};
elseif(num_cases==3)
cases_nms = {'OneSec','Bench','Theta05','Theta2'};
elseif(num_cases==2)
cases_nms = {'OneSec','Bench','Theta05'};
elseif(numcases==1)
cases_nms = {'OneSec','Bench'};
end

for i=1:length(cases_nms)
    DDY0_vs_tauS.(cases_nms{i}).debt = DDY0_vs_tauS_debt;
    DDY0_vs_tauS.(cases_nms{i}).taus = DDY0_vs_tauS_taus;
    if i > 1
        dDRER_vs_tauS.(cases_nms{i}).table = ...
            zeros(size(sum_stats,1),2);
    end
    
   
end

DDY0_vs_tauS.OneSec.debt = OneSec_elast_DDY0;
DDY0_vs_tauS.OneSec.taus = OneSec_elast_tauS;



disp(TABLE4)

%% pre-define table for storing cross-country data set (inputs into table 4
%  and other figures for paper.
%===================================================================================================

colnames1 = {'iso','year1','year2','dpdata','dydata','DDY0data','LN0data','LNTdata',...
    'drer','dy','DDY0','tauScost','SoverY','IoverY','LN0','LN1','LN2','LNT',...
    'avgtauL','maxdLN','yN0','yN1','yNT','maxdyN'};
rownames1 = sum_stats.iso;
ccdataset = array2table(NaN(length(rownames1),length(colnames1)),...
    'VariableNames',colnames1,'RowNames',rownames1);

cfnms = {'tau_S','tau_K','tau_SK','psi','psiL','tau_L','tau_KL','tau_SKL','psiSKL'};

for i=1:length(cases_nms)
    
    %original calibration dataset
    crosscountrydata.(cases_nms{i}).calib=ccdataset;
      %fill in data columns
    crosscountrydata.(cases_nms{i}).calib.iso=sum_stats.iso;
    crosscountrydata.(cases_nms{i}).calib.year1=sum_stats.year1;
    crosscountrydata.(cases_nms{i}).calib.year2=sum_stats.yearT;
    crosscountrydata.(cases_nms{i}).calib.dpdata=sum_stats.dpdata;
    crosscountrydata.(cases_nms{i}).calib.dydata=sum_stats.dydata;
    crosscountrydata.(cases_nms{i}).calib.DDY0data=sum_stats.DD_Y0; 
    crosscountrydata.(cases_nms{i}).calib.LN0data=sum_stats.LN0; 
    crosscountrydata.(cases_nms{i}).calib.LNTdata=sum_stats.LNT;   

    crosscountry_laborflows.(cases_nms{i}).calib=array2table(NaN(deep_params.T,length(sum_stats.iso)));

    %table for each counterfactual
    for k=1:length(cfnms)
    crosscountrydata.(cases_nms{i}).(cfnms{k})=ccdataset;
    
    crosscountry_laborflows.(cases_nms{i}).(cfnms{k})=array2table(NaN(deep_params.T,length(sum_stats.iso)));

    %fill in data columns
    crosscountrydata.(cases_nms{i}).(cfnms{k}).iso=sum_stats.iso;
    crosscountrydata.(cases_nms{i}).(cfnms{k}).year1=sum_stats.year1;
    crosscountrydata.(cases_nms{i}).(cfnms{k}).year2=sum_stats.yearT;
    crosscountrydata.(cases_nms{i}).(cfnms{k}).dpdata=sum_stats.dpdata;
    crosscountrydata.(cases_nms{i}).(cfnms{k}).dydata=sum_stats.dydata;
    crosscountrydata.(cases_nms{i}).(cfnms{k}).DDY0data=sum_stats.DD_Y0; 
    crosscountrydata.(cases_nms{i}).(cfnms{k}).LN0data=sum_stats.LN0; 
    crosscountrydata.(cases_nms{i}).(cfnms{k}).LNTdata=sum_stats.LNT; 
    end
 
    %matrices for recovered calibrated TauLgrids
        %didn't store calibrated tauLgrids for all calibrations
        recovered_tauLgrids.(cases_nms{i})=array2table(NaN(deep_params.T,length(sum_stats.iso)));
    
end

        temp_recovered_tauLgrid=zeros(deep_params.T,length(sum_stats.iso));


%===================================================================================================
%%







%%          RUN COUNTERFACTUALS
% ==================================================================================================
for i = 1:length(cases)
    
    filename = ['Calibration_Results_',cases{i},'.mat'];
    cd(mainpath_source)
    
    if i == 1
        case_name = 'Bench';
    elseif i == 2
        case_name = 'Theta05';
    elseif i == 3
        case_name = 'Theta2';
    elseif i == 4
        case_name = 'omc05';
    end
    
    
    load(filename)
    load(pathsfilename)
    %% fix some variables that could have been overwritten
    globals.N_lbdagrid = N_lbdagrid;
    
    
    cd(folders_destination.matlab)
      
    
    TABLE4.min_tauScost(case_name) = min(RESULTS.Calibrated_Table.tauScost_cal);
    TABLE4.max_tauScost(case_name) = max(RESULTS.Calibrated_Table.tauScost_cal);
    TABLE4.avg_tauScost(case_name) = mean(RESULTS.Calibrated_Table.tauScost_cal);
    
    
    %fill in crosscountry data with calibration results
    %colnames1 = {'iso','year1','year2','dydata','DDY0data','LN0data','LNTdata',...
    %'dy','DDY0','tauScost','SoverY','IoverY','SIcorr','LN0','LN1','LN2','LNT','maxdLN',...
    %'avgtauL','morevars'};
    %cfnms = {'tau_S','tau_K','tau_SK','psi','psiL','tau_L','tau_KL','tau_SKL','psiSKL'};

    crosscountrydata.(case_name).calib.drer=Calibrated_Table.drer_cal;
    crosscountrydata.(case_name).calib.dy=Calibrated_Table.dy_cal;
    crosscountrydata.(case_name).calib.DDY0=Calibrated_Table.DDY0_cal;
    crosscountrydata.(case_name).calib.tauScost=Calibrated_Table.tauScost_cal;
    crosscountrydata.(case_name).calib.SoverY=Calibrated_Table.SoverY_cal;
    crosscountrydata.(case_name).calib.IoverY=Calibrated_Table.IoverY_cal;
    crosscountrydata.(case_name).calib.LN0=sum_stats.LN0;
    crosscountrydata.(case_name).calib.LNT=sum_stats.LNT;
    
    for cntry=1:length(sum_stats.iso)
        cntrytauL=Calibrated_tauLgrid(:,cntry);
        %calibration and counterfactuals using calibrated Tau_L
        crosscountrydata.(case_name).calib.avgtauL(cntry)=mean(cntrytauL);
        crosscountrydata.(case_name).tau_S.avgtauL(cntry)=mean(cntrytauL);
        crosscountrydata.(case_name).tau_K.avgtauL(cntry)=mean(cntrytauL);
        crosscountrydata.(case_name).tau_SK.avgtauL(cntry)=mean(cntrytauL);
        crosscountrydata.(case_name).psi.avgtauL(cntry)=mean(cntrytauL);
        %counterfactuals holding Tau_L constant at initial level
        crosscountrydata.(case_name).tau_L.avgtauL(cntry)=cntrytauL(1);
        crosscountrydata.(case_name).psiL.avgtauL(cntry)=cntrytauL(1);
        crosscountrydata.(case_name).tau_KL.avgtauL(cntry)=cntrytauL(1);
        crosscountrydata.(case_name).tau_SKL.avgtauL(cntry)=cntrytauL(1);
        crosscountrydata.(case_name).psiSKL.avgtauL(cntry)=cntrytauL(1);
    end

    
    globals.sum_stats = sum_stats;
    
    globals_counter = globals;

    % setting options for parallel computing at the Bank of Canada cluster
    %setenv('MATLAB_WORKER_ARGS','-p long --time=0-96:00:00 --mem=40GB');
    %ppddy0cf2=parpool('edith92',27)      
    parfor i_country = 1:size(sum_stats,1)
        
	disp(sum_stats.iso(i_country))

       
        ci_initial_guess = Calibrated_Guesses(:,i_country);
        ci_AT            = Calibrated_AT(:,i_country);
        ci_gTgrid        = Calibrated_gTgrid(:,i_country);
        ci_AN            = Calibrated_AN(:,i_country);
        ci_ANAT          = Calibrated_ANAT(:,i_country);
        ci_tauKgrid_init = Calibrated_tauKgrid(:,i_country);
        ci_tauSgrid_init = Calibrated_tauSgrid(:,i_country);
        ci_ngrid         = initial_st_state.nbar(i_country)*ones(deep_params.T,1);
        ci_kT_0          = initial_st_state.kT_0(i_country);
        ci_kN_0          = initial_st_state.kN_0(i_country);
        ci_d_0           = initial_st_state.d_0(i_country);
        ci_ellN_path     = RESULTS.exo_paths.ellN(:,i_country);
        ci_ellT_path     = 1 - ci_ellN_path;
        ci_TT            = sum_stats.yearT(i_country) - sum_stats.year1(i_country) + 1;
        
        
        
        %% recover tau_ell
        % ----------------------------------------------------------------------------------
        
        f_rec_tauL = @(y) eqm_path(y,...
            globals,...
            ci_AT, ci_gTgrid, ci_AN, ci_ANAT, ci_tauKgrid_init, ci_tauSgrid_init,...
            ci_ngrid, ci_kT_0, ci_kN_0, ci_d_0, ci_ellN_path, ci_ellT_path);
        
        [residual,eq_paths] = f_rec_tauL(ci_initial_guess);
        
        ci_tauLgrid_init = eq_paths.tauLgrid;
        mean_tauL_vs_dy(i_country,i) = mean(ci_tauLgrid_init(1:ci_TT,1));
        % ----------------------------------------------------------------------------------
        
        temp_recovered_tauLgrid(:,i_country)=eq_paths.tauLgrid;
        
        
        %% ORIGINAL - REPLICATION OF THE CALIBRATION RESULTS -
        %  checked, it works, as long as the same eqm_path is called as the one used in calibration
        counter_type = {'ORIGINAL'};
        counterfactual_results_ORIGINAL = counterfactual(ci_TT,ci_initial_guess,ci_AT,...
            ci_gTgrid, ci_AN, ci_ANAT, ci_tauKgrid_init, ci_tauSgrid_init, ci_tauLgrid_init,...
            ci_ngrid, ci_kT_0, ci_kN_0, ci_d_0, ci_ellN_path, ci_ellT_path,...
            globals,counter_type);
        
        %% counter-factual setting tau_S = 0
        counter_type = {'tauS'};
        counterfactual_results_tau_S = counterfactual(ci_TT,ci_initial_guess,ci_AT,...
            ci_gTgrid, ci_AN, ci_ANAT, ci_tauKgrid_init, ci_tauSgrid_init, ci_tauLgrid_init,...
            ci_ngrid, ci_kT_0, ci_kN_0, ci_d_0, ci_ellN_path, ci_ellT_path,...
            globals,counter_type);
        
        
        %% counter-factual setting tau_K = 0
        counter_type = {'tauK'};
        counterfactual_results_tau_K = counterfactual(ci_TT,ci_initial_guess,ci_AT,...
            ci_gTgrid, ci_AN, ci_ANAT, ci_tauKgrid_init, ci_tauSgrid_init, ci_tauLgrid_init,...
            ci_ngrid, ci_kT_0, ci_kN_0, ci_d_0, ci_ellN_path, ci_ellT_path,...
            globals,counter_type);
        
        %% counter-factual setting tau_S = tau_K = 0
        counter_type = {'tauS','tauK'};
        counterfactual_results_tau_SK = counterfactual(ci_TT,ci_initial_guess,ci_AT,...
            ci_gTgrid, ci_AN, ci_ANAT, ci_tauKgrid_init, ci_tauSgrid_init, ci_tauLgrid_init,...
            ci_ngrid, ci_kT_0, ci_kN_0, ci_d_0, ci_ellN_path, ci_ellT_path,...
            globals,counter_type);
        
        %% counter-factual setting tau_ell = constant
        counter_type = {'tauL'};
        counterfactual_results_tau_L = counterfactual(ci_TT,ci_initial_guess,ci_AT,...
            ci_gTgrid, ci_AN, ci_ANAT, ci_tauKgrid_init, ci_tauSgrid_init, ci_tauLgrid_init,...
            ci_ngrid, ci_kT_0, ci_kN_0, ci_d_0, ci_ellN_path, ci_ellT_path,...
            globals,counter_type);
        
        %% counter-factual setting ppsi = 0
        counter_type = {'psi'};
        counterfactual_results_psi = counterfactual(ci_TT,ci_initial_guess,ci_AT,...
            ci_gTgrid, ci_AN, ci_ANAT, ci_tauKgrid_init, ci_tauSgrid_init, ci_tauLgrid_init,...
            ci_ngrid, ci_kT_0, ci_kN_0, ci_d_0, ci_ellN_path, ci_ellT_path,...
            globals,counter_type);
        
        %% counter-factual setting ppsi = 0 and tauL = constant
        counter_type = {'psi','tauL'};
        counterfactual_results_psiL = counterfactual(ci_TT,ci_initial_guess,ci_AT,...
            ci_gTgrid, ci_AN, ci_ANAT, ci_tauKgrid_init, ci_tauSgrid_init, ci_tauLgrid_init,...
            ci_ngrid, ci_kT_0, ci_kN_0, ci_d_0, ci_ellN_path, ci_ellT_path,...
            globals,counter_type);
        
        %% counter-factual setting tau_ell = constant AND tau_K = 0
        counter_type = {'tauK','tauL'};
        counterfactual_results_tau_KL = counterfactual(ci_TT,ci_initial_guess,ci_AT,...
            ci_gTgrid, ci_AN, ci_ANAT, ci_tauKgrid_init, ci_tauSgrid_init, ci_tauLgrid_init,...
            ci_ngrid, ci_kT_0, ci_kN_0, ci_d_0, ci_ellN_path, ci_ellT_path,...
            globals,counter_type);
        
        %% counter-factual setting tau_ell = constant AND tau_S = tau_K = 0
        counter_type = {'tauS','tauK','tauL'};
        counterfactual_results_tau_SKL = counterfactual(ci_TT,ci_initial_guess,ci_AT,...
            ci_gTgrid, ci_AN, ci_ANAT, ci_tauKgrid_init, ci_tauSgrid_init, ci_tauLgrid_init,...
            ci_ngrid, ci_kT_0, ci_kN_0, ci_d_0, ci_ellN_path, ci_ellT_path,...
            globals,counter_type);
        
        %% counter-factual --- ALL
        counter_type = {'tauS','tauK','tauL','psi'};
        counterfactual_results_psiSKL = counterfactual(ci_TT,ci_initial_guess,ci_AT,...
            ci_gTgrid, ci_AN, ci_ANAT, ci_tauKgrid_init, ci_tauSgrid_init, ci_tauLgrid_init,...
            ci_ngrid, ci_kT_0, ci_kN_0, ci_d_0, ci_ellN_path, ci_ellT_path,...
            globals,counter_type);
        
        %% record results for each country
        
        % Delta Debt vs tauS for different counterfactuals
        DDY0_vs_tauS_debt(:,i_country) = counterfactual_results_tau_S.DDY0_vs_tauS(:,2);
        DDY0_vs_tauS_taus(:,i_country) = counterfactual_results_tau_S.DDY0_vs_tauS(:,1);
        
        
        % \Delta Debt over initial GDP
        temp_table_S_DD_Y0(i_country)       = counterfactual_results_tau_S.DDY0;
        temp_table_K_DD_Y0(i_country)       = counterfactual_results_tau_K.DDY0;
        temp_table_SK_DD_Y0(i_country)      = counterfactual_results_tau_SK.DDY0;
        temp_table_psi_DD_Y0(i_country)     = counterfactual_results_psi.DDY0;
        temp_table_psiL_DD_Y0(i_country)    = counterfactual_results_psiL.DDY0;
        temp_table_L_DD_Y0(i_country)       = counterfactual_results_tau_L.DDY0;
        temp_table_KL_DD_Y0(i_country)      = counterfactual_results_tau_KL.DDY0;
        temp_table_SKL_DD_Y0(i_country)     = counterfactual_results_tau_SKL.DDY0;
        temp_table_psiSKL_DD_Y0(i_country)     = counterfactual_results_psiSKL.DDY0;
        temp_table_ORIG_DD_Y0(i_country)    = counterfactual_results_ORIGINAL.DDY0;
        
        % cost of distortion relative to GDP (average over the time horizon considered)
        temp_table_S_tauScost(i_country)    = counterfactual_results_tau_S.tauScost;
        temp_table_K_tauScost(i_country)    = counterfactual_results_tau_K.tauScost;
        temp_table_SK_tauScost(i_country)   = counterfactual_results_tau_SK.tauScost;
        temp_table_psi_tauScost(i_country)  = counterfactual_results_psi.tauScost;
        temp_table_psiL_tauScost(i_country) = counterfactual_results_psiL.tauScost;
        temp_table_L_tauScost(i_country)    = counterfactual_results_tau_L.tauScost;
        temp_table_KL_tauScost(i_country)   = counterfactual_results_tau_KL.tauScost;
        temp_table_SKL_tauScost(i_country)  = counterfactual_results_tau_SKL.tauScost;
        temp_table_psiSKL_tauScost(i_country)  = counterfactual_results_psiSKL.tauScost;
        temp_table_ORIG_tauScost(i_country) = counterfactual_results_ORIGINAL.tauScost;
        
        temp_table_S_SoverY(i_country)      = counterfactual_results_tau_S.SoverY;
        temp_table_K_SoverY(i_country)      = counterfactual_results_tau_K.SoverY;
        temp_table_SK_SoverY(i_country)     = counterfactual_results_tau_SK.SoverY;
        temp_table_psi_SoverY(i_country)    = counterfactual_results_psi.SoverY;
        temp_table_psiL_SoverY(i_country)   = counterfactual_results_psiL.SoverY;
        temp_table_L_SoverY(i_country)      = counterfactual_results_tau_L.SoverY;
        temp_table_KL_SoverY(i_country)     = counterfactual_results_tau_KL.SoverY;
        temp_table_SKL_SoverY(i_country)    = counterfactual_results_tau_SKL.SoverY;
        temp_table_psiSKL_SoverY(i_country)    = counterfactual_results_psiSKL.SoverY;
        temp_table_ORIG_SoverY(i_country)   = counterfactual_results_ORIGINAL.SoverY;
        
        temp_table_S_IoverY(i_country)      = counterfactual_results_tau_S.IoverY;
        temp_table_K_IoverY(i_country)      = counterfactual_results_tau_K.IoverY;
        temp_table_SK_IoverY(i_country)     = counterfactual_results_tau_SK.IoverY;
        temp_table_psi_IoverY(i_country)    = counterfactual_results_psi.IoverY;
        temp_table_psiL_IoverY(i_country)   = counterfactual_results_psiL.IoverY;
        temp_table_L_IoverY(i_country)      = counterfactual_results_tau_L.IoverY;
        temp_table_KL_IoverY(i_country)     = counterfactual_results_tau_KL.IoverY;
        temp_table_SKL_IoverY(i_country)    = counterfactual_results_tau_SKL.IoverY;
        temp_table_psiSKL_IoverY(i_country)    = counterfactual_results_psiSKL.IoverY;
        temp_table_ORIG_IoverY(i_country)   = counterfactual_results_ORIGINAL.IoverY;
        
        temp_table_S_dy(i_country)      = counterfactual_results_tau_S.dy;
        temp_table_K_dy(i_country)      = counterfactual_results_tau_K.dy;
        temp_table_SK_dy(i_country)     = counterfactual_results_tau_SK.dy;
        temp_table_psi_dy(i_country)    = counterfactual_results_psi.dy;
        temp_table_psiL_dy(i_country)   = counterfactual_results_psiL.dy;
        temp_table_L_dy(i_country)      = counterfactual_results_tau_L.dy;
        temp_table_KL_dy(i_country)     = counterfactual_results_tau_KL.dy;
        temp_table_SKL_dy(i_country)    = counterfactual_results_tau_SKL.dy;
        temp_table_psiSKL_dy(i_country)    = counterfactual_results_psiSKL.dy;
        temp_table_ORIG_dy(i_country)   = counterfactual_results_ORIGINAL.dy;
        
        temp_table_S_drer(i_country)      = counterfactual_results_tau_S.drer;
        temp_table_K_drer(i_country)      = counterfactual_results_tau_K.drer;
        temp_table_SK_drer(i_country)     = counterfactual_results_tau_SK.drer;
        temp_table_psi_drer(i_country)    = counterfactual_results_psi.drer;
        temp_table_psiL_drer(i_country)   = counterfactual_results_psiL.drer;
        temp_table_L_drer(i_country)      = counterfactual_results_tau_L.drer;
        temp_table_KL_drer(i_country)     = counterfactual_results_tau_KL.drer;
        temp_table_SKL_drer(i_country)    = counterfactual_results_tau_SKL.drer;
        temp_table_psiSKL_drer(i_country)    = counterfactual_results_psiSKL.drer;
        temp_table_ORIG_drer(i_country)   = counterfactual_results_ORIGINAL.drer;
        
        
        temp_table_S_ellNpath(:,i_country)      = counterfactual_results_tau_S.ellNpath;
        temp_table_K_ellNpath(:,i_country)      = counterfactual_results_tau_K.ellNpath;
        temp_table_SK_ellNpath(:,i_country)     = counterfactual_results_tau_SK.ellNpath;
        temp_table_psi_ellNpath(:,i_country)    = counterfactual_results_psi.ellNpath;
        temp_table_psiL_ellNpath(:,i_country)   = counterfactual_results_psiL.ellNpath;
        temp_table_L_ellNpath(:,i_country)      = counterfactual_results_tau_L.ellNpath;
        temp_table_KL_ellNpath(:,i_country)     = counterfactual_results_tau_KL.ellNpath;
        temp_table_SKL_ellNpath(:,i_country)    = counterfactual_results_tau_SKL.ellNpath;
        temp_table_psiSKL_ellNpath(:,i_country)    = counterfactual_results_psiSKL.ellNpath;
        temp_table_ORIG_ellNpath(:,i_country)   = counterfactual_results_ORIGINAL.ellNpath;

        temp_table_S_ellNpath(:,i_country)      = counterfactual_results_tau_S.ellNpath;
        temp_table_K_ellNpath(:,i_country)      = counterfactual_results_tau_K.ellNpath;
        temp_table_SK_ellNpath(:,i_country)     = counterfactual_results_tau_SK.ellNpath;
        temp_table_psi_ellNpath(:,i_country)    = counterfactual_results_psi.ellNpath;
        temp_table_psiL_ellNpath(:,i_country)   = counterfactual_results_psiL.ellNpath;
        temp_table_L_ellNpath(:,i_country)      = counterfactual_results_tau_L.ellNpath;
        temp_table_KL_ellNpath(:,i_country)     = counterfactual_results_tau_KL.ellNpath;
        temp_table_SKL_ellNpath(:,i_country)    = counterfactual_results_tau_SKL.ellNpath;
        temp_table_psiSKL_ellNpath(:,i_country)    = counterfactual_results_psiSKL.ellNpath;
        temp_table_ORIG_ellNpath(:,i_country)   = counterfactual_results_ORIGINAL.ellNpath;
        
        temp_table_S_yNshare_path(:,i_country)      = counterfactual_results_tau_S.yNshare_path;
        temp_table_K_yNshare_path(:,i_country)      = counterfactual_results_tau_K.yNshare_path;
        temp_table_SK_yNshare_path(:,i_country)     = counterfactual_results_tau_SK.yNshare_path;
        temp_table_psi_yNshare_path(:,i_country)    = counterfactual_results_psi.yNshare_path;
        temp_table_psiL_yNshare_path(:,i_country)   = counterfactual_results_psiL.yNshare_path;
        temp_table_L_yNshare_path(:,i_country)      = counterfactual_results_tau_L.yNshare_path;
        temp_table_KL_yNshare_path(:,i_country)     = counterfactual_results_tau_KL.yNshare_path;
        temp_table_SKL_yNshare_path(:,i_country)    = counterfactual_results_tau_SKL.yNshare_path;
        temp_table_psiSKL_yNshare_path(:,i_country)    = counterfactual_results_psiSKL.yNshare_path;
        temp_table_ORIG_yNshare_path(:,i_country)   = counterfactual_results_ORIGINAL.yNshare_path;
        
        
    end
    
    % delete local variable from parallel computing
    %delete(ppddy0cf2);
    % end of the loop over countries
    
    %% Record results in temporary table for ease of putting them in Table 4
    
    % #countries x 9 table with \Delta Debt / initial GDP in each counter-factual
    temp_table_DD_Y0.tau_S   = temp_table_S_DD_Y0;
    temp_table_DD_Y0.tau_K   = temp_table_K_DD_Y0;
    temp_table_DD_Y0.tau_SK  = temp_table_SK_DD_Y0;
    temp_table_DD_Y0.psi     = temp_table_psi_DD_Y0;
    temp_table_DD_Y0.psiL    = temp_table_psiL_DD_Y0;
    temp_table_DD_Y0.tau_L   = temp_table_L_DD_Y0;
    temp_table_DD_Y0.tau_KL  = temp_table_KL_DD_Y0;
    temp_table_DD_Y0.tau_SKL = temp_table_SKL_DD_Y0;
    temp_table_DD_Y0.psiSKL = temp_table_psiSKL_DD_Y0;
    
    % #countries x 9 table with dy in each counter-factual
    temp_table_dy.tau_S   = temp_table_S_dy;
    temp_table_dy.tau_K   = temp_table_K_dy;
    temp_table_dy.tau_SK  = temp_table_SK_dy;
    temp_table_dy.psi     = temp_table_psi_dy;
    temp_table_dy.psiL    = temp_table_psiL_dy;
    temp_table_dy.tau_L   = temp_table_L_dy;
    temp_table_dy.tau_KL  = temp_table_KL_dy;
    temp_table_dy.tau_SKL = temp_table_SKL_dy;
    temp_table_dy.psiSKL = temp_table_psiSKL_dy;

    % #countries x 9 table with drer in each counter-factual
    temp_table_drer.tau_S   = temp_table_S_drer;
    temp_table_drer.tau_K   = temp_table_K_drer;
    temp_table_drer.tau_SK  = temp_table_SK_drer;
    temp_table_drer.psi     = temp_table_psi_drer;
    temp_table_drer.psiL    = temp_table_psiL_drer;
    temp_table_drer.tau_L   = temp_table_L_drer;
    temp_table_drer.tau_KL  = temp_table_KL_drer;
    temp_table_drer.tau_SKL = temp_table_SKL_drer;
    temp_table_drer.psiSKL = temp_table_psiSKL_drer;    
   
    % #countries x 9 table with cost of distortion in each counter-factual (do we even need it?)
    % should we use GDP growth, \Delta RER, IoverY, and/or SoverY instead?
    temp_table_tauScost.tau_S = temp_table_S_tauScost;
    temp_table_tauScost.tau_K = temp_table_K_tauScost;
    temp_table_tauScost.tau_SK = temp_table_SK_tauScost;
    temp_table_tauScost.psi = temp_table_psi_tauScost;
    temp_table_tauScost.psiL = temp_table_psiL_tauScost;
    temp_table_tauScost.tau_L = temp_table_L_tauScost;
    temp_table_tauScost.tau_KL = temp_table_KL_tauScost;
    temp_table_tauScost.tau_SKL = temp_table_SKL_tauScost;
    temp_table_tauScost.psiSKL = temp_table_psiSKL_tauScost;
    
    temp_table_SoverY.tau_S = temp_table_S_SoverY;
    temp_table_SoverY.tau_K = temp_table_K_SoverY;
    temp_table_SoverY.tau_SK = temp_table_SK_SoverY;
    temp_table_SoverY.psi = temp_table_psi_SoverY;
    temp_table_SoverY.psiL = temp_table_psiL_SoverY;
    temp_table_SoverY.tau_L = temp_table_L_SoverY;
    temp_table_SoverY.tau_KL = temp_table_KL_SoverY;
    temp_table_SoverY.tau_SKL = temp_table_SKL_SoverY;
    temp_table_SoverY.psiSKL = temp_table_psiSKL_SoverY;
    
    temp_table_IoverY.tau_S = temp_table_S_IoverY;
    temp_table_IoverY.tau_K = temp_table_K_IoverY;
    temp_table_IoverY.tau_SK = temp_table_SK_IoverY;
    temp_table_IoverY.psi = temp_table_psi_IoverY;
    temp_table_IoverY.psiL = temp_table_psiL_IoverY;
    temp_table_IoverY.tau_L = temp_table_L_IoverY;
    temp_table_IoverY.tau_KL = temp_table_KL_IoverY;
    temp_table_IoverY.tau_SKL = temp_table_SKL_IoverY;
    temp_table_IoverY.psiSKL = temp_table_psiSKL_IoverY;

    temp_table_maxdellN.tau_S   = max(diff(temp_table_S_ellNpath))';
    temp_table_maxdellN.tau_K   = max(diff(temp_table_K_ellNpath))';
    temp_table_maxdellN.tau_SK  = max(diff(temp_table_SK_ellNpath))';
    temp_table_maxdellN.psi     = max(diff(temp_table_psi_ellNpath))';
    temp_table_maxdellN.psiL    = max(diff(temp_table_psiL_ellNpath))';
    temp_table_maxdellN.tau_L   = max(diff(temp_table_L_ellNpath))';
    temp_table_maxdellN.tau_KL  = max(diff(temp_table_KL_ellNpath))';
    temp_table_maxdellN.tau_SKL = max(diff(temp_table_SKL_ellNpath))';
    temp_table_maxdellN.psiSKL = max(diff(temp_table_psiSKL_ellNpath))';
    
    temp_table_maxdyN.tau_S   = max(diff(temp_table_S_yNshare_path))';
    temp_table_maxdyN.tau_K   = max(diff(temp_table_K_yNshare_path))';
    temp_table_maxdyN.tau_SK  = max(diff(temp_table_SK_yNshare_path))';
    temp_table_maxdyN.psi     = max(diff(temp_table_psi_yNshare_path))';
    temp_table_maxdyN.psiL    = max(diff(temp_table_psiL_yNshare_path))';
    temp_table_maxdyN.tau_L   = max(diff(temp_table_L_yNshare_path))';
    temp_table_maxdyN.tau_KL  = max(diff(temp_table_KL_yNshare_path))';
    temp_table_maxdyN.tau_SKL = max(diff(temp_table_SKL_yNshare_path))';
    temp_table_maxdyN.psiSKL = max(diff(temp_table_psiSKL_yNshare_path))';
    
    temp_table_LN0.tau_S   = transpose(temp_table_S_ellNpath(1,:));
    temp_table_LN0.tau_K   = transpose(temp_table_K_ellNpath(1,:));
    temp_table_LN0.tau_SK  = transpose(temp_table_SK_ellNpath(1,:));
    temp_table_LN0.psi     = transpose(temp_table_psi_ellNpath(1,:));
    temp_table_LN0.psiL    = transpose(temp_table_psiL_ellNpath(1,:));
    temp_table_LN0.tau_L   = transpose(temp_table_L_ellNpath(1,:));
    temp_table_LN0.tau_KL  = transpose(temp_table_KL_ellNpath(1,:));
    temp_table_LN0.tau_SKL = transpose(temp_table_SKL_ellNpath(1,:));
    temp_table_LN0.psiSKL = transpose(temp_table_psiSKL_ellNpath(1,:));

    temp_table_LN1.tau_S   = transpose(temp_table_S_ellNpath(2,:));
    temp_table_LN1.tau_K   = transpose(temp_table_K_ellNpath(2,:));
    temp_table_LN1.tau_SK  = transpose(temp_table_SK_ellNpath(2,:));
    temp_table_LN1.psi     = transpose(temp_table_psi_ellNpath(2,:));
    temp_table_LN1.psiL    = transpose(temp_table_psiL_ellNpath(2,:));
    temp_table_LN1.tau_L   = transpose(temp_table_L_ellNpath(2,:));
    temp_table_LN1.tau_KL  = transpose(temp_table_KL_ellNpath(2,:));
    temp_table_LN1.tau_SKL = transpose(temp_table_SKL_ellNpath(2,:));
    temp_table_LN1.psiSKL = transpose(temp_table_psiSKL_ellNpath(2,:));   

    temp_table_LN2.tau_S   = transpose(temp_table_S_ellNpath(3,:));
    temp_table_LN2.tau_K   = transpose(temp_table_K_ellNpath(3,:));
    temp_table_LN2.tau_SK  = transpose(temp_table_SK_ellNpath(3,:));
    temp_table_LN2.psi     = transpose(temp_table_psi_ellNpath(3,:));
    temp_table_LN2.psiL    = transpose(temp_table_psiL_ellNpath(3,:));
    temp_table_LN2.tau_L   = transpose(temp_table_L_ellNpath(3,:));
    temp_table_LN2.tau_KL  = transpose(temp_table_KL_ellNpath(3,:));
    temp_table_LN2.tau_SKL = transpose(temp_table_SKL_ellNpath(3,:));
    temp_table_LN2.psiSKL = transpose(temp_table_psiSKL_ellNpath(3,:));  

    temp_table_LNT.tau_S   = transpose(temp_table_S_ellNpath(37,:));
    temp_table_LNT.tau_K   = transpose(temp_table_K_ellNpath(37,:));
    temp_table_LNT.tau_SK  = transpose(temp_table_SK_ellNpath(37,:));
    temp_table_LNT.psi     = transpose(temp_table_psi_ellNpath(37,:));
    temp_table_LNT.psiL    = transpose(temp_table_psiL_ellNpath(37,:));
    temp_table_LNT.tau_L   = transpose(temp_table_L_ellNpath(37,:));
    temp_table_LNT.tau_KL  = transpose(temp_table_KL_ellNpath(37,:));
    temp_table_LNT.tau_SKL = transpose(temp_table_SKL_ellNpath(37,:));
    temp_table_LNT.psiSKL = transpose(temp_table_psiSKL_ellNpath(37,:));   
    
    
    temp_table_yN0.tau_S   = transpose(temp_table_S_yNshare_path(1,:));
    temp_table_yN0.tau_K   = transpose(temp_table_K_yNshare_path(1,:));
    temp_table_yN0.tau_SK  = transpose(temp_table_SK_yNshare_path(1,:));
    temp_table_yN0.psi     = transpose(temp_table_psi_yNshare_path(1,:));
    temp_table_yN0.psiL    = transpose(temp_table_psiL_yNshare_path(1,:));
    temp_table_yN0.tau_L   = transpose(temp_table_L_yNshare_path(1,:));
    temp_table_yN0.tau_KL  = transpose(temp_table_KL_yNshare_path(1,:));
    temp_table_yN0.tau_SKL = transpose(temp_table_SKL_yNshare_path(1,:));
    temp_table_yN0.psiSKL = transpose(temp_table_psiSKL_yNshare_path(1,:));
    
    temp_table_yN1.tau_S   = transpose(temp_table_S_yNshare_path(2,:));
    temp_table_yN1.tau_K   = transpose(temp_table_K_yNshare_path(2,:));
    temp_table_yN1.tau_SK  = transpose(temp_table_SK_yNshare_path(2,:));
    temp_table_yN1.psi     = transpose(temp_table_psi_yNshare_path(2,:));
    temp_table_yN1.psiL    = transpose(temp_table_psiL_yNshare_path(2,:));
    temp_table_yN1.tau_L   = transpose(temp_table_L_yNshare_path(2,:));
    temp_table_yN1.tau_KL  = transpose(temp_table_KL_yNshare_path(2,:));
    temp_table_yN1.tau_SKL = transpose(temp_table_SKL_yNshare_path(2,:));
    temp_table_yN1.psiSKL = transpose(temp_table_psiSKL_yNshare_path(2,:));    
    
    temp_table_yNT.tau_S   = transpose(temp_table_S_yNshare_path(37,:));
    temp_table_yNT.tau_K   = transpose(temp_table_K_yNshare_path(37,:));
    temp_table_yNT.tau_SK  = transpose(temp_table_SK_yNshare_path(37,:));
    temp_table_yNT.psi     = transpose(temp_table_psi_yNshare_path(37,:));
    temp_table_yNT.psiL    = transpose(temp_table_psiL_yNshare_path(37,:));
    temp_table_yNT.tau_L   = transpose(temp_table_L_yNshare_path(37,:));
    temp_table_yNT.tau_KL  = transpose(temp_table_KL_yNshare_path(37,:));
    temp_table_yNT.tau_SKL = transpose(temp_table_SKL_yNshare_path(37,:));
    temp_table_yNT.psiSKL = transpose(temp_table_psiSKL_yNshare_path(37,:));    


    temp_table_maxdellN.ORIG=max(diff(temp_table_ORIG_ellNpath))';
    temp_table_yN0.ORIG=transpose(temp_table_ORIG_yNshare_path(1,:));
    temp_table_yN1.ORIG=transpose(temp_table_ORIG_yNshare_path(2,:));
    temp_table_yNT.ORIG=transpose(temp_table_ORIG_yNshare_path(37,:));
    temp_table_maxdyN.ORIG = max(diff(temp_table_ORIG_yNshare_path))';
    

    %store the benchmark data
    crosscountrydata.Bench.calib.maxdLN=temp_table_maxdellN.ORIG;   
    crosscountrydata.Bench.calib.yN0=temp_table_yN0.ORIG;
    crosscountrydata.Bench.calib.yN1=temp_table_yN1.ORIG;
    crosscountrydata.Bench.calib.yNT=temp_table_yNT.ORIG;
    crosscountrydata.Bench.calib.maxdyN=temp_table_maxdyN.ORIG;

       
    %store into cross country data set
    %colnames1 = {'iso','year1','year2','dpdata','dydata','DDY0data','LN0data','LNTdata',...
    %'drer','dy','DDY0','tauScost','SoverY','IoverY','LN0','LN1','LN2','LNT',...
    %'avgtauL','maxdLN','yN0','yN1','yNT','maxdyN'};
    %cfnms = {'tau_S','tau_K','tau_SK','psi','psiL','tau_L','tau_KL','tau_SKL','psiSKL'}

for k=1:length(cfnms)
    crosscountrydata.(case_name).(cfnms{k}).drer=temp_table_drer.(cfnms{k});
    crosscountrydata.(case_name).(cfnms{k}).dy=temp_table_dy.(cfnms{k});
    crosscountrydata.(case_name).(cfnms{k}).DDY0=temp_table_DD_Y0.(cfnms{k});
    crosscountrydata.(case_name).(cfnms{k}).tauScost=temp_table_tauScost.(cfnms{k});
    crosscountrydata.(case_name).(cfnms{k}).SoverY=temp_table_SoverY.(cfnms{k});
    crosscountrydata.(case_name).(cfnms{k}).IoverY=temp_table_IoverY.(cfnms{k});
    
    crosscountrydata.(case_name).(cfnms{k}).LN0=temp_table_LN0.(cfnms{k});
    crosscountrydata.(case_name).(cfnms{k}).LN1=temp_table_LN1.(cfnms{k});
    crosscountrydata.(case_name).(cfnms{k}).LN2=temp_table_LN2.(cfnms{k});
    crosscountrydata.(case_name).(cfnms{k}).LNT=temp_table_LNT.(cfnms{k});
    crosscountrydata.(case_name).(cfnms{k}).maxdLN=temp_table_maxdellN.(cfnms{k});
    
    crosscountrydata.(case_name).(cfnms{k}).yN0=temp_table_yN0.(cfnms{k});
    crosscountrydata.(case_name).(cfnms{k}).yN1=temp_table_yN1.(cfnms{k});
    crosscountrydata.(case_name).(cfnms{k}).yNT=temp_table_yNT.(cfnms{k});    
    crosscountrydata.(case_name).(cfnms{k}).maxdyN=temp_table_maxdyN.(cfnms{k});

end  
   
    crosscountry_laborflows.(case_name).calib   = temp_table_ORIG_ellNpath;
    crosscountry_laborflows.(case_name).tau_S   = temp_table_S_ellNpath;
    crosscountry_laborflows.(case_name).tau_K   = temp_table_K_ellNpath;
    crosscountry_laborflows.(case_name).tau_SK  = temp_table_SK_ellNpath;
    crosscountry_laborflows.(case_name).psi     = temp_table_psi_ellNpath;
    crosscountry_laborflows.(case_name).psiL    = temp_table_psiL_ellNpath;
    crosscountry_laborflows.(case_name).tau_L   = temp_table_L_ellNpath;
    crosscountry_laborflows.(case_name).tau_KL  = temp_table_KL_ellNpath;
    crosscountry_laborflows.(case_name).tau_SKL = temp_table_SKL_ellNpath;
    crosscountry_laborflows.(case_name).psiSKL = temp_table_psiSKL_ellNpath;  

          
    recovered_tauLgrids.(case_name)=temp_recovered_tauLgrid;
  
    
    %% Record cross-country statistics in Table 4
    % ==============================================================================================
    
    % Calibrated Results (benchmark, different theta's, and flexible K)
    % ------------------------------------------------------------------------
    TABLE4.min_DDY0(case_name) = min(temp_table_ORIG_DD_Y0);
    TABLE4.max_DDY0(case_name) = max(temp_table_ORIG_DD_Y0);
    TABLE4.cor_DDY0_dy(case_name) = corr(temp_table_ORIG_DD_Y0,sum_stats.dydata);
    TABLE4.sse_DDY0(case_name) = sum((sum_stats.DD_Y0 - temp_table_ORIG_DD_Y0).^2)/size(sum_stats,1);
    
    TABLE4.min_tauScost(case_name) = min(temp_table_ORIG_tauScost);
    TABLE4.max_tauScost(case_name) = max(temp_table_ORIG_tauScost);
    TABLE4.avg_tauScost(case_name) = mean(temp_table_ORIG_tauScost);
    TABLE4.corr_tauScost_dy(case_name) = corr(temp_table_ORIG_tauScost,sum_stats.dydata);
    TABLE4.corr_SI(case_name) = corr(temp_table_ORIG_SoverY,temp_table_ORIG_IoverY);
    TABLE4.max_dellN(case_name) = max(max(diff(temp_table_ORIG_ellNpath)));
    TABLE4.max_dyN(case_name)   = max(max(diff(temp_table_ORIG_yNshare_path)));
    
    % ------------------------------------------------------------------------
    
    
    % for each case of calibrated results, loop over counter-factuals
    % ------------------------------------------------------------------------
    for ic = 1:length(counter_types)
        counter_type = counter_types{ic};
        colname = strcat('counter_',counter_types{ic});
        rowname = strcat(case_name,'_',colname);
        
        TABLE4.min_DDY0(rowname) = min(temp_table_DD_Y0.(counter_type));
        TABLE4.max_DDY0(rowname) = max(temp_table_DD_Y0.(counter_type));
        TABLE4.cor_DDY0_dy(rowname) = corr(temp_table_DD_Y0.(counter_type),sum_stats.dydata);
        TABLE4.sse_DDY0(rowname) = sum((sum_stats.DD_Y0 - temp_table_DD_Y0.(counter_type)).^2)/size(sum_stats,1);
        
        TABLE4.min_tauScost(rowname) = min(temp_table_tauScost.(counter_type));
        TABLE4.max_tauScost(rowname) = max(temp_table_tauScost.(counter_type));
        TABLE4.avg_tauScost(rowname) = mean(temp_table_tauScost.(counter_type));
        TABLE4.corr_tauScost_dy(rowname) = corr(temp_table_tauScost.(counter_type),sum_stats.dydata);
        TABLE4.corr_SI(rowname) = corr(temp_table_SoverY.(counter_type),temp_table_IoverY.(counter_type));
        
        TABLE4.max_dellN(rowname) = max(temp_table_maxdellN.(counter_type));
        TABLE4.max_dyN(rowname)   = max(temp_table_maxdyN.(counter_type));
        
    end
    % ------------------------------------------------------------------------
    
    
    
    DDY0_vs_tauS.(case_name).debt = DDY0_vs_tauS_debt;
    DDY0_vs_tauS.(case_name).taus = log(DDY0_vs_tauS_taus);
    dDRER_vs_tauS.(case_name).table(:,1) = temp_table_S_drer - temp_table_ORIG_drer;
    dDRER_vs_tauS.(case_name).table(:,2) = Calibrated_Table.tauS_cal;
    
    
    % ==============================================================================================
    
end


TABLE4.max_dellN('Data')    = TABLE4.max_dellN('Bench');
TABLE4.max_dyN('Data')      = TABLE4.max_dyN('Bench');


% Table 4 - without private flows and without median sq error (column 5)
TABLE4_PAPER = TABLE4({'Data','OneSec','OneSec_counter_tauS','OneSec_counter_tauSK',...
    'Bench','Bench_counter_tau_S','Bench_counter_tau_SK',...
    'Bench_counter_psiL','Bench_counter_psiSKL'},...
    {'min_DDY0','max_DDY0','cor_DDY0_dy','sse_DDY0','corr_SI','max_dellN'});



% Table 5 - columns 1-3 and the last column only
TABLE5_PAPER = TABLE4({'Data','Bench','Bench_counter_psi',...
    'Bench_counter_tau_L','Bench_counter_psiL'},...
    {'min_DDY0','max_DDY0','sse_DDY0','max_dellN'});






cd(mainpath_source)
filename = 'results_counterfactuals.mat';
save(filename)



disp('   Part of Table 4 - without private flows, and w/out median sq error')
disp(' ----------------------------------------------------------------------------------------- ')
disp(' ')
disp(TABLE4_PAPER)
disp(' ')
disp(' ')


disp('   Part of Table 5 -  columns 1-3 and the last column only')
disp(' ----------------------------------------------------------------------------------------- ')
disp(' ')
disp(TABLE5_PAPER)
disp(' ')
disp(' ')












