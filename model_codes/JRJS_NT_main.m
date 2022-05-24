%{

Replication file - calibration
"Non-Traded Goods, Factor Markets Frictions, and International Capital Flows"
by Jacek Rothert (USNA), and Jacob Short (BoC)
Review of Economic Dynamics, 2022 (accepted)
(c) JR&JS, September 2021


INPUTS: 
- JRJS_updateddata.xlsx

OUTPUTS: 
- Calibration_Results_teta01_ppsi1000.mat           
- Calibration_Results_teta05_ppsi1000.mat           
- Calibration_Results_teta2_ppsi1000.mat            
- Calibration_Results_teta01_ppsi1000_omc05.mat     
NOTE: only one output per single run of this script; need to adjust key parameters on lines 71-73


DEPENDENCIES:
a) initial_calibration.m
	(a-i) calls eqm_stst.m
b) eqm_path.m
c) main_calibration.m
	(c-i) calls eqm_path.m

EXPECTED RUN TIME: XYZ hours, (parallel computing, using XYZ nodes)

%}

%%

clear all
mainpath = 'C:\Users\R2D2\Dropbox (UW)\research\JRJS_CapitalFlows\submissions\RED\Data_Computing_Codes\model_codes';
% option to have separate folders for excel, matlab, and results files
ExcelFolder     = mainpath;     
MatlabFolder    = mainpath;     
ResultsFolder   = mainpath;     
folders.excel   = ExcelFolder;
folders.matlab  = MatlabFolder;
folders.results = ResultsFolder;


% set names for input files
% -------------------------------------------------------------------------
InputFile = 'JRJS_updateddata.xlsx';   % input file with business cycle statistics
DataSheet = 'onlyLN';
ParamSheet= 'parameters';
% -------------------------------------------------------------------------




% read input files with deep parameters and countries' stats
% -------------------------------------------------------------------------
cd(ExcelFolder)
sum_stats   = readtable(InputFile,'Sheet',DataSheet);  % country level stats
deep_params = readtable(InputFile,'Sheet',ParamSheet); % deep parameters
sum_stats.Properties.RowNames = sum_stats.iso;
model_varnames = {'iso','sse','gT','gN','dy','dp',...
    'DD_Y0','LN0','LNT','dLN',...
    'SoverY','IoverY','si_diff','si_level'};
% -------------------------------------------------------------------------




%%
%% ________________________ ADJUST KEY PARAMETERS ______________________
deep_params.teta = 0.1;  % elasticity of substitution
deep_params.ppsi = 1000.00;  % penalty on non-negative investment
deep_params.omc  = 0.20;  % penalty on non-negative investment


% benchmark - ppsi = 1000 (putty-clay, i.e. non-negative gross investment at the sectoral level)
% alternative - ppsi = 0 (putty-putty, i.e. capital easily moved between sectors)

filename = ['Calibration_Results_teta',strrep(num2str(deep_params.teta),'.',''),...
    '_ppsi',strrep(num2str(deep_params.ppsi),'.',''),'.mat'];
if abs(deep_params.omc - 0.5) < 0.001
filename = ['Calibration_Results_teta',strrep(num2str(deep_params.teta),'.',''),...
    '_ppsi',strrep(num2str(deep_params.ppsi),'_omc05.',''),'.mat'];    
end

oldresultsfile='oldresults.mat';

% --------------------------------------------------------------------------


deep_params.T   = 70;       % T = 70 seems to be sufficient (transition + 40/50) 
deep_params.rts = 0.97;     % decreasing return to scale

opcje_iter= optimset('Display','iter','MaxIter',3); 
% ONLY WHEN TESTING (change in the parfor loop in fminsearch('main_calibration'))

opcje_off = optimset('Display','off');

%options for fsolve routine on outer calibration
opfsolve = optimset('Algorithm','levenberg-marquardt','Display','off','MaxFunEvals',15000,'MaxIter',15000);

%==0 then do not read in initial guess for calibration, otherwise read in guess from previous calibration
readinitguess_cal=0;

%%


% world TFP frontier
% ---------------------------------------------------------
Astar = deep_params.gstar.^[0:1:deep_params.T]';
globals.Astar = Astar;
% ---------------------------------------------------------

parameters.deep = deep_params;


% pre-define initial calibration parameters
% -------------------------------------------------------------------------
parameters.calibration.AT_grid      = Astar;
parameters.calibration.AN_grid      = Astar;
parameters.calibration.gT_grid      = Astar(2:length(Astar),1)./Astar(1:length(Astar)-1,1);
parameters.calibration.gN_grid      = Astar(2:length(Astar),1)./Astar(1:length(Astar)-1,1);
parameters.calibration.tauK_grid    = ones(deep_params.T-1,1);
parameters.calibration.tauS_grid    = ones(deep_params.T-1,1);
% -------------------------------------------------------------------------


% calibrate initial steady-states
% -------------------------------------------------------------------------
globals.opcje_off  = opcje_off;
globals.opcje_iter = opcje_iter;
globals.opfsolve   = opfsolve;
globals.stats      = sum_stats;
globals.parameters = parameters;
colnames = {'iso','res_ini_path','update_ini_path','fval_eqm','flag_eqm','fval_inical','flag_inical',...
    'AN_0','dbar','LN_0','nbar','qN','kT_0','kN_0','D0_Y0_data','D0_Y0_model','NT_GDP_data','NT_GDP_model'};
globals.calibration_stst_output = array2table(zeros(size(sum_stats,1),size(colnames,2)),'VariableNames',colnames,...
    'RowNames',sum_stats.iso);
globals.calibration_stst_output.iso = sum_stats.iso;
cd(MatlabFolder)


initial_guesses.qN_path = zeros(deep_params.T,size(sum_stats,1));
initial_guesses.kN_path = zeros(deep_params.T,size(sum_stats,1));
initial_guesses.kT_path = zeros(deep_params.T,size(sum_stats,1));
initial_guesses.LN_path = zeros(deep_params.T,size(sum_stats,1));


stst_colnames = {'kN_0','kT_0','AN_0','LN_0','d_0','nbar'};
initial_st_state = array2table(zeros(size(sum_stats,1),length(stst_colnames)),'VariableNames',stst_colnames);


colnames = {'iso','match_dellN_flag','gT_cal','gN_cal','tauS_cal','tauK_cal','sse_cal','tauScost_cal','SIcorr_cal','SoverY_cal','IoverY_cal','sse_debt_cal','dNTshare_cal'};
Calibrated_Table =  array2table(zeros(size(sum_stats,1),length(colnames)),'VariableNames',colnames);
Calibrated_Table.iso = sum_stats.iso;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      INITIAL STEADY-STATE AND EQUILIBRIUM PATH WITH MATCHED REALLOCATION OF LABOR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for i_country = 1:size(sum_stats,1)
    
    % country and initial guess for AN0 and debt
    country = sum_stats.iso{i_country};
    
    disp([ ' steady-state calibration + labor reallocation for    -------    ', country ])
    
    x_00 = zeros(2,1);
    
    globals.eqm_paths.(country).initial_calibration = array2table(zeros(deep_params.T,size(colnames,2)),'VariableNames',colnames);
    
    %% calibration of AN0 and initial debt
    %  together with initial kN_0 and kT_0
    % -------------------------------------------------------------------------------
    globals.solution_method = 'fmin';
    f_init = @(y) initial_calibration(y,i_country,globals);
    [x_00,fmin_val,fmin_flag] = fminsearch(f_init,x_00,globals.opcje_off);
    
    globals.solution_method = 'fsolve';
    f_init = @(y) initial_calibration(y,i_country,globals);
    [x01,fvalue,flagga] = fsolve(f_init,x_00,globals.opcje_off);
    
    [residual,globals] = f_init(x01);
    
    globals.calibration_stst_output.fval_inical(i_country) = norm(fvalue);
    globals.calibration_stst_output.flag_inical(i_country) = flagga;
    
    initial_st_state.kT_0(i_country) = globals.calibration_stst_output.kT_0(i_country);
    initial_st_state.kN_0(i_country) = globals.calibration_stst_output.kN_0(i_country);
    initial_st_state.d_0(i_country)  = globals.calibration_stst_output.dbar(i_country);
    initial_st_state.AN_0(i_country) = globals.calibration_stst_output.AN_0(i_country);
    initial_st_state.LN_0(i_country) = globals.calibration_stst_output.LN_0(i_country);
    initial_st_state.nbar(i_country) = globals.calibration_stst_output.nbar(i_country);
    
    % -------------------------------------------------------------------------------
    
    
    
    
    %% check that the initial steady-state solve the eqm path along the BGP
    % -------------------------------------------------------------------------------
    
    ellN_path = globals.calibration_stst_output.LN_0(i_country)*ones(deep_params.T,1);
    ellT_path = 1 - ellN_path;
    kT_0 = globals.calibration_stst_output.kT_0(i_country);
    kN_0 = globals.calibration_stst_output.kN_0(i_country);
    d_0  = globals.calibration_stst_output.dbar(i_country);
    ngrid = globals.calibration_stst_output.nbar(i_country)*ones(deep_params.T,1);
    
    AN = Astar*globals.calibration_stst_output.AN_0(i_country);
    AT = Astar;     ANAT = AN./AT;
    gTgrid = AT(2:deep_params.T+1,1) ./ AT(1:deep_params.T,1);
    
    initial_guess = [...
        log(globals.calibration_stst_output.qN(i_country))*ones(deep_params.T,1);   ...
        log(globals.calibration_stst_output.kN_0(i_country))*ones(deep_params.T,1); ...
        log(globals.calibration_stst_output.kT_0(i_country))*ones(deep_params.T,1)  ];
    
    
    
    eqpath_resid = eqm_path(initial_guess,globals,...
        AT,gTgrid,AN,ANAT,ones(deep_params.T,1),ones(deep_params.T,1),...
        ngrid,kT_0,kN_0,d_0,ellN_path,ellT_path);
    globals.calibration_stst_output.res_ini_path(i_country) = norm(eqpath_resid);
    
    
    f_init2 = @(y) eqm_path(y,globals,...
        AT,gTgrid,AN,ANAT,ones(deep_params.T,1),ones(deep_params.T,1),...
        ngrid,kT_0,kN_0,d_0,ellN_path,ellT_path);
    [initial_guess_,fff,fflaga ]= fsolve(f_init2,initial_guess,globals.opcje_off);
    if fflaga ~= 1 % display if there is an issue with a particular country
        disp(['eqm path issue   -----  ', num2str(i_country)])
    end
    globals.calibration_stst_output.update_ini_path(i_country) = norm(initial_guess_ - initial_guess);
    if globals.calibration_stst_output.update_ini_path(i_country) > 0.001
        disp(['initial path diff from st-state in country   -----  ', num2str(i_country)])
    end
    % -------------------------------------------------------------------------------
    
    
    
    
    %% calibrate the initial eqm path with labor reallocation
    % -------------------------------------------------------------------------------
    TT      = globals.stats.yearT(i_country) - globals.stats.year1(i_country);
    d_ell   = (globals.stats.LNT(i_country) - globals.stats.LN0(i_country))/TT;
    
    ellN_initial = globals.calibration_stst_output.LN_0(i_country)*ones(deep_params.T,1);
    ellN_target  = ellN_initial;
    for t = 2:TT+1
        ellN_target(t,1) = ellN_target(t-1,1) + d_ell;
    end
    for t = TT+1:deep_params.T-1
        ellN_target(t+1,1) = ellN_target(t,1);
    end
    
    lbda_grid = linspace(0,1,101);
    i_lbda = 1;
    for i_lbda = 1:length(lbda_grid)
        ellN_path = lbda_grid(i_lbda) * ellN_target + ...
            (1-lbda_grid(i_lbda)) * ellN_initial;
        ellT_path = 1 - ellN_path;
        f_eqm_path = @(y) eqm_path(y,globals,...
            AT,gTgrid,AN,ANAT,ones(deep_params.T,1),ones(deep_params.T,1),...
            ngrid,kT_0,kN_0,d_0,ellN_path,ellT_path);
        [initial_guess,fff,fflaga ]= fsolve(f_eqm_path,initial_guess,globals.opcje_off);
        if fflaga ~= 1 % display if there is an issue with a particular country
            disp(['reallocation issue   -----  ', num2str(i_country)])
        end
    end
    
    Calibrated_Table.match_dellN_flag(i_country) = fflaga;
    
    qN = initial_guess(                1:1*deep_params.T,1);
    kN = initial_guess(  deep_params.T+1:2*deep_params.T,1);
    kT = initial_guess(2*deep_params.T+1:3*deep_params.T,1);
    
    initial_guesses.qN_path(:,i_country) = qN;
    initial_guesses.kN_path(:,i_country) = kN;
    initial_guesses.kT_path(:,i_country) = kT;
    initial_guesses.LN_path(:,i_country) = ellN_path;
    
    globals.eqm_paths.(country).initial_calibration.qN = exp(qN);
    globals.eqm_paths.(country).initial_calibration.kN = exp(kN);
    globals.eqm_paths.(country).initial_calibration.kT = exp(kT);
    
    
    % -------------------------------------------------------------------------------
    
    
    
end

gT_cal          = zeros(size(sum_stats,1),1);
gN_cal          = zeros(size(sum_stats,1),1);
tauS_cal        = zeros(size(sum_stats,1),1);
tauK_cal        = zeros(size(sum_stats,1),1);
sse_cal         = zeros(size(sum_stats,1),1);
tauScost_cal    = zeros(size(sum_stats,1),1);
SIcorr_cal      = zeros(size(sum_stats,1),1);
SoverY_cal      = zeros(size(sum_stats,1),1);
IoverY_cal      = zeros(size(sum_stats,1),1);
sse_debt_cal    = zeros(size(sum_stats,1),1);
dNTshare_cal    = zeros(size(sum_stats,1),1);

dy_cal          = zeros(size(sum_stats,1),1);
drer_cal        = zeros(size(sum_stats,1),1);
DDY0_cal        = zeros(size(sum_stats,1),1);

fval_cal        = zeros(size(sum_stats,1),1);
fexitflag_cal        = zeros(size(sum_stats,1),1);

Calibrated_debt         = zeros(deep_params.T,size(sum_stats,1));
Calibrated_gdp          = zeros(deep_params.T,size(sum_stats,1));
Calibrated_rer          = zeros(deep_params.T,size(sum_stats,1));
Calibrated_tauLgrid     = zeros(deep_params.T,size(sum_stats,1));
%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                 M A I N    C A L I B R A T I O N
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('  ')
disp(' ---- calibration running ---- ')


gstar = deep_params.gstar;

Calibrated_Guesses = NaN(...
    length(...
        [   initial_guesses.qN_path(:,1);...
            initial_guesses.kN_path(:,1);...
            initial_guesses.kT_path(:,1)        ]...
            ),...
    size(sum_stats,1)...
    );

%% OPTIONAL - read in guesses for gT, gN, tauS, tauK
% ----------------------------------------------
if readinitguess_cal == 1
cd(ResultsFolder)
load(oldresultsfile,'RESULTS')
cd(MatlabFolder)
gTgNtaus_guess = [...
    log(RESULTS.Calibrated_Table.gT_cal),...
    log(RESULTS.Calibrated_Table.gN_cal),...
    log(RESULTS.Calibrated_Table.tauS_cal),...
    log(RESULTS.Calibrated_Table.tauK_cal)...
    ];
end
% ----------------------------------------------
%%


%% the two lines below set options for parallel computing on the Bank of Canada cluster
%setenv('MATLAB_WORKER_ARGS','-p long --time=0-96:00:00 --mem=16GB');
%ppbenchv2=parpool('edithsecond',36)
%%

%% START MAIN LOOP OVER COUNTRIES
parfor i_country = 1:size(sum_stats,1)
    
    country = sum_stats.iso{i_country};
    
    ini_guess = [...
        initial_guesses.qN_path(:,i_country);...
        initial_guesses.kN_path(:,i_country);...
        initial_guesses.kT_path(:,i_country)];
    
    TT = sum_stats.yearT(i_country) - sum_stats.year1(i_country)+1;
    
    ellN_path   = initial_guesses.LN_path(:,i_country);
    ellT_path   = 1-ellN_path;
    kT_0        = initial_st_state.kT_0(i_country);
    kN_0        = initial_st_state.kN_0(i_country);
    d_0         = initial_st_state.d_0(i_country);
    AN_0        = initial_st_state.AN_0(i_country);
    LN_0        = initial_st_state.LN_0(i_country);
    nbar        = initial_st_state.nbar(i_country);
    
    data_targets = [sum_stats.dydata(i_country);...
        sum_stats.dpdata(i_country);...
        sum_stats.DD_Y0(i_country);...
        sum_stats.IoverY(i_country)];
    
    f_main_cal = @(y) main_calibration(y,data_targets,ini_guess,globals,nbar,kT_0,kN_0,d_0,AN_0,TT,ellN_path,ellT_path);
    
    if readinitguess_cal==0
        %initial guess of : g_y, g_RER, taus, tauk
        y_cal_guess = [log(gstar);log(gstar);0.05;0.03];
    else
        y_cal_guess = gTgNtaus_guess(i_country,:)';
    end
   
    
    
    
    disp([ ' -------- MAIN CALIBRATION ------------- ', country ])
    
    %call fminsearh (nelder-mead)
    [y_cal_guess,calfval,calexitflag] = fminsearch(f_main_cal,y_cal_guess,globals.opcje_off);
    %call fsolve (levenberg-marquardt algorithm)
    %[y_cal_guess,calfval,calexitflag] = fsolve(f_main_cal,y_cal_guess,globals.opfsolve);
    
    [res,results] = f_main_cal(y_cal_guess);
    
    
    % RECORD RESULTS
    gT_cal(i_country)       = exp(y_cal_guess(1));
    gN_cal(i_country)       = exp(y_cal_guess(2));
    tauS_cal(i_country)     = exp(y_cal_guess(3));
    tauK_cal(i_country)     = exp(y_cal_guess(4));
    sse_cal(i_country)      = norm(res);
    tauScost_cal(i_country) = results.taus_over_GDP_avg;
    SIcorr_cal(i_country)   = results.SIcorr;
    SoverY_cal(i_country)   = results.SoverY;
    IoverY_cal(i_country)   = results.IoverY;
    sse_debt_cal(i_country) = results.sqerr_debt;
    dNTshare_cal(i_country) = results.dNTshare;
        
    dy_cal(i_country)       = results.calibration_vector(1) + data_targets(1);
    drer_cal(i_country)     = results.calibration_vector(2) + data_targets(2);
    DDY0_cal(i_country)     = results.calibration_vector(3) + data_targets(3);

    fval_cal(i_country)      = calfval;
    %change to the following if using fsolve
    %fval_cal(i_country)      = norm(calfval);
    fexitflag_cal(i_country) = calexitflag;
    
    Calibrated_debt(:,i_country)        = results.eq_path.debt;
    Calibrated_gdp(:,i_country)         = results.eq_path.gdp;
    Calibrated_rer(:,i_country)         = results.eq_path.rer;
    Calibrated_tauLgrid(:,i_country)    = results.eq_path.tauLgrid;    
    
    
    Calibrated_Guesses(:,i_country)     = results.initial_guess; % this is the updated guess
    Calibrated_AT(:,i_country)          = results.AT;
    Calibrated_gTgrid(:,i_country)      = results.gTgrid;
    Calibrated_AN(:,i_country)          = results.AN;
    Calibrated_ANAT(:,i_country)        = results.ANAT;
    Calibrated_tauKgrid(:,i_country)    = results.tauKgrid;
    Calibrated_tauSgrid(:,i_country)    = results.tauSgrid;
    
    
end

% the line below is from computation on the Bank of Canada cluster
%delete(ppbenchv2);
%%

Calibrated_Table.gT_cal       = gT_cal;
Calibrated_Table.gN_cal       = gN_cal;
Calibrated_Table.tauS_cal     = tauS_cal;
Calibrated_Table.tauK_cal     = tauK_cal;
Calibrated_Table.sse_cal      = sse_cal;
Calibrated_Table.tauScost_cal = tauScost_cal;
Calibrated_Table.SIcorr_cal   = SIcorr_cal;
Calibrated_Table.SoverY_cal   = SoverY_cal;
Calibrated_Table.IoverY_cal   = IoverY_cal;
Calibrated_Table.sse_debt_cal = sse_debt_cal;
Calibrated_Table.dNTshare_cal = dNTshare_cal;

Calibrated_Table.dy_cal   = dy_cal;
Calibrated_Table.drer_cal = drer_cal;
Calibrated_Table.DDY0_cal = DDY0_cal;
Calibrated_Table.fval_cal = fval_cal;
Calibrated_Table.fexitflag_cal = fexitflag_cal;    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RESULTS.Calibrated_Table        = Calibrated_Table;
RESULTS.Calibrated_Guesses      = Calibrated_Guesses;
RESULTS.Initial_Steady_States   = initial_st_state;
RESULTS.exo_paths.ellN          = initial_guesses.LN_path;

RESULTS.exo_paths.AT        = Calibrated_AT;
RESULTS.exo_paths.gTgrid    = Calibrated_gTgrid;
RESULTS.exo_paths.AN        = Calibrated_AN;
RESULTS.exo_paths.ANAT      = Calibrated_ANAT;
RESULTS.exo_paths.tauKgrid  = Calibrated_tauKgrid;
RESULTS.exo_paths.tauSgrid  = Calibrated_tauSgrid;

RESULTS.endo_paths.Calibrated_debt = Calibrated_debt;
RESULTS.endo_paths.Calibrated_gdp = Calibrated_gdp;
RESULTS.endo_paths.Calibrated_rer = Calibrated_rer;
RESULTS.endo_paths.Calibrated_tauLgrid = Calibrated_tauLgrid;
%%

cd(ResultsFolder)
save(filename)




disp('  ')
disp(' ---- calibration finished ---- ')
