function [yout_,results] = calibration_1sector(xin_,stats,deep_params,globals)


xin_ = exp(xin_);
g    = xin_(1);
tauS = xin_(2);
tauK = xin_(3);

%% deep parameters
T     = deep_params.T;
TT    = stats.yearT - stats.year1 + 1;
gstar = deep_params.gstar;
Rstar = deep_params.Rstar;
dlta  = deep_params.dlta;
alfa  = deep_params.alfa;
nbar  = 1+stats.n;
Astar = globals.Astar;

if Rstar*tauK-1+dlta < 0
    yout_ = 1000*ones(length(xin_),1);
    if globals.use_fmin == 1
        yout_ = norm(yout_);
    end
    return
else
    kstar = ( alfa / (Rstar*tauK-1+dlta) )^(1/(1-alfa));
end
ystar  = kstar^alfa;
xstar  = (gstar*nbar - 1 + dlta)*kstar;
dbar   = stats.d0_over_gdp0 * ystar;
nxstar = dbar*(Rstar-gstar*nbar);
cstar  = ystar - xstar - nxstar;
if cstar < 0
    yout_ = 1000*ones(length(xin_),1);
    if globals.use_fmin == 1
        yout_ = norm(yout_);
    end
    return
end

cguess = log(cstar)*ones(T,1);



%% calibrated inputs targets
calibrated_inputs.Ngrid = nbar.^[0:1:T]';
calibrated_inputs.tauK = tauK;


tauSgrid_trgt = ones(T,1);
tauSgrid_trgt(1:TT,1) = tauS;

TFPgrid_trgt   = g.^[0:1:TT-1]';
TFPgrid_trgt   = [TFPgrid_trgt;...
    TFPgrid_trgt(length(TFPgrid_trgt),1)*gstar.^[1:1:T+1-TT]'];

if length(TFPgrid_trgt) ~= length(Astar)
    error('problem with length of TFPgrid and Astar')
end

lbda_grid = linspace(0,1,51);


for i_lbda = 1:length(lbda_grid)
    lbda = lbda_grid(i_lbda);
    calibrated_inputs.TFPgrid  = lbda*TFPgrid_trgt + (1-lbda)*Astar;
    calibrated_inputs.tauSgrid = lbda*tauSgrid_trgt + (1-lbda)*ones(T,1);
    f_eqm_path = @(y) eqm_path_1sector(y,deep_params,stats,calibrated_inputs);
    [cguess,~,flaga] = fsolve(f_eqm_path,cguess,globals.opcje_off);
end
if flaga == 1
[~,eq_paths] = f_eqm_path(cguess);
else
    yout_ = 1000*ones(length(xin_),1);
    if globals.use_fmin == 1
        yout_ = norm(yout_);
    end
    return
end

tauSgrid = calibrated_inputs.tauSgrid;

%% model moments
dy_model   = diff(log(eq_paths.gdp)); dy_model = mean(dy_model(1:TT,1));
xy_model   = mean(eq_paths.IoverY(1:TT,1));
DDY0_model = (eq_paths.DEBT(TT) - eq_paths.DEBT(1)) / eq_paths.GDP(1);
%DDY0_model = (eq_paths.DEBT(TT) - eq_paths.DEBT(1)) / eq_paths.GDP(TT);
model_mmts = [dy_model;
    xy_model;
    DDY0_model];


%% data targets
data_mmts = [   stats.dydata;
    stats.IoverY;
    stats.DD_Y0];

%% residual
yout_ = model_mmts - data_mmts;
yout_(1) = 100*yout_(1);
yout_(2) = 10*yout_(2);

if globals.use_fmin == 1
yout_ = norm(yout_);
end

if nargout == 2
results.flag_eqpath  = flaga;
results.dy           = model_mmts(1);
results.IoverY       = xy_model;
results.SoverY       = mean(eq_paths.SoverY(1:TT,1));
results.DD_Y0        = model_mmts(3);
results.dy_err       = model_mmts(1) - data_mmts(1);
results.IoverY_err   = model_mmts(2) - data_mmts(2);
results.DD_Y0_err    = model_mmts(3) - data_mmts(3);
results.g            = g;
results.tauS         = tauS;
results.tauK         = tauK;

% this is the same as in GJ
results.tauScost     = (tauSgrid(1:TT,1) - 1) .* Rstar .* (eq_paths.DEBT(1:TT,1)) ./ eq_paths.GDP(1:TT,1);
results.tauScost     = mean(abs(results.tauScost));

end

end