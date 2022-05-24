function [yout,steady_state] = eqm_stst(xin,parameters)

% parameters either from data or calibrated
AN_0 = parameters.stst.AN_0;
dbar = parameters.stst.dbar;
LN_0 = parameters.stst.LN_0;
LT_0 = 1 - LN_0;
nbar = parameters.stst.nbar;


% deep parameters, assumed, common to all countries
gstar = parameters.deep.gstar;
Rstar = parameters.deep.Rstar;
betta = parameters.deep.betta;
sgma  = parameters.deep.sgma;
dlta  = parameters.deep.dlta;
alfaT = parameters.deep.alfaT;
alfaN = parameters.deep.alfaN;
rts   = parameters.deep.rts;
omc   = parameters.deep.omc;
omx   = parameters.deep.omx;
thta  = parameters.deep.teta;

kN_0 = exp(xin(1));
kT_0 = exp(xin(2));

% sectoral outputs
yT = kT_0.^alfaT .* LT_0.^(rts-alfaT);
yN = kN_0.^alfaN .* LN_0.^(rts-alfaN) .* (AN_0).^(1-alfaN);

% MPKs
mpk_T  = alfaT .* yT ./ kT_0;
mpk_N  = alfaN .* yN ./ kN_0;

% price of non-tradeables
qN = mpk_T/mpk_N;

% investment goods
x  = ( nbar.*gstar - 1 + dlta) * (kN_0 + kT_0);
xN = x .* (omx/(1-omx).*qN).^(-omx);
xT = (omx/(1-omx)) * xN .* qN;

% consumption
cN = max(yN - xN,0.00001);
cT = cN .* ( (omc/(1-omc))*qN ).^thta ;

% trade balance and debt
%NX = yT - cT - xT;
%residual1 = dbar * (gstar * nbar - Rstar) + NX;

NX = dbar * (Rstar - gstar*nbar);
residual1 = yT - cT - xT - NX;

qq = (xT + qN.*xN) ./ x;
%R_N = (1 - dlta).*qq + mpk_N.*qN;
R_T = (1 - dlta).*qq + mpk_T;
% -----------------------------------------------------


% Euler equation for kT
% -----------------------------------------------------------------
residual2 = qq .* gstar.^sgma - betta*R_T;
% -----------------------------------------------------------------

% % Euler equation for kN
% % -----------------------------------------------------------------
% residual3 = R_T - R_N;
% % -----------------------------------------------------------------


yout = [residual1; residual2];% residual3];

if nargout == 2
    steady_state.qN   = qN;
    steady_state.kN_0 = kN_0;
    steady_state.kT_0 = kT_0;
    steady_state.D_over_GDP  = dbar  / (qN*yN + yT);
    steady_state.yN_over_GDP = qN*yN / (qN*yN + yT);
end

end