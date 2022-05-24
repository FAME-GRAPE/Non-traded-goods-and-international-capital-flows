function [yout,eq_path]= eqm_path_flex_labor(xin,...
    globals,...
    AT,gTgrid,AN,ANAT,tauKgrid,tauSgrid,tauLgrid,...
    ngrid,kT_0,kN_0,dbar)

T       = length(xin) / 3;

% deep parameters, assumed, common to all countries
gstar = globals.parameters.deep.gstar;
Rstar = globals.parameters.deep.Rstar;
betta = globals.parameters.deep.betta;
sgma  = globals.parameters.deep.sgma;
dlta  = globals.parameters.deep.dlta;
ppsi  = globals.parameters.deep.ppsi;
alfaT = globals.parameters.deep.alfaT;
alfaN = globals.parameters.deep.alfaN;
rts   = globals.parameters.deep.rts;
omc   = globals.parameters.deep.omc;
omx   = globals.parameters.deep.omx;
thta  = globals.parameters.deep.teta;

yout = ones(length(xin),1);

ellN =  exp(xin(    1:1*T,1)) ./ ( 1 + exp(xin(    1:1*T,1)) );
kN   = [kN_0;exp(xin(  T+1:2*T,1))];
kT   = [kT_0;exp(xin(2*T+1:3*T,1))];

ellT = 1 - ellN;

% pre-define some vectors
% -----------------------------------------------------------------

% investment expenditures
iN = kN(2:T+1,:).*( ngrid.*gTgrid ) - (1-dlta)*kN(1:T,1); %+ adjN
iT = kT(2:T+1,:).*( ngrid.*gTgrid ) - (1-dlta)*kT(1:T,1); %+ adjT

% non-negativity investment penalty function derivative
pnltyN = ppsi*min(zeros(T,1),iN).^2;
pnltyT = ppsi*min(zeros(T,1),iT).^2;

% sectoral outputs
yT = kT(1:T,1).^alfaT .* ellT.^(rts-alfaT);
yN = kN(1:T,1).^alfaN .* ellN.^(rts-alfaN) .* (AN(1:T,1)./AT(1:T,1)).^(1-alfaN);


%% price of non-tradeables
mplT = (rts-alfaT)*yT./ellT;
mplN = (rts-alfaN)*yN./ellN;
qN = mplT ./ mplN ./ tauLgrid;

% investment goods
x  = max(iT + iN,0.000001);
xN = x .* (omx/(1-omx).*qN).^(-omx);
xT = (omx/(1-omx)) * xN .* qN;

% consumption
cN = max(yN - xN,0.00001);
cT = cN .* ( (omc/(1-omc))*qN ).^thta ;
if abs(thta-1) < 0.01
    c = cT.^omc .* cN.^(1-omc);
else
    c = ( omc.*cT.^(1-1/thta) + (1-omc).*cN.^(1-1/thta) ).^(thta/(thta-1));
end

% trade balance and debt
NX = yT - cT - xT;
debt = dbar*ones(T+1,1);
ca = zeros(T,1);
for t = 2:T+1
    debt(t,1) = (Rstar*debt(t-1,1) - NX(t-1,1))/gTgrid(t-1,1)/ngrid(t-1,1);
    ca(t-1,1) = NX(t-1,1) - (Rstar-1)*debt(t-1,1);
end

% partials of consumption and investment aggregates
% -----------------------------------------------------
U_c  = c.^(-sgma);
G_T  =   omc   * (c./cT).^(1/thta);
G_N  = (1-omc) * (c./cN).^(1/thta);
H_T  =   omx   * (x./xT);
H_N  = (1-omx) * (x./xN);
% -----------------------------------------------------

% short-cuts
% -----------------------------------------------------
Uct  = U_c(1:T-1,1);
GTt  = G_T(1:T-1,1);
HTt  = H_T(1:T-1,1);
GNt  = G_N(1:T-1,1);
HNt  = H_N(1:T-1,1);

Uct1  = U_c(2:T,1);
GTt1  = G_T(2:T,1);
HTt1  = H_T(2:T,1);
GNt1  = G_N(2:T,1);
HNt1  = H_N(2:T,1);

zTt  = pnltyT(1:T-1,1);
zTt1 = betta * pnltyT(2:T,1);
zNt  = pnltyN(1:T-1,1);
zNt1 = betta * pnltyN(2:T,1);

mpk_T  = alfaT .* yT(1:T,1) ./ kT(1:T,1);
mpk_N  = alfaN .* yN(1:T,1) ./ kN(1:T,1);

qq = (xT + qN.*xN) ./ x;
R_N = (1 - dlta).*qq + mpk_N.*qN;
R_T = (1 - dlta).*qq + mpk_T;

% -----------------------------------------------------

% Euler equation for debt
% -----------------------------------------------------------------
LHS2 = Uct .* (gTgrid(1:T-1,1)).^sgma .* GTt;
RHS2 = betta * Uct1 .* GTt1 .*Rstar .* tauSgrid(1:T-1,1);
yout(1:T-1,1) = LHS2 - RHS2;
yout(  T  ,1) = debt(T+1,1) - debt(T,1);
% -----------------------------------------------------------------

% Euler equation for kT
% -----------------------------------------------------------------
LHS3 =        Uct .* (gTgrid(1:T-1,1)).^sgma .* (GTt./ HTt) .* ...
    (1 - zTt./(Uct.*(GTt./ HTt)));
RHS3 = betta* Uct1 .* (GTt1./HTt1).*( mpk_T(2:T,1).*HTt1 + (1-dlta) .* ...
    (1 - zTt1./(Uct1.*(GTt1./ HTt1))) );
yout(T+1:2*T-1,1)  = LHS3 - RHS3;
yout(     2*T  ,1) =  kT(T+1,1) - kT(T,1);

LHS3_ = GTt .* Uct .* qq(1:T-1,1) .* (gTgrid(1:T-1,1)).^sgma;
RHS3_ = betta* Uct1 .* GTt1 .* R_T(2:T,1) .* tauKgrid(1:T-1,1)  ...
    + ( gTgrid(1:T-1,1).*pnltyT(1:T-1,1) - betta*(1-dlta).*pnltyT(2:T,1) );
yout(T+1:2*T-1,1)  = LHS3_ - RHS3_;
% -----------------------------------------------------------------

% Euler equation for kN
% -----------------------------------------------------------------
LHS4 =        Uct .* (gTgrid(1:T-1,1)).^sgma .* (GNt./ HNt) .* ...
    (1 - zNt./(Uct.*(GNt./ HNt)));
RHS4 = betta* Uct1 .* (GNt1./HNt1).*( mpk_N(2:T,1).*HNt1 + (1-dlta) .* ...
    (1 - zNt1./(Uct1.*(GNt1./ HNt1))) );
yout(2*T+1:3*T-1,1) = LHS4 - RHS4;
yout(      3*T  ,1) = kN(T+1,1) - kN(T,1);

LHS4_ = GTt .* Uct .* qq(1:T-1) .* (gTgrid(1:T-1,1)).^sgma;
RHS4_ = betta* Uct1 .* GTt1 .* R_N(2:T,1) .* tauKgrid(1:T-1,1)  ...
    + ( gTgrid(1:T-1,1).*pnltyN(1:T-1,1) - betta*(1-dlta).*pnltyN(2:T,1) );
yout(2*T+1:3*T-1,1) = LHS4_ - RHS4_;
% -----------------------------------------------------------------


if nargout == 2
    colnames = {'debt','gdp','rer','IoverY'};
    eq_path = array2table(zeros(T,size(colnames,2)),'VariableNames',colnames);
    eq_path.debt        = debt(1:T,1).*AT(1:T,1);
    eq_path.gdp         = (qN(1).*yN + yT).*AT(1:T,1);
    eq_path.gdp_nominal = (qN.*yN + yT).*AT(1:T,1);
    eq_path.NTshare     = qN.*yN ./ (qN.*yN + yT);
    eq_path.rer         = (qN.*cN + cT)./c;
    eq_path.SoverY      = (qN.*yN + yT - qN.*cN - cT)./(qN.*yN + yT);
    eq_path.IoverY      = (qN.*xN + xT)./(qN.*yN + yT);
    eq_path.taus_Rdebt  = Rstar * abs((tauSgrid-1).* eq_path.debt(1:T)) ./ eq_path.gdp_nominal;
    eq_path.taus_onetauk_q_RkN_RkT  = ...
        (tauSgrid(1:T,1)-1).*tauKgrid(1:T,1) ...
        .* qq(1:T,1) ...
        .* ( R_N(1:T,1).*kN(1:T,1).*AT(1:T,1)  +  R_T(1:T,1) .* kT(1:T,1).*AT(1:T,1) );
    eq_path.taus_onetauk_q_RkN_RkT  = abs(eq_path.taus_onetauk_q_RkN_RkT./eq_path.gdp_nominal);
    eq_path.taus_total  = eq_path.taus_Rdebt;
    
    eq_path.tauScost    = (tauSgrid(1:T,1)-1) .* Rstar .* eq_path.debt(1:T) ;
    eq_path.tauScost = abs(eq_path.tauScost)./eq_path.gdp_nominal;
    eq_path.ellNpath = ellN;
    
end

end