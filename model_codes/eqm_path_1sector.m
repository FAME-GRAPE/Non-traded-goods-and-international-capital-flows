function [yout,eq_path] = eqm_path_1sector(xin,deep_params,data_inputs,calibrated_inputs)

%% parameters
T     = deep_params.T;
alfa  = deep_params.alfa;
dlta  = deep_params.dlta;
sgma  = deep_params.sgma;
betta = deep_params.betta;
Rstar = deep_params.Rstar;

tauK        = calibrated_inputs.tauK;
tauSgrid    = calibrated_inputs.tauSgrid;
A           = calibrated_inputs.TFPgrid;
N           = calibrated_inputs.Ngrid;

if length(A) ~= T+1
    error('lengths of TFP and T are inconsistent')
end

if length(xin) ~= T
    error('lengths of xin and T are inconsistent')
end
ggrid = A(2:T+1,1)./A(1:T,1);
ngrid = N(2:T+1,1)./N(1:T,1);

kstar = ( alfa  / (Rstar*tauK-1+dlta) )^(1/(1-alfa));
kpath = ( alfa ./ (Rstar*tauK./tauSgrid-1+dlta) ).^(1/(1-alfa));
k           = [kstar;kpath];
x           = ngrid.*ggrid.*k(2:T+1,1) - (1-dlta)*k(1:T,1);
y           = k.^alfa;
d0          = data_inputs.d0_over_gdp0 * y(1);
debt        = d0*ones(length(k),1);

c  = exp(xin);
nx = y(1:T,1) - x - c;

for t = 1:T
    debt(t+1) = (debt(t)*Rstar - nx(t))/ggrid(t)/ngrid(t);
end

yout = [...
    c(1:T-1,1).^(-sgma) .* (ggrid(1:T-1,1)).^sgma - betta*Rstar*tauSgrid(1:T-1,1).*c(2:T,1).^(-sgma);...
    debt(T+1) - debt(T)...
    ];

if nargout == 2
    eq_path.y    = y;
    eq_path.c    = c;
    eq_path.x    = x;    
    eq_path.k    = k;    
    eq_path.nx   = nx;
    eq_path.debt = debt;
    eq_path.IoverY =                x ./ y(1:T,1);
    eq_path.SoverY = ( y(1:T,1) - c ) ./ y(1:T,1);
    eq_path.gdp  =  y        .* A;
    eq_path.GDP  =  y   .* N .* A;
    eq_path.DEBT = debt .* N .* A;    
    eq_path.K    =  k   .* N .* A;    
end
end