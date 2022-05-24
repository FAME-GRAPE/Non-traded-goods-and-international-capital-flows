function [yout,output] = main_calibration(xin,...
    data_targets,initial_guess,globals,...
    nbar,kT_0,kN_0,d_0,AN_0,TT,ellN_path,ellT_path)


T       = globals.parameters.deep.T;
gstar   = globals.parameters.deep.gstar;
Rstar   = globals.parameters.deep.Rstar;

gT = exp(xin(1));
gN = exp(xin(2));
if length(xin) > 2
    tauS = exp(xin(3));
    tauK = exp(xin(4));
else
    tauS = 1;    
    tauK = 1;
end

AT_init     = globals.Astar;
AN_init     = AN_0 * globals.Astar;
tauK_init   = ones(T,1);
tauS_init   = ones(T,1);


AT_trgt     = gT.^[0:1:TT-1]';
AN_trgt     = AN_0 * gN.^[0:1:TT-1]';
AT_trgt     = [AT_trgt;AT_trgt(length(AT_trgt),1)*gstar.^[1:1:T+1-TT]'];
AN_trgt     = [AN_trgt;AN_trgt(length(AN_trgt),1)*gstar.^[1:1:T+1-TT]'];
tauK_trgt   = tauK * ones(T,1);
tauS_trgt   = tauS * ones(T,1); tauS_trgt(TT+1:T,1) = 1;



lbda_grid = linspace(0.0,1.0,51);

x_00 = initial_guess;

for i_lbda = 1:length(lbda_grid)
    
    lbda = lbda_grid(i_lbda);
    
    tauKgrid    = lbda * tauK_trgt  + (1-lbda)*tauK_init ;
    tauSgrid    = lbda * tauS_trgt  + (1-lbda)*tauS_init ;
    
    AT          = lbda * AT_trgt    + (1-lbda)*AT_init ;
    AN          = lbda * AN_trgt    + (1-lbda)*AN_init ;
    
    ANAT    = AN./AT;
    gTgrid  = AT(2:T+1,1)./AT(1:T,1);
    

    f_eqm_path = @(y) eqm_path(y,...
        globals,...
        AT,gTgrid,AN,ANAT,tauKgrid,tauSgrid,...
        nbar*ones(T,1),kT_0,kN_0,d_0,ellN_path,ellT_path);
    
        [x_00,~,extflaga] = fsolve(f_eqm_path,x_00,globals.opcje_off);        
        
        if extflaga <= 0
            break
        end
end

if extflaga >= 1
[residual,equilibrium_path] = f_eqm_path(x_00);

model_moments = [...
    mean(diff(log(equilibrium_path.gdp(1:TT)))); ...
    mean(diff(log(equilibrium_path.rer(1:TT)))); ...
    (equilibrium_path.debt(TT) - equilibrium_path.debt(1)) / equilibrium_path.gdp(1);...
    mean(equilibrium_path.IoverY(1:TT))];

if(data_targets(2)>=-0.016) && (data_targets(2)<0.0)
%yout = model_moments - data_targets;
yout(1)=2500*(model_moments(1)-data_targets(1));
yout(2)=4000*(model_moments(2)-data_targets(2));
yout(3)=40*(model_moments(3)-data_targets(3));
yout(4)=1000*(model_moments(4)-data_targets(4));
%only use the following if using fminsearch
yout = norm(yout);
elseif(data_targets(2)<-0.016)
%yout = model_moments - data_targets;
yout(1)=3000*(model_moments(1)-data_targets(1));
yout(2)=8000*(model_moments(2)-data_targets(2));
yout(3)=40*(model_moments(3)-data_targets(3));
yout(4)=1000*(model_moments(4)-data_targets(4));
%only use the following if using fminsearch
yout = norm(yout);
else
%yout = mo
yout(1)=2500*(model_moments(1)-data_targets(1));
yout(2)=3000*(model_moments(2)-data_targets(2));
yout(3)=30*(model_moments(3)-data_targets(3));
yout(4)=1000*(model_moments(4)-data_targets(4));
%only use the following if using fminsearch
yout = norm(yout);
end



else
yout = 10000;    
end

if nargout == 2
    output.initial_guess            = x_00;
    output.calibration_vector       = model_moments - data_targets;
    output.eq_path                  = equilibrium_path;    
    output.AT                       = AT;
    output.gTgrid                   = gTgrid;
    output.AN                       = AN;
    output.ANAT                     = ANAT;
    output.tauKgrid                 = tauKgrid;
    output.tauSgrid                 = tauSgrid;
    output.gT                       = gT;
    output.gN                       = gN;
    output.tauS                     = tauS;
    output.tauK                     = tauK;
    output.residual                 = norm(residual);   % double check the equilibrium path is correctly solved
    output.taus_over_GDP_avg_debt   = mean(equilibrium_path.taus_Rdebt(1:TT,1));
    output.taus_over_GDP_avg_RkNT   = mean(equilibrium_path.taus_onetauk_q_RkN_RkT(1:TT,1));
    output.taus_over_GDP_avg        = mean(equilibrium_path.taus_total(1:TT,1));    
    output.IoverY                   = mean(equilibrium_path.IoverY(1:TT,1));
    output.SoverY                   = mean(equilibrium_path.SoverY(1:TT,1));
    output.SIcorr = corr([equilibrium_path.IoverY(1:TT,1),equilibrium_path.SoverY(1:TT,1)]);
    output.SIcorr = output.SIcorr(1,2);    
    output.dNTshare = mean(diff(equilibrium_path.NTshare(1:TT,1)));
    output.sqerr_debt    = (model_moments(3) - data_targets(3))^2;
end

end
