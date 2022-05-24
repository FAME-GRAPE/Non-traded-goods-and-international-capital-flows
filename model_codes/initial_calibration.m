function [yout,globals] = initial_calibration(xxin,ii,globals)

parameters = globals.parameters;

% country-specific parameters, need to be always reset
% ------------------------------------------------------------
parameters.stst.AN_0  = exp( xxin(1) );
parameters.stst.dbar  =      xxin(2);
parameters.stst.LN_0  = globals.stats.LN0(ii);
parameters.stst.nbar = 1 + globals.stats.n(ii);
% ------------------------------------------------------------


f_stst = @(y) eqm_stst(y,parameters);

x0 = -0.5*ones(2,1);
[x0,fval,flaga] = fsolve(f_stst,x0,globals.opcje_off);
% if norm(fval) > 0.0001
%     disp(['OOOOPS,  ---  ',num2str(ii)])    
% end

[~,stst] = f_stst(x0);


yout = [stst.D_over_GDP  - globals.stats.d0_over_gdp0(ii);
    stst.yN_over_GDP - globals.stats.pN0yN0_over_gdp0(ii)];

if strcmp(globals.solution_method,'fmin') == 1
    yout = norm(yout);
end

if nargout == 2
    globals.calibration_stst_output.fval_eqm(ii) = norm(fval);
    globals.calibration_stst_output.flag_eqm(ii) = flaga;
    globals.calibration_stst_output.AN_0(ii) = parameters.stst.AN_0;
    globals.calibration_stst_output.dbar(ii) = parameters.stst.dbar;
    globals.calibration_stst_output.LN_0(ii) = parameters.stst.LN_0;
    globals.calibration_stst_output.nbar(ii) = parameters.stst.nbar;
    globals.calibration_stst_output.qN(ii)   = stst.qN;
    globals.calibration_stst_output.kT_0(ii) = stst.kT_0;
    globals.calibration_stst_output.kN_0(ii) = stst.kN_0;
    globals.calibration_stst_output.D0_Y0_data(ii) = globals.stats.d0_over_gdp0(ii);
    globals.calibration_stst_output.D0_Y0_model(ii) = stst.D_over_GDP;
    globals.calibration_stst_output.NT_GDP_data(ii) = globals.stats.pN0yN0_over_gdp0(ii);
    globals.calibration_stst_output.NT_GDP_model(ii) = stst.yN_over_GDP;    
end

end