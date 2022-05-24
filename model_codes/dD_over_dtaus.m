%{

Replication file - sensitivity
"Non-Traded Goods, Factor Markets Frictions, and International Capital Flows"
by Jacek Rothert (USNA), and Jacob Short (BoC)
Review of Economic Dynamics, 2022 (accepted)
(c) JR&JS, September 2021


INPUTS: 
- OneSec.mat                                    - ( output of JRJS_NT_counterfactuals.m )
- results_counterfactuals.mat                   - ( output of JRJS_NT_counterfactuals.m )

OUTPUTS: 
- prints Table6 in the MATLAB Command Window

EXPECTED RUN TIME: 1-2 seconds

%}

clear all
clc

% inputs
filename = 'OneSec.mat';
load(filename)
filename = 'results_counterfactuals.mat';
load(filename)

Ncountry = size(DDY0_vs_tauS.OneSec.taus,2);
Ntauvals = size(DDY0_vs_tauS.OneSec.taus,1);

for i = 1:Ncountry
    if abs(DDY0_vs_tauS.OneSec.taus(1,i)) < abs(DDY0_vs_tauS.OneSec.taus(Ntauvals,i))
        DDY0_vs_tauS.OneSec.taus(:,i) = sort(DDY0_vs_tauS.OneSec.taus(:,i),'descend');
        temp = DDY0_vs_tauS.OneSec.debt(:,i);
        for t = 1:Ntauvals
            DDY0_vs_tauS.OneSec.debt(t,i) = temp(Ntauvals-t+1,1);
        end
    end
end

cases = {'Bench','Theta05','Theta2','omc05','OneSec'};
for ic = 1:length(cases)
    case_name = cases{ic};
    taus = DDY0_vs_tauS.(case_name).taus;
    debt = DDY0_vs_tauS.(case_name).debt;
    dtaus = diff(taus);
    dD    = diff(debt);
    dD    = dD          ./ abs(debt(1:Ntauvals-1,:)); % dtaus is already in \% terms
    dD_dtaus.(case_name) = -dD(1 ,:) ./ dtaus(1 ,:);  % first row means evaluating at the empirical value
    dD_dtaus.(case_name) = dD_dtaus.(case_name)';     % transpose
end


TABLE6 = array2table(zeros(1,4),'VariableNames',cases(:,1:4));

for ic = 1:length(cases)-1
    case_name = cases{ic};
    TABLE6.(case_name) = median(dD_dtaus.(case_name)) ./ median(dD_dtaus.OneSec);
end


disp(' TABLE6 - Median elasticity of cumulative inflows w.r.t. tau_S (one-sector = 1)  ') 
disp(' ------------------------------------------------------------------------------- ')
disp(' ')
disp(TABLE6)
