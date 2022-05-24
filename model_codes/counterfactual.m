%function counter_output = counterfactual(counterfactual_inputs,globals,counter_type)

function counter_output = counterfactual(TT,...
        initial_guess,...
        AT              ,...
        gTgrid          ,...
        AN              ,...
        ANAT            ,...
        tauKgrid_init   ,...
        tauSgrid_init   ,...
        tauLgrid_init   ,...
        ngrid           ,...
        kT_0            ,...
        kN_0            ,...
        d_0             ,...
        ellN_path       ,...
        ellT_path       ,...
    globals,counter_type)

        globals_local   = globals;
        psi_init        = globals.parameters.deep.ppsi;
if psi_init > 1
    psi_trgt = 0;
elseif psi_init < 0.001
    psi_trgt = 1000;
end

        tauSgrid        = tauSgrid_init;
        tauKgrid        = tauKgrid_init;
        tauLgrid        = tauLgrid_init;
        
        lbda_grid = linspace(0,1,globals.N_lbdagrid);
        
        
            if max(ismember(counter_type,'tauL'))
            initial_guess(1:length(ellN_path),1) = log(ellN_path./(1-ellN_path));
            end        
        
            if max(ismember(counter_type,'tauS'))                
            DDY0_vs_tauS = rand(length(lbda_grid),2);
            end
            
            
        for i_lbda = 1:length(lbda_grid)            
            
            lbda = lbda_grid(i_lbda);
            
            
            if max(ismember(counter_type,'tauS'))                
            tauSgrid = lbda*ones(length(tauSgrid_init),1) + (1-lbda)*tauSgrid_init;
            end
            
            if max(ismember(counter_type,'tauK'))
            tauKgrid = lbda*ones(length(tauKgrid_init),1) + (1-lbda)*tauKgrid_init;                
            end

            if max(ismember(counter_type,'tauL'))
            tauLgrid = lbda*tauLgrid_init(1)*ones(length(tauLgrid_init),1) + (1-lbda)*tauLgrid_init;                
            end
            
            if max(ismember(counter_type,'psi'))
            ppsi = lbda*psi_trgt + (1-lbda)*psi_init;                
            globals_local.parameters.deep.ppsi = ppsi;
            end
            
            
            if max(ismember(counter_type,'tauL'))
            f_init2 = @(y) eqm_path_flex_labor(y,globals_local,...
                AT,gTgrid,AN,ANAT,tauKgrid,tauSgrid,tauLgrid,...
                ngrid,kT_0,kN_0,d_0);            
            else
            f_init2 = @(y) eqm_path(y,globals_local,...
                AT,gTgrid,AN,ANAT,tauKgrid,tauSgrid,...
                ngrid,kT_0,kN_0,d_0,ellN_path,ellT_path);                            
            end
            
            if i_lbda == 1
            [initial_guess_,~,fflaga ]= fsolve(f_init2,initial_guess,globals.opcje_off); 
            counter_output.ini_update = norm(initial_guess - initial_guess_);
            counter_output.ini_flag = fflaga;
            initial_guess = initial_guess_;                        
            else
            [initial_guess,~,fflaga ]= fsolve(f_init2,initial_guess,globals.opcje_off);                                                    
            end
            
            if max(ismember(counter_type,'tauS'))
                [~,counter_path] = f_init2(initial_guess);
                DDY0_vs_tauS(i_lbda,1) = tauSgrid(1);
                DDY0_vs_tauS(i_lbda,2) = (counter_path.debt(TT) - counter_path.debt(1)) / counter_path.gdp(1);
                %DDY0_vs_tauS(i_lbda,2) = (counter_path.debt(TT) - counter_path.debt(1)) / counter_path.gdp(TT);
            end
            
            
        end
        counter_output.final_flag = fflaga;
        
            [~,counter_path] = f_init2(initial_guess);  
                    
            counter_output.tauScost     = mean(counter_path.tauScost(1:TT,1));                
            counter_output.dy           = mean(diff(log(counter_path.gdp(1:TT))));
            counter_output.drer         = mean(diff(log(counter_path.rer(1:TT))));
            counter_output.DDY0         = (counter_path.debt(TT) - counter_path.debt(1)) / counter_path.gdp(1);
           %counter_output.DDY0         = (counter_path.debt(TT) - counter_path.debt(1)) / counter_path.gdp(TT);
            counter_output.IoverY       = mean(counter_path.IoverY(1:TT));
            counter_output.SoverY       = mean(counter_path.SoverY(1:TT));
            counter_output.ellNpath     = counter_path.ellNpath;
            counter_output.yNshare_path = counter_path.NTshare;
        
            if max(ismember(counter_type,'tauS'))
                diff_DDY0_vs_tauS = diff(DDY0_vs_tauS);
                elast = diff_DDY0_vs_tauS(:,2)./diff_DDY0_vs_tauS(:,1);
                elast = mean(elast);
                counter_output.DDY0_vs_tauS = DDY0_vs_tauS;
                counter_output.mean_DebtTauS_elast = elast;
            end  
            
end