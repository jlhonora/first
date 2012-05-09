% SpeciesSeparation2DMAIN;
% save('Optim31brainnoiseSNR10', 'm_total', 'phi_total', 'R2_total', 'vobj', 'm_ideal', 'phi_ideal', 'R2_map', 'R2_fitted');
% clear all;

% SpeciesSeparation2DMAIN;
% save('Optim32brainnoiseSNR10', 'm_total', 'phi_total', 'R2_total', 'vobj', 'm_ideal', 'phi_ideal', 'R2_map', 'R2_fitted');
% clear all;
% 
% SpeciesSeparation2DMAIN;
% save('Optim33brainnoiseSNR10', 'm_total', 'phi_total', 'R2_total', 'vobj', 'm_ideal', 'phi_ideal', 'R2_map', 'R2_fitted');
% clear all;
% 
% SpeciesSeparation2DMAIN;
% save('Optim34brainnoiseSNR10', 'm_total', 'phi_total', 'R2_total', 'vobj', 'm_ideal', 'phi_ideal', 'R2_map', 'R2_fitted');
% clear all;
% 
% SpeciesSeparation2DMAIN;
% save('Optim35brainnoiseSNR10', 'm_total', 'phi_total', 'R2_total', 'vobj', 'm_ideal', 'phi_ideal', 'R2_map', 'R2_fitted');
% clear all;
% 
% SpeciesSeparation2DMAIN;
% save('Optim36brainnoiseSNR10', 'm_total', 'phi_total', 'R2_total', 'vobj', 'm_ideal', 'phi_ideal', 'R2_map', 'R2_fitted');
% clear all;
% 
% SpeciesSeparation2DMAIN;
% save('Optim37brainnoiseSNR10', 'm_total', 'phi_total', 'R2_total', 'vobj', 'm_ideal', 'phi_ideal', 'R2_map', 'R2_fitted');
% clear all;

% SpeciesSeparation2DMAIN;
% save('Optim38brainnoiseSNR10', 'm_total', 'phi_total', 'R2_total', 'vobj', 'm_ideal', 'phi_ideal', 'R2_map', 'R2_fitted');
% clear all;
% 
% SpeciesSeparation2DMAIN;
% save('Optim39brainnoiseSNR10', 'm_total', 'phi_total', 'R2_total', 'vobj', 'm_ideal', 'phi_ideal', 'R2_map', 'R2_fitted');
% clear all;
% 
% SpeciesSeparation2DMAIN;
% save('Optim40brainnoiseSNR10', 'm_total', 'phi_total', 'R2_total', 'vobj', 'm_ideal', 'phi_ideal', 'R2_map', 'R2_fitted');
% clear all;

% SpeciesSeparation2DMAIN;
% save('Optim41brainnoiseSNR10', 'm_total', 'phi_total', 'R2_total', 'vobj', 'm_ideal', 'phi_ideal', 'R2_map', 'R2_fitted');
% clear all;
% 
% SpeciesSeparation2DMAIN;
% save('Optim42brainnoiseSNR10', 'm_total', 'phi_total', 'R2_total', 'vobj', 'm_ideal', 'phi_ideal', 'R2_map', 'R2_fitted');
% clear all;

% SpeciesSeparation2DMAIN;
% save('Optim43brainnoiseSNR10', 'm_total', 'phi_total', 'R2_total', 'vobj', 'm_ideal', 'phi_ideal', 'R2_map', 'R2_fitted');
% clear all;
% SpeciesSeparation2DMAIN;
% save('Optim44brainnoiseSNR10', 'm_total', 'phi_total', 'R2_total', 'vobj', 'm_ideal', 'phi_ideal', 'R2_map', 'R2_fitted');
% clear all;
% SpeciesSeparation2DMAIN;
% save('Optim45brainnoiseSNR10', 'm_total', 'phi_total', 'R2_total', 'vobj', 'm_ideal', 'phi_ideal', 'R2_map', 'R2_fitted');
% clear all;
% SpeciesSeparation2DMAIN;
% save('Optim46brainnoiseSNR10', 'm_total', 'phi_total', 'R2_total', 'vobj', 'm_ideal', 'phi_ideal', 'R2_map', 'R2_fitted');
% clear all;
% SpeciesSeparation2DMAIN;
% save('Optim47brainnoiseSNR10', 'm_total', 'phi_total', 'R2_total', 'vobj', 'm_ideal', 'phi_ideal', 'R2_map', 'R2_fitted');
% clear all;

% SpeciesSeparation2DMAIN;
% save('Optim48brainnoiseSNR10', 'm_total', 'phi_total', 'R2_total', 'vobj', 'm_ideal', 'phi_ideal', 'R2_map', 'R2_fitted');
% clear all;

range_sim = 46:90;
for ii = range_sim
    str = sprintf('Optim%dphantomSNR10T2_780', ii);
    SpeciesSeparation2DMAIN;
    save(str, 'm_total', 'phi_total', 'R2_total', 'vobj', 'm_ideal', 'phi_ideal', 'R2_map', 'R2_fitted');
    clear all;
end

