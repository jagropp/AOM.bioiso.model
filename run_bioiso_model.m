function [d13C_CH4,d13C_DIC,dD_CH4] = ...
                        run_bioiso_model(x,so4case)

% SET CONSTANTS
R13VPDB   = 0.011202;
R2VSMOW   = 155.76e-6;

% SET INITIAL CONDITIONS
d13C_DIC_initial = [-15.4 -15.4]; % permil
dD_H2O           = -54;           % permil
d13C_CH4_initial = -37.4;         % permil
dD_CH4_initial   = -143;          % permil
DIC_initial      = 10.5;          % mM
CH4_initial      = 1.5;           % mM
sulfate_initial  = [1.15 10.8];   % mM
sulfide_initial  = [0.2 0.3];     % mM

% SET PARAMETERS
fdG = x(1:4);

% SET ISOTOPE FRACTIONATION FACTORS

% CARBON
% (1) CH4<->CH3-SCoM; (2) CH3-SCoM<->CHO-MFR; (3) CHO-MFR<->DIC
a13eqnet = 0.9409;  % Net for CH4/DIC, (Mook 1974, EFF for CO2(g)-HCO3 at 50oc = 0.99542)
a13eq(1) = 0.99922;
a13eq(2) = 0.9540; 
a13eq(3) = (a13eqnet/(a13eq(1)*a13eq(2)));

a13kf(1) = 0.9623;
a13kf(2) = x(5);
a13kf(3) = x(6);

% HYDROGEN
% (1) CH4<->H2O
% (2) CH3-SCOM<->H2O 
% (3) CH4<->CH3-SCOM
% (4) HSCoB<->H2O
% (5) CH4<->HSCoB
% (6) CH3-SCOM<->CHOMFR
% (7) CHO-MFR<->H2O

a2eq(1) = 0.8370;
a2eq(3) = 0.9568;
a2eq(4) = 0.4686;
a2eq(6) = 0.9764;
a2eq(2) = a2eq(1)/a2eq(3);
a2eq(5) = a2eq(1)/a2eq(4);
a2eq(7) = a2eq(1)/a2eq(3)/a2eq(6);

a2kf(1) = NaN;
a2kf(2) = x(7);
a2kf(3) = 0.8547; % Scheller 2013
a2kf(4) = x(8);
a2kf(5) = 0.4098; % Scheller 2013
a2kf(6) = x(9);
a2kf(7) = x(10);

a13kr   = a13kf.*a13eq;
a2kr    = a2kf.*a2eq;

% RUN AOM MODEL
[R2CH4_dat,R13CH4_dat,R13DIC_dat,t_mod] = ...
    AOM_ODEs_solver(so4case,R13VPDB,R2VSMOW,...
    d13C_DIC_initial(so4case),dD_H2O,d13C_CH4_initial,dD_CH4_initial,...
    DIC_initial,CH4_initial,sulfate_initial(so4case),...
    sulfide_initial(so4case),a13kf,a13kr,a2kf,a2kr,fdG);

days_values = [3 6 10 14 17 20 24 27 32 40 47];

day_position = zeros(1,length(days_values));
for day = 1:length(days_values)
    day_position(day) = find(t_mod>=days_values(day),1,'first');
end
day_position = [2 day_position];

d13C_CH4 = (R13CH4_dat(day_position)./R13VPDB-1).*1000;
d13C_DIC = (R13DIC_dat(day_position)./R13VPDB-1).*1000;
dD_CH4   = (R2CH4_dat(day_position)./R2VSMOW-1).*1000;

