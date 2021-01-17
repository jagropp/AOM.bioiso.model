function [] = run_AOM_model()

% Best fit model parameters obtained with fmincon (tables S5 and S6)

% dGr allocation between reactions 1-4
xe(1,1:4) = [0.01 0.1201 0.8699 0.2];
xe(2,1:4) = [0.01 0.9800 0.0100 0.5];
% Kinetic fractionation factors
xe(1:2,5:10) = repmat([0.979 0.9847 0.8881 1 1 0.7089],2,1);

d13C_CH4 = zeros(2,12);
d13C_DIC = zeros(2,12);
dD_CH4   = zeros(2,12);

for i = 1:2
    [d13C_CH4(i,:),d13C_DIC(i,:),dD_CH4(i,:)] = ...
                            run_bioiso_model(xe(i,:),i);
end

% Load data from experiments
load('fmincon_exp_dat.mat','hig_so4_dat','low_so4_dat')

plot_AOM_results
