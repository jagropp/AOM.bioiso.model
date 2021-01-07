clear 

xe(1,1:4) = [0.01 0.1201 0.8699 0.2];
xe(2,1:4) = [0.01 0.9800 0.0100 0.5];
xe(1:2,5:10) = repmat([0.979 0.9847 0.8881 1 1 0.7089],2,1);
xe(1:2,5:10) = repmat([0.979 0.9847 0.8881 1 1 0.7089],2,1);

d13C_CH4 = zeros(2,12);
d13C_DIC = zeros(2,12);
dD_CH4   = zeros(2,12);

for i = 1:2
    [d13C_CH4(i,:),d13C_DIC(i,:),dD_CH4(i,:)] = ...
                            run_bioiso_model(xe(i,:),i);
end

load fmincon_exp_dat.mat

fig = [1 0];
plot_AOM_results
