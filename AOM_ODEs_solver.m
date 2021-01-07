function [R2CH4dat,R13CH4dat,R13DICdat,tmod] = ...
    AOM_ODEs_solver(so4case,R13VPDB,R2VSMOW,d13C_DIC_initial,dD_H2O,...
    d13C_CH4_initial,dD_CH4_initial,DIC_initial,CH4_initial,...
    sulfate_initial,sulfide_initial,a13kf,a13kr,a2kf,a2kr,fdG)

% CALCULATE CH4 AND DIC ISOTOPE COMPOSITION EVOLUTION WITH TIME

CH4     = CH4_initial;
DIC     = DIC_initial;
sulfate = sulfate_initial;
sulfide = sulfide_initial;
R13CH4  = (d13C_CH4_initial/1000 + 1)*R13VPDB;
R2CH4   = (dD_CH4_initial/1000 + 1)*R2VSMOW;
R13DIC  = (d13C_DIC_initial/1000 + 1)*R13VPDB;

t         = 0;
tmax      = 50; % Model run time in days
tmod_i    = NaN;

while t < tmax
    X = [DIC(end),R13DIC(end),CH4(end),R13CH4(end),R2CH4(end)];
    % Calculate the total dGr of the pathway from activities
    aH2S = 0.55*sulfide(end);
    aDIC = 0.55*DIC(end);
    aCH4 = 1*CH4(end);
    aSO4 = 0.11*sulfate(end);
    dGr0 = -35.8e3;
    R    = 8.31446;
    Tk   = 273.15 + 50;
    dGrtot = (dGr0 + R.*Tk.*log(aDIC.*aH2S./aCH4./aSO4))./1000;
    
    dXdt = CO2CH4Eqs(X,a13kf,a13kr,a2kf,...
        a2kr,dD_H2O,R2VSMOW,so4case,CH4_initial,dGrtot,fdG);
    fracChange = X./dXdt;
    dt = min(1e-3.*abs(fracChange(X>0)));
    
    X = X + dXdt.*dt;
    t = t + dt;
    
    DIC     = cat(1,DIC,X(1));
    R13DIC  = cat(1,R13DIC,X(2));
    CH4     = cat(1,CH4,X(3));
    R13CH4  = cat(1,R13CH4,X(4));
    R2CH4   = cat(1,R2CH4,X(5));
    sulfate = cat(1,sulfate,sulfate(end) + dXdt(3)*dt);
    sulfide = cat(1,sulfide,sulfide(end) + dXdt(1)*dt);
    
    fCH4 = min(X(3)/CH4(1));
    if fCH4 <= 0.001 || sulfate(end) <= 0.001 || dGrtot >= 0
        t = tmax;
    end
    tmod_i    = [tmod_i t];
end

d13C_DIC = (R13DIC./R13VPDB - 1).*1000;
d13C_CH4 = (R13CH4./R13VPDB - 1).*1000;
dD_CH4   = (R2CH4./R2VSMOW - 1).*1000;

cond = ~isnan(d13C_DIC) & ~isinf(d13C_DIC) & ...
    ~isnan(d13C_CH4) & ~isinf(d13C_CH4) & ...
    ~isnan(dD_CH4) & ~isinf(dD_CH4);

R2CH4dat   = R2CH4(cond)';
R13CH4dat  = R13CH4(cond)';
R13DICdat  = R13DIC(cond)';
tmod       = tmod_i;

function dXdt = CO2CH4Eqs(X,a13kf_vals,a13kr_vals,...
                           a2kf_vals,a2kr_vals,d2H_H2O,R2VSMOW,...
                           so4case,CH4i,dGr_tot,fdG)
dXdt = zeros(size(X));

R2H2O = ((d2H_H2O./1000)+1).*R2VSMOW;

DIC     = X(1);
R13DIC  = X(2);
CH4     = X(3);
R13CH4  = X(4);
R2CH4   = X(5);

fCH4 = CH4/CH4i;
switch so4case
    case 1 % Low SO42=
        Jnet = 0.1507*fCH4 - 0.0349;
    case 2 % High SO42=
        Jnet = 0.1998*fCH4;
end

dGu      = 0.253*dGr_tot; % The dGr for reactions 1-3
dGr_res  = dGr_tot-(-5)-dGu;
dGr(1:3) = fdG(1:3).*dGu;
dGr(4)   = fdG(4).*dGr_res;

R  = 8.31446e-3;
Tk = 50 + 273.15; % K
f  = exp(dGr./R./Tk);
Jf = Jnet./(1-f);
Jr = Jf - Jnet;

% Set KFFs
a13kf = [a13kf_vals(1) a13kf_vals(2) a13kf_vals(3)];
a13kr = [a13kr_vals(1) a13kr_vals(2) a13kr_vals(3)];
a2kfp = [a2kf_vals(1) a2kf_vals(2) a2kf_vals(4) a2kf_vals(5) a2kf_vals(7)];
a2krp = [a2kr_vals(1) a2kr_vals(2) a2kr_vals(4) a2kr_vals(5) a2kr_vals(7)];
a2kfs = [a2kf_vals(3) a2kf_vals(6)];
a2krs = [a2kr_vals(3) a2kr_vals(6)];

% Steady solution for 13CHO-MFR
R13CHO = (Jf(2)*a13kf(2)*Jf(1)*a13kf(1)*R13CH4/(Jr(1)*a13kr(1)+Jf(2)*a13kf(2))+Jr(3)*a13kr(3)*R13DIC)/...
    (Jr(2)*a13kr(2)+Jf(3)*a13kf(3)-Jr(2)*a13kr(2)*Jf(2)*a13kf(2)/(Jr(1)*a13kr(1)+Jf(2)*a13kf(2)));
% Steady solution for 13CH3-S-CoM
R13CH3 = (Jf(1)*a13kf(1)*R13CH4+Jr(2)*a13kr(2)*R13CHO)/(Jr(1)*a13kr(1)+Jf(2)*a13kf(2));          

% Steady solution for CH2D-S-CoM        
R2CH3 = (3*R2CH4*a2kfs(1)*Jf(1)+R2H2O*(2*a2krp(2)*Jr(2)+...
    (a2krs(2)*Jr(2)*a2krp(5)*Jr(3)/(a2krs(2)*Jr(2)+a2kfp(5)*Jf(3)))))/...
    (3*a2krs(1)*Jr(1)+a2kfs(2)*Jf(2)+2*a2kfp(2)*Jf(2)-...
    (a2kfs(2)*Jf(2)*a2krs(2)*Jr(2)/(a2krs(2)*Jr(2)+a2kfp(5)*Jf(3))));

R2COB = (Jr(4)*R2H2O*a2krp(3) + Jf(1)*R2CH4*a2kfp(4))./...
              (Jr(1)*a2krp(4) + Jf(4)*a2kfp(3));

dXdt(1) = Jf(3) - Jr(3);
dXdt(2) = 1/DIC*(Jf(3)*R13CHO*a13kf(3) - Jr(3)*R13DIC*a13kr(3) - R13DIC*dXdt(1));

% Equation for CH4
dXdt(3) = Jr(1) - Jf(1);

% Equation for R13CH4
dXdt(4) = 1/CH4*(Jr(1)*R13CH3*a13kr(1) - Jf(1)*R13CH4*a13kf(1) - R13CH4*dXdt(3));

% Equation for R2CH4
dXdt(5) = 1/CH4*(Jr(1)*(3/4*R2CH3*a2krs(1) + 1/4*R2COB*a2krp(4)) - ...
                 Jf(1)*R2CH4*(3/4*a2kfs(1) + 1/4*a2kfp(4)) - ...
                 R2CH4*dXdt(3));

