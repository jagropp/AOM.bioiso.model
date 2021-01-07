
function [R2CH4dat,R13CH4dat,R13DICdat,tmod] = ...
    AOM_ODEs_solver_twobox(so4case,R13VPDB,R2VSMOW,d13C_DIC_initial,dD_H2O,...
    d13C_CH4_initial,dD_CH4_initial,DIC_initial,CH4_initial,...
    a13kf,a13kr,a2kf,a2kr)

% CALCULATE CH4 AND DIC ISOTOPE COMPOSITION EVOLUTION WITH TIME
CH4     = CH4_initial;
DIC     = DIC_initial;
R13CH4  = (d13C_CH4_initial/1000 + 1)*R13VPDB;
R2CH4   = (dD_CH4_initial/1000 + 1)*R2VSMOW;
R13DIC  = (d13C_DIC_initial/1000 + 1)*R13VPDB;

t         = 0;
tmax      = 50; % Model run time in days
tmod_i    = NaN;

while t < tmax
    X = [DIC(end),R13DIC(end),CH4(end),R13CH4(end),R2CH4(end)];
    
    dXdt = CO2CH4Eqs(X,a13kf,a13kr,a2kf,a2kr,dD_H2O,R2VSMOW,so4case,CH4_initial);
    fracChange = X./dXdt;
    dt = min(5e-3.*abs(fracChange(X>0)));
    
    X = X + dXdt.*dt;
    t = t + dt;
    
    DIC     = cat(1,DIC,X(1));
    R13DIC  = cat(1,R13DIC,X(2));
    CH4     = cat(1,CH4,X(3));
    R13CH4  = cat(1,R13CH4,X(4));
    R2CH4   = cat(1,R2CH4,X(5));
    
    fCH4 = min(X(3)/CH4(1));
    if fCH4 <= 0.01
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

function dXdt = CO2CH4Eqs(X,a13kf,a13kr,a2kf,a2kr,d2H_H2O,R2VSMOW,so4case,CH4i)
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

switch so4case
    case 1
        f = 0.15;
    case 2
        f = 0.04;
end

Jf = Jnet/(1-f(1));
Jr = Jf - Jnet;

dXdt(1) = Jf - Jr;
dXdt(2) = 1/DIC*(Jf*R13CH4*a13kf - Jr*R13DIC*a13kr - R13DIC*dXdt(1));

% Equation for CH4
dXdt(3) = Jr - Jf;

% Equation for R13CH4
dXdt(4) = 1/CH4*(Jr*R13DIC*a13kr - Jf*R13CH4*a13kf - R13CH4*dXdt(3));
dXdt(5) = 1/CH4*(Jr*R2H2O*a2kr - Jf*R2CH4*a2kf - R2CH4*dXdt(3));

