function [] = run_AOM_model_twobox()

% SET CONSTANTS
R13VPDB   = 0.011202;
R2VSMOW   = 155.76e-6;

% SET INITIAL CONDITIONS
d13C_DIC_initial = -15.3; % permil
dD_H2O           = -54;        % permil
d13C_CH4_initial = -37.4;      % permil
dD_CH4_initial   = -143;       % permil
DIC_initial      = 11;         % mM
CH4_initial      = 1.5;        % mM

% SET PARAMETERS
so4case   = [1 2];

% SET ISOTOPE FRACTIONATION FACTORS

% CARBON
a13eq = 0.9409;  % Net for CH4/DIC, (Mook 1974, EFF for CO2(g)-HCO3 at 50oc = 0.99542)
a13kf = 0.9623;
a13kr = a13kf.*a13eq;

% HYDROGEN
a2eq = 0.8370;
a2kf = 0.85;
a2kr = a2kf.*a2eq;

% RUN AOM MODEL
for i = 1:2
    [R2CH4dat{i},R13CH4dat{i},R13DICdat{i},tmod{i}] = ...
        AOM_ODEs_solver_twobox(so4case(i),R13VPDB,R2VSMOW,...
        d13C_DIC_initial,dD_H2O,d13C_CH4_initial,dD_CH4_initial,...
        DIC_initial,CH4_initial,a13kf,a13kr,a2kf,a2kr);
end

%%PLOT

% Load experimetal data
load('fmincon_exp_dat.mat')

% Low SO4
d13C_DICdat_l = (R13DICdat{1}./R13VPDB-1).*1000;
d13C_CH4dat_l = (R13CH4dat{1}./R13VPDB-1).*1000;
dD_CH4dat_l   = (R2CH4dat{1}./R2VSMOW-1).*1000;
% High SO4
d13C_DICdat_h = (R13DICdat{2}./R13VPDB-1).*1000;
d13C_CH4dat_h = (R13CH4dat{2}./R13VPDB-1).*1000;
dD_CH4dat_h   = (R2CH4dat{2}./R2VSMOW-1).*1000;

xl = tmod{1}; xl(1) = 0;
xh = tmod{2}; xh(1) = 0;

%%Figure for isotopes

FigFaceA     = 0.3;
xlimval      = [0 50];
AxesFontSize = 15;
TitleXpos    = -0.1;
TitleYpos    = 1.0;
FigColors    = [0 0 0; 1 0 0];
FigMarkSize  = 9;
% FigColors    = [55,126,184;77,175,74]./255;

set(figure,'Units','Centimeters','Position',[2 2 38 10])

subplot(1,3,1)
plot(xl,d13C_CH4dat_l,'Color',FigColors(1,:),'LineWidth',1); hold on
plot(xh,d13C_CH4dat_h,'Color',FigColors(2,:),'LineWidth',1); hold on
errorbar(low_so4_dat.t_exp,low_so4_dat.d13C_CH4,low_so4_dat.d13C_CH4_err,...
    'o','MarkerFaceColor','w','LineWidth',1,'CapSize',0,...
    'MarkerEdgeColor',FigColors(1,:),'Color',FigColors(1,:),'MarkerSize',FigMarkSize)
hold on
errorbar(hig_so4_dat.t_exp,hig_so4_dat.d13C_CH4,hig_so4_dat.d13C_CH4_err,...
    'o','MarkerFaceColor','w','LineWidth',1,'CapSize',0,...
    'MarkerEdgeColor',FigColors(2,:),'Color',FigColors(2,:),'MarkerSize',FigMarkSize)
xlim(xlimval)
ylim([-80 35])
xlabel('Time (days)')
ylabel(['\delta^{13}C-CH_4 (' char(8240) ')'])
set(gca,'FontSize',AxesFontSize)
t = title('A');
t.Units = 'normalized';
t.Position(1) = TitleXpos;
t.Position(2) = TitleYpos;

subplot(1,3,2)
plot(xl,dD_CH4dat_l,'Color',FigColors(1,:),'LineWidth',1); hold on
plot(xh,dD_CH4dat_h,'Color',FigColors(2,:),'LineWidth',1); hold on
errorbar(low_so4_dat.t_exp,low_so4_dat.dD_CH4,low_so4_dat.dD_CH4_err,...
    'o','MarkerFaceColor','w','LineWidth',1,'CapSize',0,...
    'MarkerEdgeColor',FigColors(1,:),'Color',FigColors(1,:),'MarkerSize',FigMarkSize)
hold on
errorbar(hig_so4_dat.t_exp,hig_so4_dat.dD_CH4,hig_so4_dat.dD_CH4_err,...
    'o','MarkerFaceColor','w','LineWidth',1,'CapSize',0,...
    'MarkerEdgeColor',FigColors(2,:),'Color',FigColors(2,:),'MarkerSize',FigMarkSize)
xlim([0 35])
ylim([-200 300])
xlabel('Time (days)')
ylabel(['\deltaD-CH_4 (' char(8240) ')'])
set(gca,'FontSize',AxesFontSize,'YTick',-200:100:200)
t = title('B');
t.Units = 'normalized';
t.Position(1) = TitleXpos;
t.Position(2) = TitleYpos;

subplot(1,3,3)
plot(xl,d13C_DICdat_l,'Color',FigColors(1,:),'LineWidth',1); hold on
plot(xh,d13C_DICdat_h,'Color',FigColors(2,:),'LineWidth',1); hold on
errorbar(low_so4_dat.t_exp,low_so4_dat.d13C_DIC,low_so4_dat.d13C_DIC_err,...
    'o','MarkerFaceColor','w','LineWidth',1,'CapSize',0,...
    'MarkerEdgeColor',FigColors(1,:),'Color',FigColors(1,:),'MarkerSize',FigMarkSize)
hold on
errorbar(hig_so4_dat.t_exp,hig_so4_dat.d13C_DIC,hig_so4_dat.d13C_DIC_err,...
    'o','MarkerFaceColor','w','LineWidth',1,'CapSize',0,...
    'MarkerEdgeColor',FigColors(2,:),'Color',FigColors(2,:),'MarkerSize',FigMarkSize)
xlim(xlimval)
xlabel('Time (days)')
ylabel(['\delta^{13}C-DIC (' char(8240) ')'])
set(gca,'FontSize',AxesFontSize,'YTick',-20:1:-16)
t = title('C');
t.Units = 'normalized';
t.Position(1) = TitleXpos;
t.Position(2) = TitleYpos;

f = get(gca,'Children');
f = flipud(f);
l = legend([f(4) f(3) f(2) f(1)],...
    'High sulfate','Low sulfate','{\itf} = 0.04','{\itf} = 0.15');
l.Location   = 'southeast';
l.FontSize   = 12;
l.NumColumns = 2;





