x   = low_so4_dat.t_exp;

x_lim_val      = [0 50];
axes_font_size = 15;
fig_colors     = [0 0 0; 1 0 0];
fig_mark_size  = 9;
line_width     = 1.5;

set(figure,'Units','Centimeters','Position',[5 5 42 9])

subplot(1,4,1)
for i = 1:2
    plot(x,d13C_CH4(i,:),'Color',fig_colors(i,:),'LineWidth',line_width)
    hold on
end
errorbar(low_so4_dat.t_exp,low_so4_dat.d13C_CH4,low_so4_dat.d13C_CH4_err,...
    'o','MarkerFaceColor','w','LineWidth',1,'CapSize',0,...
    'MarkerEdgeColor',fig_colors(1,:),'Color',fig_colors(1,:),'MarkerSize',fig_mark_size)
hold on
errorbar(hig_so4_dat.t_exp,hig_so4_dat.d13C_CH4,hig_so4_dat.d13C_CH4_err,...
    'o','MarkerFaceColor','w','LineWidth',1,'CapSize',0,...
    'MarkerEdgeColor',fig_colors(2,:),'Color',fig_colors(2,:),'MarkerSize',fig_mark_size)
xlim(x_lim_val)
ylim([-80 35])
xlabel('Time (days)')
ylabel(['\delta^{13}C-CH_4 (' char(8240) ')'])
set(gca,'FontSize',axes_font_size)

subplot(1,4,2)
for i = 1:2
    plot(x,dD_CH4(i,:),'Color',fig_colors(i,:),'LineWidth',line_width)
    hold on
end
errorbar(low_so4_dat.t_exp,low_so4_dat.dD_CH4,low_so4_dat.dD_CH4_err,...
    'o','MarkerFaceColor','w','LineWidth',1,'CapSize',0,...
    'MarkerEdgeColor',fig_colors(1,:),'Color',fig_colors(1,:),'MarkerSize',fig_mark_size)
hold on
errorbar(hig_so4_dat.t_exp,hig_so4_dat.dD_CH4,hig_so4_dat.dD_CH4_err,...
    'o','MarkerFaceColor','w','LineWidth',1,'CapSize',0,...
    'MarkerEdgeColor',fig_colors(2,:),'Color',fig_colors(2,:),'MarkerSize',fig_mark_size)
xlim([0 35])
ylim([-220 480])
xlabel('Time (days)')
ylabel(['\deltaD-CH_4 (' char(8240) ')'])
set(gca,'FontSize',axes_font_size,'YTick',-200:200:400)

subplot(1,4,3)
for i = 1:2
    plot(x,d13C_DIC(i,:),'Color',fig_colors(i,:),'LineWidth',line_width)
    hold on
end
errorbar(low_so4_dat.t_exp,low_so4_dat.d13C_DIC,low_so4_dat.d13C_DIC_err,...
    'o','MarkerFaceColor','w','LineWidth',1,'CapSize',0,...
    'MarkerEdgeColor',fig_colors(1,:),'Color',fig_colors(1,:),'MarkerSize',fig_mark_size)
hold on
errorbar(hig_so4_dat.t_exp,hig_so4_dat.d13C_DIC,hig_so4_dat.d13C_DIC_err,...
    'o','MarkerFaceColor','w','LineWidth',1,'CapSize',0,...
    'MarkerEdgeColor',fig_colors(2,:),'Color',fig_colors(2,:),'MarkerSize',fig_mark_size)
xlim(x_lim_val)
xlabel('Time (days)')
ylabel(['\delta^{13}C-DIC (' char(8240) ')'])
set(gca,'FontSize',axes_font_size,'YTick',[-20 -18 -16])

subplot(1,4,4)
for i = 1:2
    hold on
    plot(dD_CH4(i,:),d13C_CH4(i,:),'Color',fig_colors(i,:),'LineWidth',line_width)
    hold on
end
errorbar(low_so4_dat.dD_CH4',low_so4_dat.d13C_CH4',...
    low_so4_dat.d13C_CH4_err,low_so4_dat.d13C_CH4_err,...
    low_so4_dat.dD_CH4_err,low_so4_dat.dD_CH4_err,...
    'o','MarkerFaceColor','w','LineWidth',1,'CapSize',0,...
    'MarkerEdgeColor',fig_colors(1,:),'Color',fig_colors(1,:),'MarkerSize',fig_mark_size)
hold on
errorbar(hig_so4_dat.dD_CH4',hig_so4_dat.d13C_CH4',...
    hig_so4_dat.d13C_CH4_err,hig_so4_dat.d13C_CH4_err,...
    hig_so4_dat.dD_CH4_err,hig_so4_dat.dD_CH4_err,...
    'o','MarkerFaceColor','w','LineWidth',1,'CapSize',0,...
    'MarkerEdgeColor',fig_colors(2,:),'Color',fig_colors(2,:),'MarkerSize',fig_mark_size)
xlim([-200 300])
ylim([-70 15])
xlabel(['\deltaD-CH_4 (' char(8240) ')'])
ylabel(['\delta^{13}C-CH_4 (' char(8240) ')'])
set(gca,'FontSize',axes_font_size)



