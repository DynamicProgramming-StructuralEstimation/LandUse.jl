xvec=tvec;%1./tauvec;
%khaki grey: [0.5,0.5,0]
%green: [0.1,0.8,0.1]
%cinnamon: [0.67,0.31,0.02]

% plot(xvec,L10(:,j)./(L10(:,j)+L20(:,j)),'Color',[0.67,0.31,0.02],'LineWidth',2);
% xlim([xvec(1) xvec(end)])
% title('Share of Employment in Agriculture','FontSize',16);
% ax = gca % Get handle to current axes.
% ax.XColor = 'r'; % Red
% ax.YColor = 'b'; % Blue

figure
yyaxis left 
plot(xvec,L10(:,j)./(L10(:,j)+L20(:,j)),'Color',[0.1,0.8,0.1],'LineWidth',2);
%title('Agricultural shares','FontSize',16)
ax = gca % Get handle to current axes.
ax.XColor = 'k'; 
ax.YColor = [0.1,0.8,0.1]; 
%xlabel('Values from 0 to 25') 
xlim([xvec(1) xvec(end)])
ylabel('Agricultural Employment','Color',[0.1,0.8,0.1],'FontSize',14)
yyaxis right 
plot(xvec,S10(:,j),'Color',[0.67,0.31,0.02],'LineWidth',2);
ax = gca % Get handle to current axes.
ax.YColor = [0.67,0.31,0.02]; 
ylabel('Agricultural Land','Color',[0.67,0.31,0.02],'FontSize',14)  
saveas(gcf, 'plots/AgrEmpLand0','epsc');  

figure
plot(xvec,phi0(:,j),'LineWidth',2);
xlim([xvec(1) xvec(end)])
%title('City Size','FontSize',16);  
saveas(gcf, 'plots/CitySize0','epsc');  

figure
plot(xvec,L20(:,j)./phi0(:,j),'Color',[1,0,0],'LineWidth',2)
xlim([xvec(1) xvec(end)])
ylim([0 50])
%title('Urban Density','FontSize',16)
saveas(gcf, 'plots/UrbanDens0','epsc');  

figure
yyaxis left 
plot(xvec,D_cbd0(:,j),'Color',[1,0,0],'LineWidth',2);
%title('Urban density','FontSize',16)
xlim([xvec(1) xvec(end)])
ylim([0 200])
ax = gca % Get handle to current axes.
ax.XColor = 'k'; 
ax.YColor = [1,0,0]; 
ylabel('Density at the City Center','Color',[1,0,0],'FontSize',14)
yyaxis right 
plot(xvec,D_fringe0(:,j),'Color',[0 0 1],'LineWidth',2);
ylim([0 20])
ylabel('Density at the City Fringe','Color',[0 0 1],'FontSize',14)
ax = gca % Get handle to current axes.
ax.YColor = [0 0 1]; 
saveas(gcf, 'plots/DensCbdFringe0','epsc');  

figure
plot(xvec,(r0(:,j).*(L10(:,j)+L20(:,j))./P0(:,j))/(r0(10,j).*(L10(10,j)+L20(10,j))./P0(10,j)),'LineWidth',2);
xlim([xvec(1) xvec(end)])
%title('Real Land Rents Index','FontSize',16);
saveas(gcf, 'plots/LandRents0','epsc');  

% figure
% plot(xvec,(r0(:,j).*(L10(:,j)+L20(:,j))./P0(:,j))/(r0(end,j)./P0(end,j)),xvec,((r0(:,j).*(L10(:,j)+L20(:,j))-q10(:,j).*S10(:,j)-q10(:,j).*Shr0(:,j))./P0(:,j))/((r0(end,j)-q10(end,j).*S10(end,j)-q10(end,j).*Shr0(end,j))./P0(end,j)),'LineWidth',2);
% legend({'Aggregate','Urban'},'Location','NorthWest','FontSize',12)
% xlim([xvec(1) xvec(end)])
% title('Real Land Rents Index','FontSize',16);

figure
plot(xvec,100*(r0(:,j).*(L10(:,j)+L20(:,j)))./INC0(:,j),'LineWidth',2);
xlim([xvec(1) xvec(end)])
%title('Land Rents (as % of aggregate income)','FontSize',16);
saveas(gcf, 'plots/LandRentsOverIncome0','epsc');

% figure
% yyaxis left 
% plot(xvec,100*(r0(:,j).*(L10(:,j)+L20(:,j)))./INC0(:,j),'LineWidth',2);
% xlim([xvec(1) xvec(end)])
% title('Land Rents (as % of aggregate income)','FontSize',16);
% ylabel('Aggregate Land Rents','FontSize',14)
% yyaxis right
% plot(xvec,100*(r0(:,j).*(L10(:,j)+L20(:,j))-q10(:,j).*S10(:,j)-q10(:,j).*Shr0(:,j))./INCpc0(:,j),'LineWidth',2);
% ylabel('Urban Land Rents','FontSize',14) 

% figure
% plot(xvec,(r0(:,j).*(L10(:,j)+L20(:,j))./P0(:,j))/(r0(end,j)./P0(end,j)),xvec,((r0(:,j).*(L10(:,j)+L20(:,j))-q10(:,j).*S10(:,j)-q10(:,j).*Shr0(:,j))./P0(:,j))/((r0(end,j)-q10(end,j).*S10(end,j)-q10(end,j).*Shr0(end,j))./P0(end,j)),'LineWidth',2);
% legend({'Aggregate','Urban'},'Location','NorthWest','FontSize',12)
% xlim([xvec(1) xvec(end)])
% title('Real Land Rents Index','FontSize',16);

figure
yyaxis left 
plot(xvec,q_cbd0(:,j)./P0(:,j),'Color',[1,0,0],'LineWidth',2);
%title('Land Real Rental Price','FontSize',16)
xlim([xvec(1) xvec(end)])
ylabel('City Center Land Real Rental Price','Color',[1,0,0],'FontSize',14)
ax = gca % Get handle to current axes.
ax.XColor = 'k'; 
ax.YColor = [1,0,0]; 
yyaxis right 
plot(xvec,q_fringe0(:,j)./P0(:,j),'Color',[0,0,1],'LineWidth',2);
ylabel('City Fringe Land Real Rental Price','Color',[0,0,1],'FontSize',14)
ax = gca % Get handle to current axes.
ax.YColor = [0,0,1]; 
saveas(gcf, 'plots/LandPriceCBDFringe0','epsc');

% figure
% plot(xvec,q10(:,j)./P0(:,j),'LineWidth',2);
% xlim([xvec(1) xvec(end)])
% title('Farm Land Real Rental Price','FontSize',16);

figure
yyaxis left 
plot(xvec,qh_cbd0(:,j)./P0(:,j),'Color',[1,0,0],'LineWidth',2);
%title('Housing Rental Price','FontSize',16)
xlim([xvec(1) xvec(end)])
ylabel('City Center Housing Real Rental Price','Color',[1,0,0],'FontSize',14)
ax = gca % Get handle to current axes.
ax.XColor = 'k'; 
ax.YColor = [1,0,0]; 
yyaxis right 
plot(xvec,qh_fringe0(:,j)./P0(:,j),'Color',[0,0,1],'LineWidth',2);
ylabel('City Fringe Housing Real Rental Price','Color',[0,0,1],'FontSize',14)
ax = gca % Get handle to current axes.
ax.YColor = [0,0,1];
saveas(gcf, 'plots/HousePriceCBDFringe0','epsc'); 

figure
plot(xvec,w20(:,j)./w10(:,j),'LineWidth',2);
xlim([xvec(1) xvec(end)])
%title('Earnings Gap: Urban over Rural','FontSize',16);
saveas(gcf, 'plots/AGP0','epsc');

figure
plot(xvec,qh_cbd0(:,j)./qh_fringe0(:,j),'LineWidth',2);
xlim([xvec(1) xvec(end)])
%title('Housing Rent Gap: CBD vs fringe','FontSize',16);
saveas(gcf, 'plots/HousePriceGap0','epsc');

figure
plot(xvec,p10(:,j),'LineWidth',2);
xlim([xvec(1) xvec(end)])
%title('Agricultural Good Relative Price','FontSize',16);
saveas(gcf, 'plots/AgrPrice0','epsc');

figure
yyaxis left 
plot(xvec,H_cbd0(:,j),'LineWidth',2);
title('Urban Housing','FontSize',16)
xlim([xvec(1) xvec(end)])
ylabel('Housing at the City Center','FontSize',14)
yyaxis right 
plot(xvec,H_fringe0(:,j),'LineWidth',2);
ylabel('Housing at the City Fringe','FontSize',14) 



