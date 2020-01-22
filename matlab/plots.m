xvec=tvec;%1./tauvec;
%khaki grey: [0.5,0.5,0]
%green: [0.1,0.8,0.1]
%cinnamon: [0.67,0.31,0.02]

% plot(xvec,L1(:,j)./(L1(:,j)+L2(:,j)),'Color',[0.67,0.31,0.02],'LineWidth',2);
% xlim([xvec(1) xvec(end)])
% title('Share of Employment in Agriculture','FontSize',16);
% ax = gca % Get handle to current axes.
% ax.XColor = 'r'; % Red
% ax.YColor = 'b'; % Blue

figure
yyaxis left 
plot(xvec,L1(:,j)./(L1(:,j)+L2(:,j)),'Color',[0.1,0.8,0.1],'LineWidth',2);
%title('Agricultural shares','FontSize',16)
ax = gca % Get handle to current axes.
ax.XColor = 'k'; 
ax.YColor = [0.1,0.8,0.1]; 
%xlabel('Values from 0 to 25') 
xlim([xvec(1) xvec(end)])
ylabel('Agricultural Employment','Color',[0.1,0.8,0.1],'FontSize',14)
yyaxis right 
plot(xvec,S1(:,j),'Color',[0.67,0.31,0.02],'LineWidth',2);
ax = gca % Get handle to current axes.
ax.YColor = [0.67,0.31,0.02]; 
ylabel('Agricultural Land','Color',[0.67,0.31,0.02],'FontSize',14)  
saveas(gcf, 'plots/AgrEmpLand','epsc');  

figure
plot(xvec,phi(:,j),'LineWidth',2);
xlim([xvec(1) xvec(end)])
%title('City Size','FontSize',16);  
saveas(gcf, 'plots/CitySize','epsc');  

figure
plot(xvec,L2(:,j)./phi(:,j),'Color',[1,0,0],'LineWidth',2)
xlim([xvec(1) xvec(end)])
ylim([0 50])
%title('Urban Density','FontSize',16)
saveas(gcf, 'plots/UrbanDens','epsc');  

figure
yyaxis left 
plot(xvec,D_cbd(:,j),'Color',[1,0,0],'LineWidth',2);
%title('Urban density','FontSize',16)
xlim([xvec(1) xvec(end)])
ylim([0 200])
ax = gca % Get handle to current axes.
ax.XColor = 'k'; 
ax.YColor = [1,0,0]; 
ylabel('Density at the City Center','Color',[1,0,0],'FontSize',14)
yyaxis right 
plot(xvec,D_fringe(:,j),'Color',[0 0 1],'LineWidth',2);
ylim([0 20])
ylabel('Density at the City Fringe','Color',[0 0 1],'FontSize',14)
ax = gca % Get handle to current axes.
ax.YColor = [0 0 1]; 
saveas(gcf, 'plots/DensCbdFringe','epsc');  

figure
plot(xvec,(r(:,j).*(L1(:,j)+L2(:,j))./P(:,j))/(r(10,j).*(L1(10,j)+L2(10,j))./P(10,j)),'LineWidth',2);
xlim([xvec(1) xvec(end)])
%title('Real Land Rents Index','FontSize',16);
saveas(gcf, 'plots/LandRents','epsc');  

figure
plot(xvec,100*(r(:,j).*(L1(:,j)+L2(:,j)))./INC(:,j),'LineWidth',2);
xlim([xvec(1) xvec(end)])
%title('Land Rents (as % of aggregate income)','FontSize',16);
saveas(gcf, 'plots/LandRentsOverIncome','epsc');

% figure
% yyaxis left 
% plot(xvec,100*(r(:,j).*(L1(:,j)+L2(:,j)))./INC(:,j),'LineWidth',2);
% xlim([xvec(1) xvec(end)])
% title('Land Rents (as % of aggregate income)','FontSize',16);
% ylabel('Aggregate Land Rents','FontSize',14)
% yyaxis right
% plot(xvec,100*(r(:,j).*(L1(:,j)+L2(:,j))-q1(:,j).*S1(:,j)-q1(:,j).*Shr(:,j))./INC(:,j),'LineWidth',2);
% ylabel('Urban Land Rents','FontSize',14) 

figure
yyaxis left 
plot(xvec,q_cbd(:,j)./P(:,j),'Color',[1,0,0],'LineWidth',2);
%title('Land Real Rental Price','FontSize',16)
xlim([xvec(1) xvec(end)])
ylabel('City Center Land Real Rental Price','Color',[1,0,0],'FontSize',14)
ax = gca % Get handle to current axes.
ax.XColor = 'k'; 
ax.YColor = [1,0,0]; 
yyaxis right 
plot(xvec,q_fringe(:,j)./P(:,j),'Color',[0,0,1],'LineWidth',2);
ylabel('City Fringe Land Real Rental Price','Color',[0,0,1],'FontSize',14)
ax = gca % Get handle to current axes.
ax.YColor = [0,0,1]; 
saveas(gcf, 'plots/LandPriceCBDFringe','epsc');

% figure
% plot(xvec,q1(:,j)./P(:,j),'LineWidth',2);
% xlim([xvec(1) xvec(end)])
% title('Farm Land Real Rental Price','FontSize',16);

figure
yyaxis left 
plot(xvec,qh_cbd(:,j)./P(:,j),'Color',[1,0,0],'LineWidth',2);
%title('Housing Rental Price','FontSize',16)
xlim([xvec(1) xvec(end)])
ylabel('City Center Housing Real Rental Price','Color',[1,0,0],'FontSize',14)
ax = gca % Get handle to current axes.
ax.XColor = 'k'; 
ax.YColor = [1,0,0]; 
yyaxis right 
plot(xvec,qh_fringe(:,j)./P(:,j),'Color',[0,0,1],'LineWidth',2);
ylabel('City Fringe Housing Real Rental Price','Color',[0,0,1],'FontSize',14)
ax = gca % Get handle to current axes.
ax.YColor = [0,0,1];
saveas(gcf, 'plots/HousePriceCBDFringe','epsc'); 

figure
plot(xvec,w2(:,j)./w1(:,j),'LineWidth',2);
xlim([xvec(1) xvec(end)])
%title('Earnings Gap: Urban over Rural','FontSize',16);
saveas(gcf, 'plots/AGP','epsc');

figure
plot(xvec,qh_cbd(:,j)./qh_fringe(:,j),'LineWidth',2);
xlim([xvec(1) xvec(end)])
%title('Housing Rent Gap: CBD vs fringe','FontSize',16);
saveas(gcf, 'plots/HousePriceGap','epsc');

figure
plot(xvec,p1(:,j),'LineWidth',2);
xlim([xvec(1) xvec(end)])
%title('Agricultural Good Relative Price','FontSize',16);
saveas(gcf, 'plots/AgrPrice','epsc');

figure
yyaxis left 
plot(xvec,H_cbd(:,j),'LineWidth',2);
title('Total Urban Housing','FontSize',16)
xlim([xvec(1) xvec(end)])
ylabel('Housing at the City Center','FontSize',14)
yyaxis right 
plot(xvec,H_fringe(:,j),'LineWidth',2);
ylabel('Housing at the City Fringe','FontSize',14) 
saveas(gcf, 'plots/TotalHousingCBDFringe','epsc');



