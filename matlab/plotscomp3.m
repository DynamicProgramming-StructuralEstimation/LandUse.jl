xvec=tvec;

% plot(xvec,L1(:,1)./(L1(:,1)+L2(:,1)),'--',xvec,L1(:,3)./(L1(:,3)+L2(:,3)),'color','r','LineWidth',2)
% xlim([xvec(1) xvec(end)])
% legend({'Constant \tau, Constant L','Decreasing \tau, Increasing L'},'Location','NorthEast','FontSize',12)
% title('Share of Employment in Agriculture','FontSize',16);

figure
yyaxis left 
plot(xvec,L1(:,1)./(L1(:,1)+L2(:,1)),'--',xvec,L1(:,3)./(L1(:,3)+L2(:,3)),'-','Color',[0.1,0.8,0.1],'LineWidth',2);
%title('Agricultural shares','FontSize',16)
ax = gca % Get handle to current axes.
ax.XColor = 'k'; 
ax.YColor = [0.1,0.8,0.1]; 
%xlabel('Values from 0 to 25') 
xlim([xvec(1) xvec(end)])
ylabel('Agricultural Employment','Color',[0.1,0.8,0.1],'FontSize',14)
legend({'Constant \tau, Constant L','Decreasing \tau, Increasing L'},'Location','NorthEast','FontSize',12)
yyaxis right 
plot(xvec,S1(:,1),'--',xvec,S1(:,3),'-','Color',[0.67,0.31,0.02],'LineWidth',2);
ax = gca % Get handle to current axes.
ax.YColor = [0.67,0.31,0.02]; 
ylabel('Agricultural Land','Color',[0.67,0.31,0.02],'FontSize',14)  
saveas(gcf, 'plots/AgrEmpLand_comp3','epsc');  

figure
plot(xvec,phi(:,1),'--',xvec,phi(:,3),'color','b','LineWidth',2);
xlim([xvec(1) xvec(end)])
legend({'Constant \tau, Constant L','Decreasing \tau, Increasing L'},'Location','NorthWest','FontSize',12)
%title('City Size','FontSize',16);  
saveas(gcf, 'plots/CitySize_comp3','epsc');  

figure
plot(xvec,L2(:,1)./phi(:,1),'--',xvec,L2(:,3)./phi(:,3),'Color',[1,0,0],'LineWidth',2)
xlim([xvec(1) xvec(end)])
ylim([0 75])
legend({'Constant \tau, Constant L','Decreasing \tau, Increasing L'},'Location','SouthWest','FontSize',12)
%title('Urban Density','FontSize',16)
saveas(gcf, 'plots/UrbanDens_comp3','epsc');  

figure
semilogy(xvec,L2(:,1)./phi(:,1),'--',xvec,L2(:,3)./phi(:,3),'Color',[1,0,0],'LineWidth',2)
xlim([xvec(1) xvec(end)])
ylim([0 75])
legend({'Constant \tau, Constant L','Decreasing \tau, Increasing L'},'Location','SouthWest','FontSize',12)
%title('Urban Density','FontSize',16)
saveas(gcf, 'plots/LogUrbanDens_comp3','epsc'); 

figure
yyaxis left 
plot(xvec,D_cbd(:,1),'--',xvec,D_cbd(:,3),'-','Color',[1,0,0],'LineWidth',2);
%title('Urban density','FontSize',16)
xlim([xvec(1) xvec(end)])
%ylim([0 400])
ax = gca % Get handle to current axes.
ax.XColor = 'k'; 
ax.YColor = [1,0,0]; 
ylabel('Density at the City Center','Color',[1,0,0],'FontSize',14)
legend({'Constant \tau, Constant L','Decreasing \tau, Increasing L'},'Location','NorthEast','FontSize',12)
yyaxis right 
plot(xvec,D_fringe(:,1),'--',xvec,D_fringe(:,3),'-','Color',[0 0 1],'LineWidth',2);
ylim([0 20])
ylabel('Density at the City Fringe','Color',[0 0 1],'FontSize',14)
ax = gca % Get handle to current axes.
ax.YColor = [0 0 1]; 
saveas(gcf, 'plots/DensCbdFringe_comp3','epsc'); 

figure
yyaxis left 
semilogy(xvec,D_cbd(:,1),'--',xvec,D_cbd(:,3),'-','Color',[1,0,0],'LineWidth',2);
%title('Urban density','FontSize',16)
xlim([xvec(1) xvec(end)])
ylim([0 220])
ax = gca % Get handle to current axes.
ax.XColor = 'k'; 
ax.YColor = [1,0,0]; 
ylabel('Density at the City Center','Color',[1,0,0],'FontSize',14)
%legend({'Constant \epsilon=4','Location-specific \epsilon'},'Location','SouthWest','FontSize',12)
yyaxis right 
semilogy(xvec,D_fringe(:,1),'--',xvec,D_fringe(:,3),'-','Color',[0 0 1],'LineWidth',2);
ylim([0 60])
ylabel('Density at the City Fringe','Color',[0 0 1],'FontSize',14)
ax = gca % Get handle to current axes.
ax.YColor = [0 0 1]; 
saveas(gcf, 'plots/LogDensCbdFringe_comp3','epsc'); 

figure
plot(xvec,(r(:,1).*(L1(:,1)+L2(:,1))./P(:,1))/(r(10,1).*(L1(10,1)+L2(10,1))./P(10,1)),'--',xvec,(r(:,3).*(L1(:,3)+L2(:,3))./P(:,3))/(r(10,3).*(L1(10,3)+L2(10,3))./P(10,3)),'-','color','b','LineWidth',2);
xlim([xvec(1) xvec(end)])
legend({'Constant \tau, Constant L','Decreasing \tau, Increasing L'},'Location','NorthWest','FontSize',12)
%title('Real Land Rents Index','FontSize',16);
saveas(gcf, 'plots/LandRents_comp3','epsc');  

figure
%title('Land Rents (as % of aggregate income)','FontSize',16);
yyaxis left 
plot(xvec,100*(r(:,1).*(L1(:,1)+L2(:,1)))./INC(:,1),'--',xvec,100*(r(:,3).*(L1(:,3)+L2(:,3)))./INC(:,3),'-','color','r','LineWidth',2);
xlim([xvec(1) xvec(end)])
ylim([0 35])
legend({'Constant \tau, Constant L','Decreasing \tau, Increasing L'},'Location','NorthEast','FontSize',12)
ylabel('Land Rents (as % of aggregate income)','Color',[1,0,0],'FontSize',12)
ax = gca % Get handle to current axes.
ax.XColor = 'k'; 
ax.YColor = [1,0,0]; 
yyaxis right 
plot(xvec,100*(r(:,1).*(L1(:,1)+L2(:,1))-q1(:,1).*S1(:,1)-q1(:,1).*Shr(:,1))./INC(:,1),'--',xvec,100*(r(:,3).*(L1(:,3)+L2(:,3))-q1(:,3).*S1(:,3)-q1(:,3).*Shr(:,3))./INC(:,3),'-','color','b','LineWidth',2);
xlim([xvec(1) xvec(end)])
legend({'Constant \tau, Constant L','Decreasing \tau, Increasing L'},'Location','North','FontSize',12)
ylabel('Urban Land Rents (as % of aggregate income)','Color','b','FontSize',12)
ax = gca % Get handle to current axes.
ax.YColor = [0,0,1]; 
saveas(gcf, 'plots/LandRentsOverIncome_comp3','epsc');

figure
yyaxis left 
plot(xvec,q_cbd(:,1)./P(:,1),'--',xvec,q_cbd(:,3)./P(:,3),'-','Color',[1,0,0],'LineWidth',2);
%title('Land Real Rental Price','FontSize',16)
xlim([xvec(1) xvec(end)])
ylabel('City Center Land Real Rental Price','Color',[1,0,0],'FontSize',14)
ax = gca % Get handle to current axes.
ax.XColor = 'k'; 
ax.YColor = [1,0,0]; 
yyaxis right 
legend({'Constant \tau, Constant L','Decreasing \tau, Increasing L'},'Location','NorthWest','FontSize',12)
plot(xvec,q_fringe(:,1)./P(:,1),'--',xvec,q_fringe(:,3)./P(:,3),'-','Color',[0,0,1],'LineWidth',2);
ylabel('City Fringe Land Real Rental Price','Color',[0,0,1],'FontSize',14)
ax = gca % Get handle to current axes.
ax.YColor = [0,0,1]; 
saveas(gcf, 'plots/LandPriceCBDFringe_comp3','epsc');

figure
yyaxis left 
plot(xvec,qh_cbd(:,1)./P(:,1),'--',xvec,qh_cbd(:,3)./P(:,3),'-','Color',[1,0,0],'LineWidth',2);
%title('Housing Rental Price','FontSize',16)
xlim([xvec(1) xvec(end)])
ylabel('City Center Housing Real Rental Price','Color',[1,0,0],'FontSize',14)
ax = gca % Get handle to current axes.
ax.XColor = 'k'; 
ax.YColor = [1,0,0]; 
legend({'Constant \tau, Constant L','Decreasing \tau, Increasing L'},'Location','NorthWest','FontSize',12)
yyaxis right 
plot(xvec,qh_fringe(:,1)./P(:,1),'--',xvec,qh_fringe(:,3)./P(:,3),'-','Color',[0,0,1],'LineWidth',2);
ylabel('City Fringe Housing Real Rental Price','Color',[0,0,1],'FontSize',14)
ax = gca % Get handle to current axes.
ax.YColor = [0,0,1];
saveas(gcf, 'plots/HousePriceCBDFringe_comp3','epsc'); 

figure
plot(xvec,w2(:,1)./w1(:,1),'--',xvec,w2(:,3)./w1(:,3),'-','color','b','LineWidth',2);
xlim([xvec(1) xvec(end)])
legend({'Constant \tau, Constant L','Decreasing \tau, Increasing L'},'Location','NorthWest','FontSize',12)
%title('Earnings Gap: Urban over Rural','FontSize',16);
saveas(gcf, 'plots/AGP_comp3','epsc');

figure
plot(xvec,qh_cbd(:,1)./qh_fringe(:,1),'--',xvec,qh_cbd(:,3)./qh_fringe(:,3),'-','color','b','LineWidth',2);
xlim([xvec(1) xvec(end)])
legend({'Constant \tau, Constant L','Decreasing \tau, Increasing L'},'Location','NorthWest','FontSize',12)
%title('Housing Rent Gap: CBD vs fringe','FontSize',16);
saveas(gcf, 'plots/HousePriceGap_comp3','epsc');

figure
plot(xvec,p1(:,1),'--',xvec,p1(:,3),'-','color','b','LineWidth',2);
xlim([xvec(1) xvec(end)])
ylim([0 1.3])
legend({'Constant \tau, Constant L','Decreasing \tau, Increasing L'},'Location','NorthEast','FontSize',12)
%title('Agricultural Good Relative Price','FontSize',16);
saveas(gcf, 'plots/AgrPrice_comp3','epsc');

figure
yyaxis left 
plot(xvec,H_cbd(:,1),'--',xvec,H_cbd(:,3),'-','LineWidth',3);
title('Total Urban Housing','FontSize',16)
xlim([xvec(1) xvec(end)])
ylabel('Housing at the City Center','FontSize',14)
legend({'Constant \tau, Constant L','Decreasing \tau, Increasing L'},'Location','NorthWest','FontSize',12)
yyaxis right 
plot(xvec,H_fringe(:,1),'--',xvec,H_fringe(:,3),'-','LineWidth',2);
ylim([0 1])
ylabel('Housing at the City Fringe','FontSize',14) 
saveas(gcf, 'plots/TotalHousingCBDFringe_comp3','epsc');


