%file to compute plots comparing constant epsilon results with epsilon-specific results
xvec=tvec;

% plot(xvec,L10(:,j)./(L10(:,j)+L20(:,j)),'--',xvec,L1(:,j)./(L1(:,j)+L2(:,j)),'color','r','LineWidth',2)
% xlim([xvec(1) xvec(end)])
% legend({'Constant \epsilon=4','Location-specific \epsilon'},'Location','NorthEast','FontSize',12)
% title('Share of Employment in Agriculture','FontSize',16);
figure
yyaxis left 
plot(xvec,L10(:,j)./(L10(:,j)+L20(:,j)),'--',xvec,L1(:,j)./(L1(:,1)+L2(:,j)),'-','Color',[0.1,0.8,0.1],'LineWidth',2);
%title('Agricultural shares','FontSize',16)
ax = gca % Get handle to current axes.
ax.XColor = 'k'; 
ax.YColor = [0.1,0.8,0.1]; 
%xlabel('Values from 0 to 25') 
xlim([xvec(1) xvec(end)])
ylabel('Agricultural Employment','Color',[0.1,0.8,0.1],'FontSize',14)
legend({'Constant \epsilon=4','Location-specific \epsilon'},'Location','NorthEast','FontSize',12)
yyaxis right 
plot(xvec,S10(:,j),'--',xvec,S1(:,j),'-','Color',[0.67,0.31,0.02],'LineWidth',2);
ax = gca % Get handle to current axes.
ax.YColor = [0.67,0.31,0.02]; 
ylabel('Agricultural Land','Color',[0.67,0.31,0.02],'FontSize',14)  
saveas(gcf, 'plots/AgrEmpLand_comp','epsc');  

figure
plot(xvec,phi0(:,j),'--',xvec,phi(:,j),'color','b','LineWidth',2);
xlim([xvec(1) xvec(end)])
legend({'Constant \epsilon=4','Location-specific \epsilon'},'Location','NorthWest','FontSize',12)
%title('City Size','FontSize',16);  
saveas(gcf, 'plots/CitySize_comp','epsc');  

figure
plot(xvec,L20(:,j)./phi0(:,j),'--',xvec,L2(:,j)./phi(:,j),'Color',[1,0,0],'LineWidth',2)
xlim([xvec(1) xvec(end)])
ylim([0 50])
legend({'Constant \epsilon=4','Location-specific \epsilon'},'Location','NorthEast','FontSize',12)
%title('Urban Density','FontSize',16)
saveas(gcf, 'plots/UrbanDens_comp','epsc');  

figure
semilogy(xvec,L20(:,j)./phi0(:,j),'--',xvec,L2(:,j)./phi(:,j),'Color',[1,0,0],'LineWidth',2)
xlim([xvec(1) xvec(end)])
ylim([0 50])
legend({'Constant \epsilon=4','Location-specific \epsilon'},'Location','NorthEast','FontSize',12)
%title('Urban Density','FontSize',16)
saveas(gcf, 'plots/LogUrbanDens_comp','epsc'); 

figure
yyaxis left 
plot(xvec,D_cbd0(:,j),'--',xvec,D_cbd(:,j),'-','Color',[1,0,0],'LineWidth',2);
%title('Urban density','FontSize',16)
xlim([xvec(1) xvec(end)])
ylim([0 140])
ax = gca % Get handle to current axes.
ax.XColor = 'k'; 
ax.YColor = [1,0,0]; 
ylabel('Density at the City Center','Color',[1,0,0],'FontSize',14)
legend({'Constant \epsilon=4','Location-specific \epsilon'},'Location','SouthWest','FontSize',12)
yyaxis right 
plot(xvec,D_fringe0(:,j),'--',xvec,D_fringe(:,j),'-','Color',[0 0 1],'LineWidth',2);
ylim([0 11.5])
ylabel('Density at the City Fringe','Color',[0 0 1],'FontSize',14)
ax = gca % Get handle to current axes.
ax.YColor = [0 0 1]; 
saveas(gcf, 'plots/DensCbdFringe_comp','epsc');  

figure
yyaxis left 
semilogy(xvec,D_cbd0(:,j),'--',xvec,D_cbd(:,j),'-','Color',[1,0,0],'LineWidth',2);
%title('Urban density','FontSize',16)
xlim([xvec(1) xvec(end)])
ylim([0 140])
ax = gca % Get handle to current axes.
ax.XColor = 'k'; 
ax.YColor = [1,0,0]; 
ylabel('Density at the City Center','Color',[1,0,0],'FontSize',14)
legend({'Constant \epsilon=4','Location-specific \epsilon'},'Location','SouthWest','FontSize',12)
yyaxis right 
semilogy(xvec,D_fringe0(:,j),'--',xvec,D_fringe(:,j),'-','Color',[0 0 1],'LineWidth',2);
ylim([0 20])
ylabel('Density at the City Fringe','Color',[0 0 1],'FontSize',14)
ax = gca % Get handle to current axes.
ax.YColor = [0 0 1]; 
saveas(gcf, 'plots/LogDensCbdFringe_comp','epsc');  

figure
plot(xvec,(r0(:,j).*(L10(:,j)+L20(:,j))./P0(:,j))/(r0(10,j).*(L10(10,j)+L20(10,j))./P0(10,j)),'--',xvec,(r(:,j).*(L1(:,j)+L2(:,j))./P(:,j))/(r(10,j).*(L1(10,j)+L2(10,j))./P(10,j)),'-','color','b','LineWidth',2);
xlim([xvec(1) xvec(end)])
legend({'Constant \epsilon=4','Location-specific \epsilon'},'Location','NorthWest','FontSize',12)
%title('Real Land Rents Index','FontSize',16);
saveas(gcf, 'plots/LandRents_comp','epsc');  

figure
%title('Land Rents (as % of aggregate income)','FontSize',16);
yyaxis left 
plot(xvec,100*(r0(:,j).*(L10(:,j)+L20(:,j)))./INC0(:,j),'--',xvec,100*(r(:,j).*(L1(:,j)+L2(:,j)))./INC(:,j),'-','color','r','LineWidth',2);
xlim([xvec(1) xvec(end)])
legend({'Constant \epsilon=4','Location-specific \epsilon'},'Location','NorthEast','FontSize',12)
ylabel('Land Rents (as % of aggregate income)','Color',[1,0,0],'FontSize',12)
ax = gca % Get handle to current axes.
ax.XColor = 'k'; 
ax.YColor = [1,0,0]; 
yyaxis right 
plot(xvec,100*(r0(:,j).*(L10(:,j)+L20(:,j))-q10(:,j).*S10(:,j)-q10(:,j).*Shr0(:,j))./INC0(:,j),'--',xvec,100*(r(:,j).*(L1(:,j)+L2(:,j))-q1(:,j).*S1(:,j)-q1(:,j).*Shr(:,j))./INC(:,j),'-','color','b','LineWidth',2);
xlim([xvec(1) xvec(end)])
legend({'Constant \epsilon=4','Location-specific \epsilon'},'Location','North','FontSize',12)
ylabel('Urban Land Rents (as % of aggregate income)','Color','b','FontSize',12)
ax = gca % Get handle to current axes.
ax.YColor = [0,0,1]; 
saveas(gcf, 'plots/LandRentsOverIncome_comp','epsc');

figure
yyaxis left 
plot(xvec,q_cbd0(:,j)./P0(:,j),'--',xvec,q_cbd(:,j)./P(:,j),'-','Color',[1,0,0],'LineWidth',2);
%title('Land Real Rental Price','FontSize',16)
xlim([xvec(1) xvec(end)])
ylabel('City Center Land Real Rental Price','Color',[1,0,0],'FontSize',14)
ax = gca % Get handle to current axes.
ax.XColor = 'k'; 
ax.YColor = [1,0,0]; 
yyaxis right 
legend({'Constant \epsilon=4','Location-specific \epsilon'},'Location','NorthWest','FontSize',12)
plot(xvec,q_fringe0(:,j)./P0(:,j),'--',xvec,q_fringe(:,j)./P(:,j),'-','Color',[0,0,1],'LineWidth',2);
ylabel('City Fringe Land Real Rental Price','Color',[0,0,1],'FontSize',14)
ax = gca % Get handle to current axes.
ax.YColor = [0,0,1]; 
saveas(gcf, 'plots/LandPriceCBDFringe_comp','epsc');

figure
yyaxis left 
plot(xvec,qh_cbd0(:,j)./P0(:,j),'--',xvec,qh_cbd(:,j)./P(:,j),'-','Color',[1,0,0],'LineWidth',2);
%title('Housing Rental Price','FontSize',16)
xlim([xvec(1) xvec(end)])
ylabel('City Center Housing Real Rental Price','Color',[1,0,0],'FontSize',14)
ax = gca % Get handle to current axes.
ax.XColor = 'k'; 
ax.YColor = [1,0,0]; 
legend({'Constant \epsilon=4','Location-specific \epsilon'},'Location','NorthWest','FontSize',12)
yyaxis right 
plot(xvec,qh_fringe0(:,j)./P0(:,j),'--',xvec,qh_fringe(:,j)./P(:,j),'-','Color',[0,0,1],'LineWidth',2);
ylabel('City Fringe Housing Real Rental Price','Color',[0,0,1],'FontSize',14)
ax = gca % Get handle to current axes.
ax.YColor = [0,0,1];
saveas(gcf, 'plots/HousePriceCBDFringe_comp','epsc'); 

figure
plot(xvec,w20(:,j)./w10(:,j),'--',xvec,w2(:,j)./w1(:,j),'-','color','b','LineWidth',2);
xlim([xvec(1) xvec(end)])
legend({'Constant \epsilon=4','Location-specific \epsilon'},'Location','NorthWest','FontSize',12)
%title('Earnings Gap: Urban over Rural','FontSize',16);
saveas(gcf, 'plots/AGP_comp','epsc');

figure
plot(xvec,qh_cbd0(:,j)./qh_fringe0(:,j),'--',xvec,qh_cbd(:,j)./qh_fringe(:,j),'-','color','b','LineWidth',2);
xlim([xvec(1) xvec(end)])
legend({'Constant \epsilon=4','Location-specific \epsilon'},'Location','NorthWest','FontSize',12)
%title('Housing Rent Gap: CBD vs fringe','FontSize',16);
saveas(gcf, 'plots/HousePriceGap_comp','epsc');

figure
plot(xvec,p10(:,j),'--',xvec,p1(:,j),'-','color','b','LineWidth',2);
xlim([xvec(1) xvec(end)])
legend({'Constant \epsilon=4','Location-specific \epsilon'},'Location','NorthEast','FontSize',12)
%title('Agricultural Good Relative Price','FontSize',16);
saveas(gcf, 'plots/AgrPrice_comp','epsc');

figure
yyaxis left 
plot(xvec,H_cbd0(:,j),'--',xvec,H_cbd(:,j),'-','LineWidth',2);
title('Total Urban Housing','FontSize',16)
xlim([xvec(1) xvec(end)])
ylabel('Housing at the City Center','FontSize',14)
yyaxis right 
plot(xvec,H_fringe0(:,j),'--',xvec,H_fringe(:,j),'-','LineWidth',2);
ylabel('Housing at the City Fringe','FontSize',14) 
saveas(gcf, 'plots/TotalHousingCBDFringe_comp','epsc');


