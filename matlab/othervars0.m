gama2=gama/(1+epsilon); %gamma adjusted for housing supply elasticity

L20(i,j)=L-L10(i,j);
w20(i,j)=theta2*L20(i,j)^eta; %wage rate urban sector with no commuting costs
w10(i,j)=psi*w20(i,j)*(1-tau*phi0(i,j)); %wage rate rural sector
Shr0(i,j)=gama2/(q10(i,j)*chi1)*(w10(i,j)+r0(i,j)-p10(i,j)*cbar+sbar)*L10(i,j); %rural land used for residential purposes (housing plus wasted when chi1<1).
Y10(i,j)=theta1*(alpha*L10(i,j).^((sigma-1)/sigma)+(1-alpha)*S10(i,j).^((sigma-1)/sigma)).^(sigma/(sigma-1)); 
Y20(i,j)=theta2.*L20(i,j).^(1+eta);

norder=64;
[gnodes,gweights]=gausslegendre(norder); %Gauss-Legendre quadrature nodes and weights
dnodes=(gnodes+1)*(phi0(i,j)-0)/2+0;
gama2=gama/(1+epsilon); %gamma adjusted for housing supply elasticity
for n=1:norder
    d=dnodes(n);
    commvec(n)=(tau*d)*1/gama2*q10(i,j)/((w20(i,j)*(1-tau*phi0(i,j))+r0(i,j)-p10(i,j)*cbar+sbar)^(1/gama2))*(w20(i,j).*(1-tau.*d)+r0(i,j)-p10(i,j)*cbar+sbar).^(1/gama2-1);
end
intcommcosts0(i,j)=(phi0(i,j)-0)/2*ones(1,norder)*(gweights.*commvec)';
GDP0(i,j)=(p10(i,j)*Y10(i,j)+Y20(i,j))/L; %value of GDP per capita
INC0(i,j)=GDP0(i,j)-intcommcosts0(i,j); %value of disposable income (GDP net of commuting costs). TO CHECK: should we multiply intcommcosts0 by L?       
EAR0(i,j)=GDP0(i,j)-intcommcosts0(i,j)-r0(i,j)*L; %value of average labor earnings (GDP net of commuting costs minus land rents). TO CHECK: should we multiply intcommcosts0 by L?

q_fringe0(i,j)=q10(i,j); %price of a unit of land at the fringe
qh_fringe0(i,j)=((1+epsilon)*q_fringe0(i,j))^(1/(1+epsilon));%price of a unit of housing at the fringe
H_fringe0(i,j)=qh_fringe0(i,j)^epsilon; %total housing at the fringe
D_fringe0(i,j)=1/gama*(qh_fringe0(i,j)^(1+epsilon))/(w20(i,j)*(1-tau*phi0(i,j))+r0(i,j)-p10(i,j)*cbar+sbar); %population density at the fringe

qh_cbd0(i,j)=qh_fringe0(i,j)*((w20(i,j)+r0(i,j)-p10(i,j)*cbar+sbar)/(w20(i,j).*(1-tau.*phi0(i,j))+r0(i,j)-p10(i,j)*cbar+sbar))^(1/gama);%price of a unit of housing at the cbd
q_cbd0(i,j)=1/(1+epsilon)*qh_cbd0(i,j)^(1+epsilon);%price of a unit of land at the cbd
H_cbd0(i,j)=qh_cbd0(i,j)^epsilon; %total housing at the cbd
D_cbd0(i,j)=1/gama*(qh_cbd0(i,j)^(1+epsilon))/(w20(i,j)+r0(i,j)-p10(i,j)*cbar+sbar); %population density at the cbd

% phim=phi0(i,j)/2;
% D_m0(i,j)=Hsupply_f(phim,w20(i,j),r0(i,j),p10(i,j),q10(i,j),phi0(i,j))/hdemand_f(phim,w20(i,j),r0(i,j),p10(i,j),q10(i,j),phi0(i,j));

if i==1;
    P0(i,j)=p10(1,j);
else
    Pgr_laspeyres0(i,j)=(p10(i,j).*Y10(i-1,j)+Y20(i-1,j))./(p10(i-1,j).*Y10(i-1,j)+Y20(i-1,j));
    Pgr_paasche0(i,j)=(p10(i,j).*Y1_f(L10(i,j),S10(i,j))+Y2_f(L20(i,j)))./(p10(i-1,j).*Y1_f(L10(i,j),S10(i,j))+Y2_f(L20(i,j)));
    Pgr0(i,j)=(Pgr_laspeyres0(i,j).^0.5).*(Pgr_paasche0(i,j).^0.5);
    P0(i,j)=P0(i-1,j)*Pgr0(i,j);
end
