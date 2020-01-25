function F=EqsysCD_pfixed_f(x) %"small open economy model"
% Equilibrium equations when relative price of agr good is exogenous 

global gama nu cbar sbar theta2 theta1 alpha lambda tau chi1 chi2 L epsilon eta psi sigma

gama2=gama/(1+epsilon); %gamma adjusted for housing supply elasticity

q1=x(1); %land price in rural sector
phi=x(2); %city radius
r=x(3); %land rents
L1=x(4); %employment in rural sector
p1=x(5); %relative price of rural-sector good
S1=x(6); %land used as farm sector input

L2=L-L1; %employment in urban sector
w2=theta2*L2^eta; %wage rate urban sector with no commuting costs
gama2=gama/(1+epsilon); %gamma adjusted for housing supply elasticity
w1=psi*w2*(1-tau*phi); %wage rate rural sector
S1h=gama2/(q1*chi1)*(w1+r-p1*cbar+sbar)*L1; %rural land used for residential purposes (housing plus wasted when chi1<1).

%F(1)=p1*alpha*theta1*(S1/L1)^(1-alpha)-w1; %rural firm labor optimization
F(1)=p1*alpha*theta1*L1^(-1/sigma)*(alpha*L1^((sigma-1)/sigma)+(1-alpha)*S1^((sigma-1)/sigma))^(1/(sigma-1))-w1; %rural firm labor optimization

%F(2)=p1*(1-alpha)*theta1*(L1/S1)^alpha-q1; %rural firm land optimization
F(2)=p1*(1-alpha)*theta1*S1^(-1/sigma)*(alpha*L1^((sigma-1)/sigma)+(1-alpha)*S1^((sigma-1)/sigma))^(1/(sigma-1))-q1; %rural firm land optimization

F(3)=1+tau*w2*L2/(chi2*q1)-((w2+r-p1*cbar+sbar)/(w2*(1-tau*phi)+r-p1*cbar+sbar))^(1/gama2); %urban workers accomodation

%F(4)=q1*(1-lambda-phi)+q1/(tau*w2)*(((w2+r-p1*cbar+sbar)^(1/gama+1))/((w2*(1-tau*phi)+r-p1*cbar+sbar)^(1/gama))-(w2*(1-tau*phi)+r-p1*cbar+sbar))-r*L; %land rents
F(4)=q1*S1+q1*chi1*S1h...
    +chi2*gama2*q1/(tau*(1+gama2)*w2)*((w2+r-p1*cbar+sbar)^(1/gama2+1)/((w2*(1-tau*phi)+r-p1*cbar+sbar)^(1/gama2))-(w2*(1-tau*phi)+r-p1*cbar+sbar))...
    -r*L; %land rents                                                                                                                                                                                                                                                     
%F(4)=r-0;

F(5)=1-lambda-phi-S1h-S1; %land market clearing   
F(6)=p1-2.3866;
% F(6)=L1*(1-gama)*(1-nu)*(w1+r-p1*cbar+sbar)...
%      +chi2*(1-gama)*(1-nu)*q1/((1+gama2)*w2*tau)*((w2+r-p1*cbar+sbar)^(1/gama2+1)/((w2*(1-tau*phi)+r-p1*cbar+sbar)^(1/gama2))-(w2*(1-tau*phi)+r-p1*cbar+sbar))...
%      +chi2*gama2*q1/((1+gama2)*tau*w2)*(((w2+r-p1*cbar+sbar)^(1/gama2+1))/((w2*(1-tau*phi)+r-p1*cbar+sbar)^(1/gama2))-(w2*(1-tau*phi)+r-p1*cbar+sbar))...
%      -phi*chi2*q1...
%      +epsilon*q1*S1h...
%      +epsilon*chi2*gama2*q1/(tau*(1+gama2)*w2)*((w2+r-p1*cbar+sbar)^(1/gama2+1)/((w2*(1-tau*phi)+r-p1*cbar+sbar)^(1/gama2))-(w2*(1-tau*phi)+r-p1*cbar+sbar))...
%      -sbar*L...
%      -theta2*L2^(1+eta); %urban good market clearing condition: here in the second term missing w2 at denominator
F

%F(3)=q1+tau/chi2*w2*L2-q1*((w2+r+p1*cbar)/(w2*(1-tau*phi)+r-p1*cbar))^(1/gama); %urban workers accomodation
%F(4)=tau*w2*q1*S1+tau*w2*q1*chi1*S1h+...
%    chi2*gama/(1+gama)*q1*(((w2+r-p1*cbar)^(1/gama+1))/((w2*(1-tau*phi)+r-p1*cbar)^(1/gama))-(w2*(1-tau*phi)+r-p1*cbar))-tau*w2*r*L; %land rents                                                                                                                                                                                                                                                     
%F(5)=q1-q1*lambda-q1*phi-gama/chi1*(w1+r-p1*cbar)*L1-q1*S1; %land market clearing                                                                          
%F(6)=tau*w2*L1*(1-gama)*(1-nu)*(w1+r-p1*cbar)...
%    +chi2*(1-gama)*(1-nu)/(1+gama)*q1*(((w2+r-p1*cbar)^(1/gama+1))/((w2*(1-tau*phi)+r-p1*cbar)^(1/gama))-(w2*(1-tau*phi)+r-p1*cbar))...
%    +tau*chi2*q1/(1+gama)*(((w2+r-p1*cbar)^(1/gama+1))/((w2*(1-tau*phi)+r-p1*cbar)^(1/gama))-(w2*(1-tau*phi)+r-p1*cbar))-tau*phi*chi2*q1...
%    -tau*w2*theta2*L2;
