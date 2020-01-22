function F=EqsysCD0_f(x)
% Equilibrium equations of simplified model without city structure (no commuting costs)

global gama nu cbar sbar theta2 theta1 alpha lambda tau chi1 chi2 L epsilon eta psi sigma

q1=x(1); %land price in rural sector
L1=x(2); %employment in rural sector

L2=L-L1; %employment in urban sector
w2=theta2; %wage rate urban sector with no commuting costs
w1=w2; %wage rate rural sector
S1=(1-alpha)/alpha*w1/q1*L1; %land used as farm sector input
r=1/L*q1*(1-lambda); %land rents
p1=w1/(alpha*theta1)*(L1/S1)^(1-alpha); %relative price of rural-sector good
phi=(L2/chi2)*(gama*(w2+r-p1*cbar)/q1); %city radius
S1h=(L1/chi1)*(gama*(w1+r-p1*cbar)/q1);%rural land used for residential purposes

F(1)=(1-gama)*(1-nu)*(w1+r-p1*cbar)-L2*theta2;
F(2)=S1+S1h+phi-(1-lambda);
