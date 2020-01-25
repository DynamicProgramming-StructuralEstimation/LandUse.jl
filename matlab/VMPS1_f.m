function F=VMPS1(p1,L1,S1)
% Value of marginal product of labor in urban sector

global gama nu cbar sbar theta2 theta1 alpha lambda tau chi1 chi2 L epsilon eta psi sigma

F=p1*(1-alpha)*theta1*S1^(-1/sigma)*(alpha*L1^((sigma-1)/sigma)+(1-alpha)*S1^((sigma-1)/sigma))^(1/(sigma-1)); 