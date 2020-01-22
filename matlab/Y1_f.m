function F=Y1_f(L1,S1)
% Value of marginal product of labor in urban sector

global gama nu cbar sbar theta2 theta1 alpha lambda tau chi1 chi2 L epsilon eta psi sigma

F=theta1*(alpha*L1.^((sigma-1)/sigma)+(1-alpha)*S1.^((sigma-1)/sigma)).^(sigma/(sigma-1)); 