function F=VMPL2_f(L2)
% Value of marginal product of labor in urban sector

global gama nu cbar sbar theta2 theta1 alpha lambda tau chi1 chi2 L epsilon eta psi sigma

F=theta2*L2^eta; %wage rate urban sector with no commuting costs
