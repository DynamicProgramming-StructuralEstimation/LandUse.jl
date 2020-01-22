function F=w1fromw2_f(w2,phi)
% Value of marginal product of labor in urban sector

global gama nu cbar sbar theta2 theta1 alpha lambda tau chi1 chi2 L epsilon eta psi sigma

F=psi*w2*(1-tau*phi); %wage rate rural sector (from indifference condition at the fringe)
