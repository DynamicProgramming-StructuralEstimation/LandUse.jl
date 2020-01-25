function F=mh_f(qh)
%computes materials (urban good) used for housing services given the housing price qh.
%note: Cobb-Douglas housing production function assumed. 

global gama nu cbar sbar theta2 theta1 alpha lambda tau chi1 chi2 L epsilon eta psi sigma

F=(epsilon/(1+epsilon))*(Hs_f(qh))^((1+epsilon)/epsilon);
