function F=q_f(qh)
%computes rental price of land at a location given the housing price qh at that location
%note: housing production function Cobb-Douglas assumed.

global gama nu cbar sbar theta2 theta1 alpha lambda tau chi1 chi2 L epsilon eta psi sigma

F=(1/(1+epsilon))*qh^(1+epsilon);