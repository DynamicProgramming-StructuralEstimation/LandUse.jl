function F=Hs_f(qh)
%computes total supply of urban housing at a given location.
%note: housing production function Cobb-Douglas assumed.

global gama nu cbar sbar theta2 theta1 alpha lambda tau chi1 chi2 L epsilon eta psi sigma

%F=chi2^(1/(1+epsilon))*qh^epsilon;
F=qh^epsilon;