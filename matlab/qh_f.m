function F=qh_f(q1)
%computes the housing rental price at a location given the land price in that location
%note: Cobb-Douglas housing production function assumed.

global gama nu cbar sbar theta2 theta1 alpha lambda tau chi1 chi2 L epsilon eta psi sigma

F=((1+epsilon)*q1)^(1/(1+epsilon)); %rural housing price