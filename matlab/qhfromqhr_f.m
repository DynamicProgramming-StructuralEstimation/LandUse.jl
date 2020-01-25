function F=qhfromqhr_f(qhr,wd,w1,r,p1)
%computes urban housing price at a given location (where earnings net of commuting costs are wd) given the rural housing price
%note: utility function Cobb-Douglas assumed.

global gama nu cbar sbar theta2 theta1 alpha lambda tau chi1 chi2 L epsilon eta psi sigma

F=qhr*((wd+r-p1*cbar+sbar)/(w1+r-p1*cbar+sbar))^(1/gama);