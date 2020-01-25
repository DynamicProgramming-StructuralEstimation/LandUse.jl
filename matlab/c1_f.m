function F=c1_f(w,r,p1)
%computes individual consumption of rural good at a given location (where earnings net of commuting costs are wd and housing price is qhd).
%note: Cobb-Douglas utility function. 

global gama nu cbar sbar theta2 theta1 alpha lambda tau chi1 chi2 L epsilon eta psi sigma

F=nu*(1-gama)/p1*(w+r-p1*cbar+sbar)+cbar;
