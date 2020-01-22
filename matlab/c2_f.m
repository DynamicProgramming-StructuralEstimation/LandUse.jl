function F=c2_f(w,r,p1)
%computes individual consumption of urban good at a given location (where earnings net of commuting costs are wd and housing price is qhd).
%note: Cobb-Douglas utility function. 

global gama nu cbar sbar theta2 theta1 alpha lambda tau chi1 chi2 L epsilon eta psi sigma

F=(1-nu)*(1-gama)*(w+r-p1*cbar+sbar)-sbar;
