function F=hd_f(w,qh,r,p1)
%computes individual demand of housing at a given location (where earnings net of commuting costs are w and housing price is qh).
%note: Cobb-Douglas utility function assumed. 

global gama nu cbar sbar theta2 theta1 alpha lambda tau chi1 chi2 L epsilon eta psi sigma

F=gama*(w+r-p1*cbar+sbar)/qh;
