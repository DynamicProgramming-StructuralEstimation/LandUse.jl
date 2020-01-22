function F=w_f(d,w2)
%computes urban worker income net of commuting costs at a given location

global gama nu cbar sbar theta2 theta1 alpha lambda tau chi1 chi2 L epsilon eta psi sigma

F=w2*(1-tau*d);