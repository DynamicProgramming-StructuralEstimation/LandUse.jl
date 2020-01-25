function F=epsi_f(d,phi,eps0,eps1)
%computes urban worker income net of commuting costs at a given location

global gama nu cbar sbar theta2 theta1 alpha lambda tau chi1 chi2 L epsilon eta psi sigma

%eps0=5;
%eps1=30;
F=exp(log(eps0)-eps1*log(1+phi-d));