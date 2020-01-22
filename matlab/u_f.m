function F=u_f(c1,c2,h)
%Computes welfare

global gama nu cbar sbar theta2 theta1 alpha lambda tau chi1 chi2 L epsilon eta psi sigma

F=(c1-cbar)^(nu*(1-gama))*(c2+sbar)^((1-nu)*(1-gama))*h^gama;