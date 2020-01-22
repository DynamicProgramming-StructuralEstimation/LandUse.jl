function F=intcc_f(qhr,w1,w2,r,p1,phi,eps0,eps1)
% Total commuting costs (integral of commuting costs at all urban locations)

global gama nu cbar sbar theta2 theta1 alpha lambda tau chi1 chi2 L epsilon eta psi sigma

norder=64;
[gnodes,gweights]=gausslegendre(norder); %Gauss-Legendre quadrature nodes and weights
dnodes=(gnodes+1)*(phi-0)/2+0;
for n=1:norder
    epsilon=epsi_f(dnodes(n),phi,eps0,eps1);
    wd=w_f(dnodes(n),w2); %net earnings in a location
    qhd=qhfromqhr_f(qhr,wd,w1,r,p1); %housing price in a location
    d(n)=Hs_f(qhd)/hd_f(wd,qhd,r,p1);
    commvec(n)=d(n)*(w2-wd); %urban good used for commuting in a urban location
end
F=(phi-0)/2*ones(1,norder)*(gweights.*commvec)';