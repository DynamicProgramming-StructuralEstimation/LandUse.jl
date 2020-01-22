function F=EqsysGC_f(x,eps0,eps1)
% Equilibrium equations when general case

global gama nu cbar sbar theta2 theta1 alpha lambda tau chi1 chi2 L epsilon eta psi sigma

q1=x(1); %land price in rural sector
phi=x(2); %city radius
r=x(3); %land rents
L1=x(4); %employment in rural sector
p1=x(5); %relative price of rural-sector good
S1=x(6); %land used as farm sector input

L2=L-L1; %employment in urban sector
w2=VMPL2_f(L2); %wage rate urban sector with no commuting costs
w1=w1fromw2_f(w2,phi); %wage rate rural sector  (from indifference condition at the fringe)

epsilon=epsi_f(phi,phi,eps0,eps1);
qhr=qh_f(q1); %housing price in rural sector
Shr=L1*hd_f(w1,qhr,r,p1)/Hs_f(qhr);; %land allocated to rural housing
c2r=c2_f(w1,r,p1); %urban-good consumption by rural workers; note: as if urban workers in cbd earning w1
Mhr=mh_f(qhr); %urban good used as housing services in rural housing (in each unit of land)

norder=64;
[gnodes,gweights]=gausslegendre(norder); %Gauss-Legendre quadrature nodes and weights
dnodes=(gnodes+1)*(phi-0)/2+0;
for n=1:norder
    epsilon=epsi_f(dnodes(n),phi,eps0,eps1);
    wd=w_f(dnodes(n),w2); %net earnings in a location
    qhd=qhfromqhr_f(qhr,wd,w1,r,p1); %housing price in a location
    mhousvec(n)=mh_f(qhd); %urban good used as housing services in urban housing in a location
    densvec(n)=Hs_f(qhd)/hd_f(wd,qhd,r,p1); %number of workers living in a urban location
    consuvec(n)=densvec(n)*c2_f(wd,r,p1); %total urban consumption in a urban location
    %commvec(n)=densvec(n)*(w2-wd); %urban good used for commuting in a urban location
    qvec(n)=q_f(qhd); %urban land rents in a location
end
intMhu=(phi-0)/2*ones(1,norder)*(gweights.*mhousvec)'; %total urban good used as housing services in urban housing (integral of materials at each location) 
intDensity=(phi-0)/2*ones(1,norder)*(gweights.*densvec)'; %integral of urban workers (integral of urban density at each location)
intc2u=(phi-0)/2*ones(1,norder)*(gweights.*consuvec)'; %total consumption of urban good in cities (integral at each urban location)
%intcc=(phi-0)/2*ones(1,norder)*(gweights.*commvec)'; %total cost of commuting (integral of commuting costs at each urban location)
intcc=intcc_f(qhr,w1,w2,r,p1,phi,eps0,eps1); %total cost of commuting (integral of commuting costs at each urban location)
intqu=(phi-0)/2*ones(1,norder)*(gweights.*qvec)'; %land rents in cities (integral of urban land rents at each location)

F(1)=w1-VMPL1_f(p1,L1,S1); %rural firm labor optimization

F(2)=q1-VMPS1_f(p1,L1,S1); %rural firm land optimization

F(3)=L2-intDensity;

F(4)=q1*S1+q1*Shr+intqu-r*L; %land rents

F(5)=1-lambda-phi-S1-Shr; %land market clearing   

F(6)=L1*c2r+Shr*Mhr+intc2u+intcc+intMhu-Y2_f(L2); %urban good market clearing condition: here in the second term missing w2 at denominator


