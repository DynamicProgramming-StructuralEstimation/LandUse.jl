                                                                                               
tvec=[1855;1870;1885;1900;1915;1930;1945;1960;1975;1990;2005;2020;2035;2050];

thetavec=[[0.32, 0.33, 0.34, 0.36, 0.38, 0.41, 0.48, 0.7, 1.35, 2.3, 3, 4.5, 5, 5.5]',...
    [0.32, 0.33, 0.34, 0.36, 0.38, 0.41, 0.48, 0.7, 1.35, 2.3, 3, 4.5, 5, 5.5]',...
    [0.32, 0.33, 0.34, 0.36, 0.38, 0.41, 0.48, 0.7, 1.35, 2.3, 3, 4.5, 5, 5.5]'];
    %linspace(0.25,5,10);%thetavec=[0.24, 0.28, 0.35, 0.45, 0.55, 0.85, 1.5, 2.25, 3, 5, 7, 8];%
%thetavec=[0.32, 0.33, 0.35, 0.37, 0.40, 0.50, 0.85, 1.65, 2.5, 3, 4.5, 5, 5.5];
%theta2=0.25; %urban sector TFP
%theta1=0.25;%farm sector TFP

%tau=8; %commuting costs
tauvec=[ones(length(tvec),1)*7,...
    [14, 13.5, 13, 12.5, 11.5, 10.5, 9.5, 8.5, 7.5, 7.4, 7.2, 7.1, 7, 7]',...
    [14, 13.5, 13, 12.5, 11.5, 10.5, 9.5, 8.5, 7.5, 7.4, 7.2, 7.1, 7, 7]'];
    %[23, 23, 22, 20, 16, 14, 12, 9 , 7.5, 7.4, 7.2, 7.1, 7, 7]'];%tauvec=linspace(8,2,10); 

epsconst=4; %housing supply elasticity
epsilonvec=epsconst*ones(1,3);
%epsilonvec=[2,4]; 

gama=1/4; %housing weight in utility function
%gamavec=[1/5,1/2]; 

sigma=0.99; %land-labor elasticity of substitution in farm production function
%sigmavec=[0.99;0.25]; 

nu=0.015; %weight of rural good consumption on consumption composite
%nuvec=[0.05,0.15]; 

eta=0; %agglomeration forces
%etavec=[0;0.04]; 

cbar=0.2; %agr cons subsistence level

sbar=0; %neg urba subsistence level
%sbarvec=[0;0.1]; 

%L=1; %total population (note: total land normalized to 1) 
%Lvec=[1, 1.01, 1.05, 1.20, 1.33, 1.36, 1.32, 1.26, 1.40, 1.61, 1.74, 1.86, 1.91, 1.92];
%Lvec=[1, 1.05, 1.20, 1.33, 1.36, 1.32, 1.26, 1.40, 1.61, 1.74, 1.86, 1.91, 1.92];
Lvec=[ones(length(tvec),1),ones(length(tvec),1),exp(linspace(log(1),log(1.92),14))']; %Note: make Lvec(1)=100 to get meaningful densities (i.e. in terms of people/km2)

alpha=0.75;%labor weight on farm sector production function

lambda=0;%useless land: non-farm, non-urban land (forests, national parks...)

psi=1; %urban ammenities relative to rural ammenities (note: creates wedge in urban-rural wage)
%psivec=linspace(1,2,10); 

chi1=1; %"easiness" to convert land into housing in rural sector
chi2=1; %"easiness" to convert land into housing in urban sector
    if chi1~=1||chi2~=1;display('need adjustments in the code when chi1 and chi2 are different than one');return;end
