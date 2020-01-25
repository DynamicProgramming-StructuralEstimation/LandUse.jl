% Cobb-Douglas preferences, flexible housing elasticity, flexible commuting costs
tic
clc
clear all

global gama nu cbar sbar theta2 theta1 alpha lambda tau chi1 chi2 L epsilon eta psi sigma

paramsCD


GC=1; %set GC=1 to run the flexible model (location-specific epsilon, location-specific tau)

% J=size(thetavec,2); %number of simulations (e.g. time-varying theta; time-varying theta,tau; time-varying theta,tau,L)
J = 1;
for j=1:J %length(epsilonvec) %Uncomment in order to compare the results for different values of sigma (or any other paramter)
    
    for i=1:length(thetavec)
%        [i,j]
       
       epsilon=epsilonvec(j);%sigma=sigmavec(j); %gama=gamavec(j); %tau=tauvec(j); %nu=nuvec(j); %chi2=chi2(j); %eta=etavec(j); %epsilon=epsilonvec(j); %sbar=sbarvec(j); %L=Lvec(j); %u=uvec(j);
    
       theta2=thetavec(i,j);  %theta2=(1+u*(Lvec(i)^(0.1)-1)); %theta1=(0.2+(1-u)*(thetavec(i)-0.2));%theta2=(0.2+u*(thetavec(i)-0.2));
       theta1=theta2;       
       tau=tauvec(i,j);
       L=Lvec(i,j);

       %%SPRING 2019 MODELS (constant epsilon)       
       options=optimoptions('fsolve','Display','none','Algorithm','trust-region-dogleg');%options=optimoptions('fsolve','Display','none','FunValCheck','on','TolFun',1e-20,'MaxFunEvals',1e+3,'MaxIter',1e+3,'StepTolerance',1e-10);%options=optimset('Display','off','FunValCheck','on','TolFun',1e-12,'MaxFunEvals',1e+30,'MaxIter',1e+30);   %'trust-region-dogleg', 'trust-region', or 'levenberg-marquardt'.
       if i==1
           stmodel; %structural change model (no commuting costs);
           x00=[q100(i,j),0.00005,r00(i,j),L100(i,j),p100(i,j),S100(i,j)];%initial guess from stmodel
           [x0,Fval0,exitflag0(i,j)]=fsolve('EqsysCD_f',x00,options);Fvalmax0(i,j)=max(Fval0);
           if exitflag0(i,j)~=1 && Fvalmax0(i,j)>1e-5;
               tau=0.5*(tauvec(1)+tauvec(end));[x0,Fval0,exitflag0(i,j)]=fsolve('EqsysCD_f',x00,options);Fvalmax0(i,j)=max(Fval0);
               tau=tauvec(1);[x0,Fval0,exitflag0(i,j)]=fsolve('EqsysCD_f',x0,options);Fvalmax0(i,j)=max(Fval0);
               if exitflag0(i,j)~=1 && Fvalmax0(i,j)>1e-5;display('exitflag negative');return;end
           end
       else
           [x0,Fval0,exitflag0(i,j)]=fsolve('EqsysCD_f',x0,options);Fvalmax0(i,j)=max(Fval0);if exitflag0(i,j)~=1 && Fvalmax0(i,j)>1e-5;display('exitflag CD negative');return;end
       end
       if isreal(x0)==0;display('Complex number in x0');break;end
       q10(i,j)=x0(1); phi0(i,j)=x0(2); r0(i,j)=x0(3); L10(i,j)=x0(4); p10(i,j)=x0(5); S10(i,j)=x0(6);
       othervars0 %computes other variables
       %soemodel %small-open economy model (p1 exogenous); vars denoted by b
       
       %urbanmodel %"urban economics" model (q1 exogenous, p1 exogenous, S1 exogenous);vars denoted by c
    end  
    display('starting values or each period:')
    [q10 phi0 r0 L10 p10 S10]
    
    if GC==1 %%FLEXIBLE MODEL (Fall 2019)
        for i=1:length(thetavec)
            fprintf('period: %d\n',i)
            
            theta2=thetavec(i,j);  %theta2=(1+u*(Lvec(i)^(0.1)-1)); %theta1=(0.2+(1-u)*(thetavec(i)-0.2));%theta2=(0.2+u*(thetavec(i)-0.2));
            theta1=theta2;
            tau=tauvec(i,j);
            L=Lvec(i,j);
            
            if i==1
                fprintf('starting epsilon search\n')
                eps1vec=linspace(0,15,10); %change of elasticity wrt to distance
                eps0=epsilonvec(j); %elasticity at the fringe
                for ii=2
                    fprintf('epsilon search step %d\n',ii)
                    x0=[q10(i,j),phi0(i,j),r0(i,j),L10(i,j),p10(i,j),S10(i,j)];%initial guess from "simple" model, first value of theta and last value of epsilon
                    [x,Fval,exitflag(i,j)]=fsolve('EqsysGC_f',x0,options,eps0,eps1vec(ii));Fvalmax(i,j)=max(Fval);if exitflag(i,j)~=1 && Fvalmax(i,j)>1e-5;display('exitflag negative');return;end
                end
                for ii=3:length(eps1vec)
                    fprintf('epsilon search step %d\n',ii)
                    [x,Fval,exitflag(i,j)]=fsolve('EqsysGC_f',x,options,eps0,eps1vec(ii));Fvalmax(i,j)=max(Fval);if exitflag(i,j)~=1 && Fvalmax(i,j)>1e-5;display('exitflag negative');return;end
                end
            else
                it(i,j)=1;
                [x,Fval,exitflag(i,j)]=fsolve('EqsysGC_f',x,options,eps0,eps1vec(end));Fvalmax(i,j)=max(Fval);%if exitflag(i,j)<1 && Fvalmax(i,j)>1e-5;display('exitflag negative');return;end
                while Fvalmax(i,j)>1e-5 && it(i,j)<=10;
                    it(i,j)=it(i,j)+1
                    [x,Fval,exitflag(i,j)]=fsolve('EqsysGC_f',x,options,eps0,eps1vec(end));
                    Fvalmax(i,j)=max(Fval);
                end
                if exitflag(i,j)~=1 && Fvalmax(i,j)>1e-5;display('exitflag negative');return;end
            end
            q1(i,j)=x(1); phi(i,j)=x(2); r(i,j)=x(3); L1(i,j)=x(4); p1(i,j)=x(5); S1(i,j)=x(6);
            othervars %computes other variables
        end
    end
end
% save('allsims')
toc

%j=1;plots0 %main plots for the "simple" model (Spring 2019)
%j=1;plots %main plots for the "flexible" model (Fall 2019)
%j=1; plotscomp %comparison of results of "simple" and "flexible" models
%plotscomp2 %comparison of results of flexible model: baseline vs falling tau
%plotscomp3 %comparison of results of flexible model: baseline vs falling tau+rising L

