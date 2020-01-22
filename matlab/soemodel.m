        x0=[q1(i,j),phi(i,j),r(i,j),L1(i,j),p1(i,j),S1(i,j)];    
        [x,Fval,exitflag]=fsolve('EqsysCD_pfixed_f',x0,options);
        if exitflag<1;display('exitflag negative');return;end%if exitflag~=1;display('exitflag different from 1');return;end
        q1b(i,j)=x(1); %land price in rural sector
        phib(i,j)=x(2); %city radius
        rb(i,j)=x(3); %land rents
        L1b(i,j)=x(4); %employment in rural sector
        p1b(i,j)=x(5); %relative price of rural-sector good
        S1b(i,j)=x(6); %land used as farm sector input
        
        L2b(i,j)=L-L1b(i,j); %employment in rural sector
        w2b(i,j)=theta2*L2b(i,j)^eta; %wage rate urban sector with no commuting costs
        w1b(i,j)=w2b(i,j)*(1-tau*phib(i,j)); %wage rate rural sector
%        S1h(i,j)=gama/(q1(i,j)*chi1)*(w1(i,j)+r(i,j)-p1(i,j)*cbar)*L1(i,j); %rural land used for residential purposes (housing plus wasted when chi1<1).
%        c1(i,j)=nu*(1-gama)*(w2(i,j)+r(i,j)-p1(i,j)*cbar)/p1(i,j)+cbar;
%        c2(i,j)=(1-nu)*(1-gama)*(w2(i,j)+r(i,j)-p1(i,j)*cbar);
        Db_fringe(i,j)=(chi2*q1b(i,j)/gama2);
        Db_m(i,j)=(chi2*q1b(i,j)/gama2)*((w2b(i,j)*(1-tau*phim)+rb(i,j)-p1b(i,j)*cbar))^(1/gama2-1)/((w2b(i,j)*(1-tau*phib(i,j))+rb(i,j)-p1b(i,j)*cbar)^(1/gama2));       
        Db_cbd(i,j)=(chi2*q1b(i,j)/gama2)*((w2b(i,j)+rb(i,j)-p1b(i,j)*cbar))^(1/gama2-1)/((w2b(i,j)*(1-tau*phib(i,j))+rb(i,j)-p1b(i,j)*cbar)^(1/gama2));
        