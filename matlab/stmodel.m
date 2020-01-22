        x000=[1,0.5];
        [x00,Fval00,exitflag00]=fsolve('EqsysCD0_f',x000,options);
        %if exitflag~=1;display('exitflag different from 1');return;end
        q100(i,j)=x00(1); %land price in rural sector
        L100(i,j)=x00(2); %employment in rural sector
        L200(i,j)=L-L100(i,j); %employment in urban sector
        w200(i,j)=theta2; %wage rate urban sector with no commuting costs
        w100(i,j)=w200(i,j); %wage rate rural sector
        S100(i,j)=(1-alpha)/alpha*w100(i,j)/q100(i,j)*L100(i,j); %farm land
        r00(i,j)=1/L*q100(i,j)*(1-lambda); %land rents
        p100(i,j)=w100(i,j)/(alpha*theta1)*(L100(i,j)/S100(i,j))^(1-alpha); %farm good relative price
        phi00(i,j)=(L200(i,j)/chi2)*(gama*(w200(i,j)+r00(i,j)-p100(i,j)*cbar)/q100(i,j)); %city size
        