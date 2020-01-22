L2(i,j)=L-L1(i,j); %employment in rural sector
w2(i,j)=VMPL2_f(L2(i,j)); %wage rate urban sector with no commuting costs
w1(i,j)=w1fromw2_f(w2(i,j),phi(i,j)); %wage rate rural sector

epsilon=epsilonvec(j);
T(i,j)=0; %lump-sum transfer from the housing tax revenues
qhr(i,j)=qh_f(q1(i,j)); %housing price in rural sector
hdr(i,j)=hd_f(w1(i,j),qhr(i,j),r(i,j),p1(i,j));
Shr(i,j)=L1(i,j)*hd_f(w1(i,j),qhr(i,j),r(i,j),p1(i,j))/Hs_f(qhr(i,j)); %rural land used for residential purposes (housing plus wasted when chi1<1).

GDP(i,j)=(p1(i,j)*Y1_f(L1(i,j),S1(i,j))+Y2_f(L2(i,j)))/L; %value of GDP per capita
INC(i,j)=GDP(i,j)-intcc_f(qhr(i,j),w1(i,j),w2(i,j),r(i,j),p1(i,j),phi(i,j),eps0,eps1vec(end)); %value of disposable income (GDP net of commuting costs) TO CHECK: should we multiply intcc by L?       
EAR(i,j)=INC(i,j)-r(i,j)*L; %value of average labor earnings  (GDP net of commuting costs minus land rents)

epsilon=epsi_f(phi(i,j),phi(i,j),eps0,eps1vec(end));
qh_fringe(i,j)=qhr(i,j);%price of a unit of housing at the fringe
q_fringe(i,j)=q_f(qh_fringe(i,j)); %price of a unit of land at the fringe
h_fringe(i,j)=hd_f(w_f(phi(i,j),w2(i,j)),qh_fringe(i,j),r(i,j),p1(i,j));
H_fringe(i,j)=Hs_f(qh_fringe(i,j));
D_fringe(i,j)=H_fringe(i,j)/h_fringe(i,j);

epsilon=epsi_f(0,phi(i,j),eps0,eps1vec(end));
qh_cbd(i,j)=qhfromqhr_f(qhr(i,j),w2(i,j),w1(i,j),r(i,j),p1(i,j));%price of a unit of housing at the cbd
q_cbd(i,j)=q_f(qh_cbd(i,j));%price of a unit of land at the cbd
h_cbd(i,j)=hd_f(w2(i,j),qh_cbd(i,j),r(i,j),p1(i,j));
H_cbd(i,j)=Hs_f(qh_cbd(i,j));
D_cbd(i,j)=H_cbd(i,j)/h_cbd(i,j);

c1r(i,j)=c1_f(w1(i,j),r(i,j),p1(i,j));
c2r(i,j)=c2_f(w1(i,j),r(i,j),p1(i,j));
hr(i,j)=hd_f(w1(i,j),qhr(i,j),r(i,j),p1(i,j));
welf(i,j)=u_f(c1r(i,j),c2r(i,j),hr(i,j));
            
if i==1
    P(1,j)=p1(1,j);
else
    Pgr_laspeyres(i,j)=(p1(i,j).*Y1_f(L1(i-1,j),S1(i-1,j))+Y2_f(L2(i-1,j)))./(p1(i-1,j).*Y1_f(L1(i-1,j),S1(i-1,j))+Y2_f(L2(i-1,j)));
    Pgr_paasche(i,j)=(p1(i,j).*Y1_f(L1(i,j),S1(i,j))+Y2_f(L2(i,j)))./(p1(i-1,j).*Y1_f(L1(i,j),S1(i,j))+Y2_f(L2(i,j)));
    Pgr(i,j)=(Pgr_laspeyres(i,j).^0.5).*(Pgr_paasche(i,j).^0.5);
    P(i,j)=P(i-1,j)*Pgr(i,j);
end
