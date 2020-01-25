plot(tvec,L1./(L1+L2))
legend('Agr Emp / Emp')
figure
plot(tvec,L1./L2)
legend('Agr Emp / Urb. Emp')
figure
%plot(thetavec,VAagr)
%legend('Agr VA share')
%figure
plot(tvec,phi)
legend('city size')
figure
plot(tvec,w2./w1)
legend('AGP')
figure
plot(tvec,r./P)
legend('land rents')
figure
plot(tvec,r./w2)
legend('land rents over urban wage')
figure
plot(tvec,q1)
legend('farm land price')
figure
plot(tvec,q1./w2)
legend('farm land price over urban wage')
figure
plot(tvec,p1)
legend('agr good relative price')
figure
plot(tvec,L2./phi)
legend('Urban density')
