N=1000;
ind=[1:N];
I=0.1+exp(-N./ind)*exp(1)*0.1;

tau_m=10;

T=tau_m*log(tau_m*I./(tau_m*I-1));
f=1000./T;

subplot(211);
plot(I,f,'-k','Linewidth',2);
hold on;
plot([0,0.1],[0,0],'-k','Linewidth',4);
hold off;
set(gca,'Fontsize',16);
xlabel('$I$','Fontsize',20);
ylabel('$f$','Fontsize',20);
axis([0,0.2,0,150]);
shg;

