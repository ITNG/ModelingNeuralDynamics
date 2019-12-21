clear; clf;

v=[-1000:500]/10;

r=(120-v)./(1+exp(-(v+15)/5))*4/25;

subplot(211);
plot(v,r,'-k','Linewidth',2);
set(gca,'Fontsize',16);
xlabel('$v$ [mV]','Fontsize',20);
ylabel('$c_\infty$','Fontsize',20);
shg;