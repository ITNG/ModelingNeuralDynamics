clear; clf;
phi=[0:100]/100;
tau_m=2; 
I=0.13;
dv=0.1;



g=atan(tan(pi*phi-pi/2)+dv/sqrt(tau_m*I-1/4))/pi+1/2-phi;


subplot(111);
plot(phi,g,'-k','Linewidth',6);
set(gca,'Fontsize',24);
xlabel('$\varphi$','Fontsize',32);
ylabel('$g$','Fontsize',32);
axis([0,1,0,1]); 
axis('square');
shg;