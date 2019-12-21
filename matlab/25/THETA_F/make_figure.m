clear; clf;
phi=[0:100]/100;
tau_m=2; 
I=0.13;
dv=0.1;

f=atan(tan(pi*(phi)-pi/2)+dv/sqrt(tau_m*I-1/4))/pi+1/2;

subplot(111);
plot(phi,f,'-k','Linewidth',6);
set(gca,'Fontsize',24);
xlabel('$\varphi$','Fontsize',32);
ylabel('$f$','Fontsize',32);
axis([0,1,0,1]); 
axis('square');
hold on;
plot([0,1],[0,1],'--k','Linewidth',2);
hold off;
shg;