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

epsilon=0.025;
plot([0.5-epsilon,0.5+epsilon],[0.5+epsilon,0.5-epsilon],'-k','Linewidth',1);
plot([0.7-epsilon,0.7+epsilon],[0.7+epsilon,0.7-epsilon],'-k','Linewidth',1);
h=text(0.55,0.43,'$0$','Fontsize',32); set(h,'Rotation',45);
h=text(0.75,0.63,'$s$','Fontsize',32); set(h,'Rotation',45);
delta=0.09;
plot([0.7,0.7-delta],[0.7,0.7+delta],'-r','Linewidth',4);
text(0.5,0.62,'$\tilde{f}(s)$','Color','red','Rotation',45,'Fontsize',32);
hold off;
shg;