clear; clf;
f_tilde=@(s) 0.5*(s+1/sqrt(2)).*(1/sqrt(2)-s);
for k=1:101,
    s(k)=-1/sqrt(2)+(k-1)/100*sqrt(2);
    u=f_tilde(s(k))/sqrt(2);
    phi(k)=0.5+s(k)/sqrt(2)-u;
    f(k)=0.5+s(k)/sqrt(2)+u;
end;

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
delta=0.15;
plot([0.7,0.7-delta],[0.7,0.7+delta],'-r','Linewidth',4);
plot([0.7,0.7-delta],[0.7,0.7],'-b','Linewidth',2);
plot([0.7-delta,0.7-delta],[0.7,0.7+delta],'-b','Linewidth',2);
text(0.49,0.75,'$u$','color','b','Fontsize',32);
text(0.59,0.67,'$u$','color','b','Fontsize',32);
text(0.67,0.8,'$\tilde{f}$','color','r','Fontsize',32);
plot(0.7-delta,0.7+delta,'.k','Markersize',70);
text(0.7-delta-0.43,0.7+delta+0.05,'$(\varphi,f(\varphi))$','Fontsize',32);
hold off;
shg;