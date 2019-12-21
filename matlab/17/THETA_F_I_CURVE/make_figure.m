N=200;
I=1/2+[0:N]/N*0.1;

f=1000*sqrt(2*I-1)/pi;

subplot(211);
plot(I,f,'-k','Linewidth',2);
hold on;
plot([0,1/2],[0,0],'-k','Linewidth',4);
hold off;
set(gca,'Fontsize',16);
xlabel('$I$','Fontsize',20);
ylabel('$f$','Fontsize',20);
axis([0.4,0.6,0,150]);
shg;

