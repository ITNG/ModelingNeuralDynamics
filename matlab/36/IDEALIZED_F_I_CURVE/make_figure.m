subplot(211);

I=(-100:0)/100; f=0*I;
plot(I,f,'-k','Linewidth',4);
hold on;
I=(0:100)/100; f=sqrt(I);
plot(I,f,'-k','Linewidth',2);
hold off;
set(gca,'Fontsize',16);
ylabel('$f$','Fontsize',20);
text(-0.02,-0.1,'$I_c$','Fontsize',20);
axis([-1,1,0,1]);
set(gca,'Xtick',[]); set(gca,'Ytick',[]);
shg;
