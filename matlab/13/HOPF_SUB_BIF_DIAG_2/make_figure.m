I=-1+[0:100]/100*0.75;

subplot(111);
plot(I,zeros(101,1),'-k','Linewidth',2);
hold on;

I=-1/4+[0:100]/100*0.25;
r0=sqrt(1/2-sqrt(1/4+I));
R0=sqrt(1/2+sqrt(1/4+I));
plot(I,zeros(101,1),'-k','Linewidth',2);
plot(I,r0,'--k','Linewidth',2);
plot(I,-r0,'--k','Linewidth',2);
plot(I,R0,'-k','Linewidth',2);
plot(I,-R0,'-k','Linewidth',2);

I=[0:100]/100;
R0=sqrt(1/2+sqrt(1/4+I));
plot(I,zeros(101,1),'--k','Linewidth',2);
plot(I,R0,'-k','Linewidth',2);
plot(I,-R0,'-k','Linewidth',2);

set(gca,'Fontsize',20);
xlabel('$I$','Fontsize',26);
ylabel('oscillation amplitude','Fontsize',26);
hold off;
axis([-1,1,-2,2]); axis('square');
shg;