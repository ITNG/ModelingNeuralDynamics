I=-1;
r=[0:100]/100*1.2;
subplot(111);
plot(r,I*r+r.^3,'-k','Linewidth',2);
hold on;
I=0;
plot(r,I*r+r.^3,'-k','Linewidth',2);
I=1;
plot(r,I*r+r.^3,'-k','Linewidth',2);
set(gca,'Fontsize',24);
xlabel('$r$','Fontsize',32); ylabel('$f$','Fontsize',32);
axis([0,1.2,-1,3]);
set(gca,'Xtick',[0:0.3:1.2]);
plot([0,1.2],[0,0],'--k');
text(0.55,1.4,'$I=1$','Fontsize',32)
text(0.95,0.70,'$I=0$','Fontsize',32)
text(0.7,-0.6,'$I=-1$','Fontsize',32)
hold off;
shg;