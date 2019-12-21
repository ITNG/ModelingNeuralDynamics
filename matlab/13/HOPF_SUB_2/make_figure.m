I=-0.2;
r=[0:100]/100*1.2;
subplot(111);
plot(r,I*r+r.^3-r.^5,'-k','Linewidth',2);
hold on;
I=0;
plot(r,I*r+r.^3-r.^5,'-k','Linewidth',2);
I=0.2;
plot(r,I*r+r.^3-r.^5,'-k','Linewidth',2);
I=-0.4;
plot(r,I*r+r.^3-r.^5,'-k','Linewidth',2);
set(gca,'Fontsize',24);
xlabel('$r$','Fontsize',32); ylabel('$f$','Fontsize',32);
axis([0,1.2,-0.5,0.5]);
set(gca,'Xtick',[0:0.3:1.2]);
plot([0,1.2],[0,0],'--k');
text(0.7,0.23,'$I=0$','Fontsize',28)
text(0.6,0.07,'$I=-0.2$','Fontsize',28)
text(0.8,0.39,'$I=0.2$','Fontsize',28)
text(0.5,-0.06, '$I=-0.4$','Fontsize',28);
set(get(gca,'ylabel'),'position',[-0.07 0 0])
hold off;
shg;