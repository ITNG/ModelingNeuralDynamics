clear; clf;

varphi=(0:1000)/1000;

subplot(111);
plot(varphi,g(varphi),'-k','Linewidth',5);
set(gca,'Fontsize',24);
xlabel('$\varphi$','Fontsize',32);
ylabel('$g(\varphi)$','Fontsize',32);
axis([0,1,-1,0]); axis('square');
shg;