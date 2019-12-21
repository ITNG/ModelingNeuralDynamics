v=[-100:50];
B=1./(1+exp(-0.062*v)/3.57);

subplot(211);
plot(v,B);

set(gca,'Fontsize',16);
xlabel('$v_{\rm post}$','Fontsize',20);
ylabel('$B$','Fontsize',20);