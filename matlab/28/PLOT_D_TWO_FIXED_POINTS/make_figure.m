psi=[0:300]/300;

subplot(111);
D=@(psi) g_0(psi)-g_0(1-psi);

plot(psi,D(psi),'-k','Linewidth',5);
set(gca,'Fontsize',24);
xlabel('$\psi$','Fontsize',32);
ylabel('$D(\psi)$','Fontsize',32);
axis([0,1,-0.1,0.1]);
axis('square');
hold on;

c=0.08;
plot([0,1],[c,c],'--r','Linewidth',3);
text(-0.075,0.08,'$c$','color','red','Fontsize',32);


for i=1:100,
    psi=0.5+0.5*i/100;
    if D(psi)>c,
        psi_0=psi;
    end;
end;

psi_left=0.5;
psi_right=psi_0;
while psi_right-psi_left>10^(-12),
    psi_c=(psi_left+psi_right)/2;
    if D(psi_c)>c,
        psi_right=psi_c;
    else
        psi_left=psi_c;
    end;
end;
psi=(psi_left+psi_right)/2;
plot([psi,psi],[-0.1,D(psi)],'--r','Linewidth',2);
plot(psi,-0.1,'or','Markersize',15,'Markerfacecolor','w');

psi_left=psi_0;
psi_right=1;
while psi_right-psi_left>10^(-12),
    psi_c=(psi_left+psi_right)/2;
    if D(psi_c)>c,
        psi_left=psi_c;
    else
        psi_right=psi_c;
    end;
end;
psi=(psi_left+psi_right)/2;
plot([psi,psi],[-0.1,D(psi)],'--r','Linewidth',2);
plot(psi,-0.1,'or','Markersize',15,'Markerfacecolor','r');

hold off;
shg;