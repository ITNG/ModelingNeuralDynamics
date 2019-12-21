tau_m=10;
I=0.12;
tau_w=100;
delta=0.01;
z_max=0.05;
period=tau_m*log(tau_m*I/(tau_m*I-1));
frequency=1000/period

dt=0.01;
dt05=dt/2;

N=100;

v=zeros(N+1,1);
w=[0:N]'/N*z_max;
phi=zeros(N,1);
t=0;

done=zeros(N+1,1);

while min(done)==0,
    t_old=t;
    v_old=v;
    w_old=w;
    v_inc=-v/tau_m+I-w.*v;
    w_inc=-w/tau_w;
    v_tmp=v+dt05*v_inc;
    w_tmp=w+dt05*w_inc;
    v_inc=-v_tmp/tau_m+I-w_tmp.*v_tmp;
    w_inc=-w_tmp/tau_w;
    v=v+dt*v_inc;
    w=w+dt*w_inc;
    t=t+dt;
    ind=find(v>1 & done==0);
    done(ind)=1;
    phi(ind)=(v(ind)-1).*w_old(ind)+(1-v_old(ind)).*w(ind);
    phi(ind)=phi(ind)./(v(ind)-v_old(ind))+delta;
end;

subplot(111);
z=[0:N]'/N*z_max;
plot(z,phi,'-k','Linewidth',2);
set(gca,'Fontsize',24);
axis([0,z_max,0,z_max]);
axis('square');
hold on;
plot([0,z_max],[0,z_max],'--k');
set(gca,'Xtick',[0,0.02,0.04]);
set(gca,'Ytick',[0,0.02,0.04]);
plot([0.02791,0.02791],[0,0.02791],':k','Linewidth',1)
plot(0.02791,0,'.r','Markersize',40);
text(0.0265,-0.003,'$z_\ast$','Fontsize',32,'Color','r');
hold off;
xlabel('$z$','Fontsize',32); ylabel('$\phi (z)$','Fontsize',32);
shg;

