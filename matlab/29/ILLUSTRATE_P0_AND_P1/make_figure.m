clear; 
clf;
tau_m=10;
I=0.12;
g=0.15;
tau_I=7;


dt=0.01;
dt05=dt/2;

subplot(211);
v(1)=0;
k=1;
while v(k)<1,
    v_inc=-v(k)/tau_m+I-g*exp(-(k-1)*dt/tau_I)*v(k);
    v_tmp=v(k)+dt05*v_inc;
    v_inc=-v_tmp/tau_m+I-g*exp(-(k-1/2)*dt/tau_I)*v_tmp;
    v(k+1)=v(k)+dt*v_inc;
    k=k+1;
end;
plot([0:k-2]*dt,v(1:k-1),'-k','Linewidth',2);
hold on;
period=(k-2)*dt*(v(k)-1)+(k-1)*dt*(1-v(k-1));
period=period/(v(k)-v(k-1));
plot([(k-2)*dt,period],[v(k-1),1],'-k','Linewidth',2);
plot([period,period],[0,1],':k','Linewidth',2);
text(period-0.5,-0.15,'$P_0$','Fontsize',20);

v(1)=1;
k=1;
while v(k)<=1,
    v_inc=-v(k)/tau_m+I-g*exp(-(k-1)*dt/tau_I)*v(k);
    v_tmp=v(k)+dt05*v_inc;
    v_inc=-v_tmp/tau_m+I-g*exp(-(k-1/2)*dt/tau_I)*v_tmp;
    v(k+1)=v(k)+dt*v_inc;
    k=k+1;
end;
plot([0:k-2]*dt,v(1:k-1),'-k','Linewidth',2);
hold on;
period=(k-2)*dt*(v(k)-1)+(k-1)*dt*(1-v(k-1));
period=period/(v(k)-v(k-1));
plot([(k-2)*dt,period],[v(k-1),1],'-k','Linewidth',2);
plot([period,period],[0,1],':k','Linewidth',2);
text(period-1,-0.15,'$P_1$','Fontsize',20);


hold off;
set(gca,'Fontsize',16);
set(gca,'Xtick',[0:15:30]);
axis([0,30,0,1]);
shg;

