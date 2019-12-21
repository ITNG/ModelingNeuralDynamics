clear;  clf;

tau_m=10; I=0.12;
T=tau_m*log(tau_m*I/(tau_m*I-1)); % firing period without inhibition


N=10;

dt=0.001;
dt05=dt/2;

subplot(311);
for i=1:N,
    v0=(1-exp(-(N-i)*T/(N*tau_m)))*tau_m*I;
    v(1)=v0;
    k=1;
    while v(k)<1,
        v_inc=-v(k)/tau_m+I;
        v_tmp=v(k)+dt05*v_inc;
        v_inc=-v_tmp/tau_m+I;
        v(k+1)=v(k)+dt*v_inc;
        k=k+1;
    end;
    plot([0:k-2]*dt,v(1:k-1),'-k','Linewidth',1);
    hold on;
    period=(k-2)*dt*(v(k)-1)+(k-1)*dt*(1-v(k-1));
    period=period/(v(k)-v(k-1));
    plot([(k-2)*dt,period],[v(k-1),1],'-k','Linewidth',1);
    plot([period,period],[0,1],':k','Linewidth',1);
end;
hold off;
set(gca,'Fontsize',14);
axis([0,30,0,1]);
shg;

g_syn=0.15; tau_I=9;

subplot(312);
for i=1:N,
    v0=(1-exp(-(N-i)*T/(N*tau_m)))*tau_m*I;
    v(1)=v0;
    k=1;
    while v(k)<1,
        v_inc=-v(k)/tau_m+I-g_syn*exp(-(k-1)*dt/tau_I)*v(k);
        v_tmp=v(k)+dt05*v_inc;
        v_inc=-v_tmp/tau_m+I-g_syn*exp(-(k-1/2)*dt/tau_I)*v_tmp;
        v(k+1)=v(k)+dt*v_inc;
        k=k+1;
    end;
    plot([0:k-2]*dt,v(1:k-1),'-k','Linewidth',1);
    hold on;
    period=(k-2)*dt*(v(k)-1)+(k-1)*dt*(1-v(k-1));
    period=period/(v(k)-v(k-1));
    plot([(k-2)*dt,period],[v(k-1),1],'-k','Linewidth',1);
    plot([period,period],[0,1],':k','Linewidth',1);
end;
hold off;
set(gca,'Fontsize',14);
axis([0,30,0,1]);
shg;

g_syn=2; tau_I=1;

subplot(313);
for i=1:N,
    v0=(1-exp(-(N-i)*T/(N*tau_m)))*tau_m*I;
    v(1)=v0;
    k=1;
    while v(k)<1,
        v_inc=-v(k)/tau_m+I-g_syn*exp(-(k-1)*dt/tau_I)*v(k);
        v_tmp=v(k)+dt05*v_inc;
        v_inc=-v_tmp/tau_m+I-g_syn*exp(-(k-1/2)*dt/tau_I)*v_tmp;
        v(k+1)=v(k)+dt*v_inc;
        k=k+1;
    end;
    plot([0:k-2]*dt,v(1:k-1),'-k','Linewidth',1);
    hold on;
    period=(k-2)*dt*(v(k)-1)+(k-1)*dt*(1-v(k-1));
    period=period/(v(k)-v(k-1));
    plot([(k-2)*dt,period],[v(k-1),1],'-k','Linewidth',1);
    plot([period,period],[0,1],':k','Linewidth',1);
end;
hold off;
set(gca,'Fontsize',14);
axis([0,30,0,1]);
xlabel('$t$','Fontsize',18)
shg;