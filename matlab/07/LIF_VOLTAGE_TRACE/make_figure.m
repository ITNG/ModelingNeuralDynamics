clear; clf;
tau_m=10;
I=0.11;

t_final=100;
dt=0.01;
dt05=dt/2;


subplot(211);
k=1;
v(1)=0;
t(1)=0;

while t(k)<t_final,
    v_inc=-v(k)/tau_m+I;
    v_tmp=v(k)+dt05*v_inc;
    v_inc=-v_tmp/tau_m+I;
    v_new=v(k)+dt*v_inc;
    if v_new<=1,
        v(k+1)=v_new;
        t(k+1)=t(k)+dt;
        k=k+1;
    else
        t_old=t(k);
        t_new=t_old+dt;
        t_spike=(v_new-1)*t_old+(1-v(k))*t_new;
        t_spike=t_spike/(v_new-v(k));
        t(k+1)=t_spike;
        v(k+1)=1;
        k=k+1;
        plot(t,v,'-k','Linewidth',2);
        hold on;
        plot([t(k),t(k)],[0,1],'--k','Linewidth',1)
        clear v;
        t_now=t(k);
        clear t;
        v(1)=0;
        t(1)=t_now;
        k=1;   
    end;
end;
plot(t,v,'-k','Linewidth',2);
hold off;
set(gca,'Fontsize',16);
axis([0,t_final,0,2]);
xlabel('$t$','Fontsize',20); ylabel('$v$','Fontsize',20);

shg;
