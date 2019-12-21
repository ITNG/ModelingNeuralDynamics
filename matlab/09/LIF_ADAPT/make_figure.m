clear;
tau_m=10;
I=0.13;
tau_a=40;
delta=0.05;

t_final=300;
dt=0.01;
dt05=dt/2;


subplot(211);
k=1;
v(1)=0;
w(1)=0;
t(1)=0;

while t(k)<t_final,
    v_inc=-v(k)/tau_m+I-w(k)*v(k);
    w_inc=-w(k)/tau_a;
    v_tmp=v(k)+dt05*v_inc;
    w_tmp=w(k)+dt05*w_inc;
    v_inc=-v_tmp/tau_m+I-w_tmp*v_tmp;
    w_inc=-w_tmp/tau_a;
    v_new=v(k)+dt*v_inc;
    w_new=w(k)+dt*w_inc;
    if v_new<=1,
        v(k+1)=v_new;
        w(k+1)=w_new;
        t(k+1)=t(k)+dt;
        k=k+1;
    else
        t_old=t(k);
        t_new=t_old+dt;
        t_spike=(v_new-1)*t_old+(1-v(k))*t_new;
        t_spike=t_spike/(v_new-v(k));
        t(k+1)=t_spike;
        v(k+1)=1;
        w(k+1)=((v_new-1)*w(k)+(1-v(k))*w_new)/(v_new-v(k));
        k=k+1;
        subplot(211);
        plot(t,v,'-k','Linewidth',2);
        hold on;
        plot([t(k),t(k)],[0,1],'--k','Linewidth',1)
        clear v;
        subplot(212);
        plot(t,w,'-k','Linewidth',2);
        hold on;
        plot([t(k),t(k)],[w(k),w(k)+delta],'--k','Linewidth',1);
        w_now=w(k);
        clear w;
        t_now=t(k);
        clear t;
        v(1)=0;
        w(1)=w_now+delta;
        t(1)=t_now;
        k=1;   
    end;
end;
subplot(211);
plot(t,v,'-k','Linewidth',2);
hold off;
set(gca,'Fontsize',16);
axis([0,t_final,0,2]);
ylabel('$v$','Fontsize',20);
axis([0,t_final,0,4]);
set(gca,'Ytick',[0,1]);
subplot(212);
plot(t,w,'-k','Linewidth',2);
set(gca,'Fontsize',16);
axis([0,t_final,0,0.1]);
xlabel('$t$','Fontsize',20); ylabel('$w$','Fontsize',20);
hold off;

shg;
