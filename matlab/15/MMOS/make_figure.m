clear; clf;
a=5; tau_n=60; I_ext=-4.2; tau_adapt=150; delta=0.2;

t_final=1000;
dt=0.01; dt05=dt/2;
m_steps=round(t_final/dt);
v=zeros(m_steps+1,1); n=v;


v(1)=-1;
n(1)=-4.75;
I_adapt=-delta;



for k=1:m_steps,
    v_inc=v(k)-v(k)^3/3-n(k)+I_ext+I_adapt;
    n_inc=(a*v(k)-n(k))/tau_n;
    v_tmp=v(k)+dt05*v_inc;
    n_tmp=n(k)+dt05*n_inc;
    v_inc=v_tmp-v_tmp^3/3-n_tmp+I_ext+I_adapt*exp(-dt05/tau_adapt);
    n_inc=(a*v_tmp-n_tmp)/tau_n;
    v(k+1)=v(k)+dt*v_inc;
    n(k+1)=n(k)+dt*n_inc;
    I_adapt=I_adapt*exp(-dt/tau_adapt);
    if v(k+1)<0 && v(k)>=0, 
        I_adapt=I_adapt-delta;
    end;
end;


hold off;

subplot(211);
t=[0:m_steps]*dt;
plot(t,v,'-k','Linewidth',2);
set(gca,'Fontsize',14);
ylabel('$v$','Fontsize',18);
xlabel('$t$','Fontsize',18);
axis([0,t_final,-3,3]);



shg;


