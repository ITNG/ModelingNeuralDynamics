function P0=P0(tau_m,J,g,tau_I);

I=J+1/tau_m;

dt=0.01;
dt05=dt/2;

v=0;
k=1;
while v<1,
    v_old=v;
    v_inc=-v/tau_m+I-g*exp(-(k-1)*dt/tau_I)*v;
    v_tmp=v+dt05*v_inc;
    v_inc=-v_tmp/tau_m+I-g*exp(-(k-1/2)*dt/tau_I)*v_tmp;
    v=v+dt*v_inc;
    k=k+1;
end;
period=(k-2)*dt*(v-1)+(k-1)*dt*(1-v_old);
period=period/(v-v_old);


P0=period;