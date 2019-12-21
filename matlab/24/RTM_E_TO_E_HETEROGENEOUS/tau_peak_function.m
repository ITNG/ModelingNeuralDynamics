function tau_peak_function=tau_peak_function(tau_d,tau_r,tau_d_q);

dt=0.01;
dt05=dt/2;

s=0;
t=0;
s_inc=exp(-t/tau_d_q)*(1-s)/tau_r-s*tau_d;
while s_inc>0,
    t_old=t;
    s_inc_old=s_inc;
    s_tmp=s+dt05*s_inc;
    s_inc_tmp=exp(-(t+dt05)/tau_d_q)*(1-s_tmp)/tau_r-s_tmp/tau_d;
    s=s+dt*s_inc_tmp;
    t=t+dt;
    s_inc=exp(-t/tau_d_q)*(1-s)/tau_r-s/tau_d;
end;
tau_peak_function=(t_old*(-s_inc)+t*s_inc_old)/(s_inc_old-s_inc);