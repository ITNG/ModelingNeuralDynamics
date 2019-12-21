function tau_d_q_function=tau_d_q_function(tau_d,tau_r,tau_hat)

tau_d_q_left=1;
while tau_peak_function(tau_d,tau_r,tau_d_q_left)>tau_hat,
    tau_d_q_left=tau_d_q_left/2;
end;


tau_d_q_right=tau_r;
while tau_peak_function(tau_d,tau_r,tau_d_q_right)<tau_hat,
    tau_d_q_right=tau_d_q_right*2;
end;

while tau_d_q_right-tau_d_q_left>10^(-12),
    tau_d_q_mid=(tau_d_q_left+tau_d_q_right)/2;
    if tau_peak_function(tau_d,tau_r,tau_d_q_mid)<=tau_hat,
        tau_d_q_left=tau_d_q_mid;
    else
        tau_d_q_right=tau_d_q_mid;
    end;
end;

tau_d_q_function=(tau_d_q_left+tau_d_q_right)/2;

%..........................................................................

function tau_peak_function=tau_peak_function(tau_d,tau_r,tau_d_q)

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