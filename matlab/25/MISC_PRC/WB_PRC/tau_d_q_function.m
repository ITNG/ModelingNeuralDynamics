function tau_d_q_function=tau_d_q_function(tau_d,tau_r,tau_hat);

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

