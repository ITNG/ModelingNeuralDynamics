clear; clf; rng('default'); rng(63806);
tic;

% Define network parameters: 

i_ext_e=0; i_ext_i=0;
g_ei=.35; g_ie=0.5; g_ii=0.15; 
v_rev_e=0; v_rev_i=-75; 
tau_r_e=0.5; tau_peak_e=0.5; tau_d_e=3; 
tau_r_i=0.5; tau_peak_i=0.5; tau_d_i=9; 
t_final=200;   % Time (in ms) simulated. 
dt=0.01;       % Time step used in solving the differential equations.

% Process network parameters a bit:

tau_dq_e=tau_d_q_function(tau_d_e,tau_r_e,tau_peak_e);
tau_dq_i=tau_d_q_function(tau_d_i,tau_r_i,tau_peak_i);
dt05=dt/2; m_steps=round(t_final/dt);

% Define periodic stimuli to the E-cells

f_main=40;
f_dist=25; 
T_main=1000/f_main;
T_dist=1000/f_dist;
alpha_main=20;
alpha_dist=0.1; 
A_main=0.5;
A_dist=0.5;
t=(0:m_steps)*dt;
N=3000;
s=(1:N)/N;
den_main=mean(exp(alpha_main*cos(pi*s).^2)-1);
den_dist=mean(exp(alpha_dist*cos(pi*s).^2)-1);
I_main=A_main*(exp(alpha_main*cos(pi*t/T_main).^2)-1)/den_main;
I_dist=A_dist*(exp(alpha_dist*cos(pi*t/T_dist).^2)-1)/den_dist; 

% initialize dynamic variables

v_e=-70; m_e=m_e_inf(v_e); h_e=h_e_inf(v_e); n_e=n_e_inf(v_e); q_e=0; s_e=0;
v_i=-75; m_i=m_i_inf(v_i); h_i=h_i_inf(v_i); n_i=n_i_inf(v_i); q_i=0; s_i=0;

% solve the system of Hodgkin-Huxley-like equations using the midpoint method

num_spikes_e=0; t_e_spikes=[]; 
num_spikes_i=0; t_i_spikes=[]; 

for k=1:m_steps,
    t_old=t(k); t_new=t(k+1); t_mid=(t_old+t_new)/2;
	
	v_e_inc=0.1*(-67-v_e)+80*n_e.^4.*(-100-v_e)+100*m_e.^3.*h_e.*(50-v_e) ...
               +(g_ie'*s_i).*(v_rev_i-v_e) ...
               +i_ext_e+I_main(k)+I_dist(k);
	n_e_inc=(n_e_inf(v_e)-n_e)./tau_n_e(v_e);
    h_e_inc=(h_e_inf(v_e)-h_e)./tau_h_e(v_e);
    q_e_inc=(1+tanh(v_e/10))/2.*(1-q_e)/0.1-q_e./tau_dq_e;
    s_e_inc=q_e.*(1-s_e)./tau_r_e-s_e./tau_d_e;
    
	v_i_inc=0.1*(-65-v_i)+9*n_i.^4.*(-90-v_i)+35*m_i.^3.*h_i.*(55-v_i) ...
               +(g_ei'*s_e).*(v_rev_e-v_i)+(g_ii'*s_i).*(v_rev_i-v_i) ...
               +i_ext_i;
	n_i_inc=(n_i_inf(v_i)-n_i)./tau_n_i(v_i);
    h_i_inc=(h_i_inf(v_i)-h_i)./tau_h_i(v_i);
    q_i_inc=(1+tanh(v_i/10))/2.*(1-q_i)/0.1-q_i./tau_dq_i;
    s_i_inc=q_i.*(1-s_i)./tau_r_i-s_i./tau_d_i;

	v_e_tmp=v_e+dt05*v_e_inc;
	n_e_tmp=n_e+dt05*n_e_inc;
	m_e_tmp=m_e_inf(v_e_tmp);
    h_e_tmp=h_e+dt05*h_e_inc;
    q_e_tmp=q_e+dt05*q_e_inc;   
    s_e_tmp=s_e+dt05*s_e_inc; 
    
	v_i_tmp=v_i+dt05*v_i_inc;
	n_i_tmp=n_i+dt05*n_i_inc;
	m_i_tmp=m_i_inf(v_i_tmp);
    h_i_tmp=h_i+dt05*h_i_inc;
    q_i_tmp=q_i+dt05*q_i_inc;   
    s_i_tmp=s_i+dt05*s_i_inc;    

	v_e_inc=0.1*(-67-v_e_tmp)+80*n_e_tmp.^4.*(-100-v_e_tmp)+100*m_e_tmp.^3.*h_e_tmp.*(50-v_e_tmp) ...
               +(g_ie'*s_i_tmp).*(v_rev_i-v_e_tmp) ...
               +i_ext_e+(I_main(k)+I_main(k+1))/2+(I_dist(k)+I_dist(k+1))/2;
	n_e_inc=(n_e_inf(v_e_tmp)-n_e_tmp)./tau_n_e(v_e_tmp);
    h_e_inc=(h_e_inf(v_e_tmp)-h_e_tmp)./tau_h_e(v_e_tmp);
    q_e_inc=(1+tanh(v_e_tmp/10))/2.*(1-q_e_tmp)/0.1-q_e_tmp./tau_dq_e;
    s_e_inc=q_e_tmp.*(1-s_e_tmp)./tau_r_e-s_e_tmp./tau_d_e;
    
	v_i_inc=0.1*(-65-v_i_tmp)+9*n_i_tmp.^4.*(-90-v_i_tmp)+35*m_i_tmp.^3.*h_i_tmp.*(55-v_i_tmp) ...
               +(g_ei'*s_e_tmp).*(v_rev_e-v_i_tmp)+(g_ii'*s_i_tmp).*(v_rev_i-v_i_tmp) ...
               +i_ext_i;
	n_i_inc=(n_i_inf(v_i_tmp)-n_i_tmp)./tau_n_i(v_i_tmp);
    h_i_inc=(h_i_inf(v_i_tmp)-h_i_tmp)./tau_h_i(v_i_tmp);
    q_i_inc=(1+tanh(v_i_tmp/10))/2.*(1-q_i_tmp)/0.1-q_i_tmp./tau_dq_i;
    s_i_inc=q_i_tmp.*(1-s_i_tmp)./tau_r_i-s_i_tmp./tau_d_i;
        
    v_e_old=v_e;
    v_i_old=v_i;

	v_e=v_e+dt*v_e_inc;
	m_e=m_e_inf(v_e); h_e=h_e+dt*h_e_inc; n_e=n_e+dt*n_e_inc; 
    q_e=q_e+dt*q_e_inc;
    s_e=s_e+dt*s_e_inc;
    
	v_i=v_i+dt*v_i_inc;
    m_i=m_i_inf(v_i); h_i=h_i+dt*h_i_inc; n_i=n_i+dt*n_i_inc; 
    q_i=q_i+dt*q_i_inc;
    s_i=s_i+dt*s_i_inc;
    

	% Determine whether the E- and I-cell spiked in the couurent step
    
    if v_e_old>-20 & v_e<=-20,
        num_spikes_e=num_spikes_e+1;
        t_e_spikes(num_spikes_e)=(t_old*(-20-v_e)+t_new*(v_e_old+20))/(v_e_old-v_e);
    end 
    
    if v_i_old>-20 & v_i<=-20,
        num_spikes_i=num_spikes_i+1;
        t_i_spikes(num_spikes_i)=(t_old*(-20-v_i)+t_new*(v_i_old+20))/(v_i_old-v_i);
    end 

end;

% plot the spike rastergram

rastergram;
toc

