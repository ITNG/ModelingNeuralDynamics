clear; clf; rng('default'); rng(63806);
tic;

% Define network parameters: 

i_ext_e=1.4; i_ext_i=0;
g_ei=0.25; g_ie=0.25;
v_rev_e=0; v_rev_i=-75; 
tau_r_e=0.5; tau_peak_e=0.5; tau_d_e=3; 
tau_r_i=0.5; tau_peak_i=0.5; tau_d_i=9; 
t_final=200;   % Time (in ms) simulated. 
dt=0.01;       % Time step used in solving the differential equations.


% Process network parameters a bit:

tau_dq_e=tau_d_q_function(tau_d_e,tau_r_e,tau_peak_e);
tau_dq_i=tau_d_q_function(tau_d_i,tau_r_i,tau_peak_i);
dt05=dt/2; m_steps=round(t_final/dt);

% initialize dynamic variables

v_e=zeros(m_steps+1,1); v_i=zeros(m_steps+1,1);
v_e(1)=-75; v_i(1)=-75;
h_e=0.1; n_e=0.1; q_e=0; s_e=0; 
h_i=0.1; n_i=0.1; q_i=0; s_i=0;


% solve the system of Hodgkin-Huxley-like equations using the midpoint method

num_spikes_e=0; 
num_spikes_i=0; 

for k=1:m_steps,
    t_old=(k-1)*dt; t_new=k*dt; t_mid=(t_old+t_new)/2;
	
	v_e_inc=0.1*(-67-v_e(k))+80*n_e.^4.*(-100-v_e(k)) ...
           +100*m_e_inf(v_e(k)).^3.*h_e.*(50-v_e(k)) ...
           +(g_ie*s_i).*(v_rev_i-v_e(k)) ...
           +i_ext_e;
	n_e_inc=(n_e_inf(v_e(k))-n_e)./tau_n_e(v_e(k));
    h_e_inc=(h_e_inf(v_e(k))-h_e)./tau_h_e(v_e(k));
    q_e_inc=(1+tanh(v_e(k)/10))/2.*(1-q_e)/0.1-q_e./tau_dq_e;
    s_e_inc=q_e.*(1-s_e)./tau_r_e-s_e./tau_d_e;
	v_i_inc=0.1*(-65-v_i(k))+9*n_i.^4.*(-90-v_i(k)) ...
           +35*m_i_inf(v_i(k)).^3.*h_i.*(55-v_i(k)) ...
           +(g_ei*s_e).*(v_rev_e-v_i(k)) ...
           +i_ext_i;
	n_i_inc=(n_i_inf(v_i(k))-n_i)./tau_n_i(v_i(k));
    h_i_inc=(h_i_inf(v_i(k))-h_i)./tau_h_i(v_i(k));
    q_i_inc=(1+tanh(v_i(k)/10))/2.*(1-q_i)/0.1-q_i./tau_dq_i;
    s_i_inc=q_i.*(1-s_i)./tau_r_i-s_i./tau_d_i;

	v_e_tmp=v_e(k)+dt05*v_e_inc;
	n_e_tmp=n_e+dt05*n_e_inc;
    h_e_tmp=h_e+dt05*h_e_inc;
    q_e_tmp=q_e+dt05*q_e_inc;   
    s_e_tmp=s_e+dt05*s_e_inc;    
	v_i_tmp=v_i(k)+dt05*v_i_inc;
	n_i_tmp=n_i+dt05*n_i_inc;
    h_i_tmp=h_i+dt05*h_i_inc;
    q_i_tmp=q_i+dt05*q_i_inc;   
    s_i_tmp=s_i+dt05*s_i_inc;    

	v_e_inc=0.1*(-67-v_e_tmp)+80*n_e_tmp.^4.*(-100-v_e_tmp) ...
           +100*m_e_inf(v_e_tmp).^3.*h_e_tmp.*(50-v_e_tmp) ...
           +(g_ie*s_i_tmp).*(v_rev_i-v_e_tmp) ...
           +i_ext_e;
	n_e_inc=(n_e_inf(v_e_tmp)-n_e_tmp)./tau_n_e(v_e_tmp);
    h_e_inc=(h_e_inf(v_e_tmp)-h_e_tmp)./tau_h_e(v_e_tmp);
    q_e_inc=(1+tanh(v_e_tmp/10))/2.*(1-q_e_tmp)/0.1-q_e_tmp./tau_dq_e;
    s_e_inc=q_e_tmp.*(1-s_e_tmp)./tau_r_e-s_e_tmp./tau_d_e;
	v_i_inc=0.1*(-65-v_i_tmp)+9*n_i_tmp.^4.*(-90-v_i_tmp) ...
           +35*m_i_inf(v_i_tmp).^3.*h_i_tmp.*(55-v_i_tmp) ...
           +(g_ei*s_e_tmp).*(v_rev_e-v_i_tmp) ...
           +i_ext_i;
	n_i_inc=(n_i_inf(v_i_tmp)-n_i_tmp)./tau_n_i(v_i_tmp);
    h_i_inc=(h_i_inf(v_i_tmp)-h_i_tmp)./tau_h_i(v_i_tmp);
    q_i_inc=(1+tanh(v_i_tmp/10))/2.*(1-q_i_tmp)/0.1-q_i_tmp./tau_dq_i;
    s_i_inc=q_i_tmp.*(1-s_i_tmp)./tau_r_i-s_i_tmp./tau_d_i;
        
	v_e(k+1)=v_e(k)+dt*v_e_inc;
	h_e=h_e+dt*h_e_inc; 
    n_e=n_e+dt*n_e_inc; 
    q_e=q_e+dt*q_e_inc;
    s_e=s_e+dt*s_e_inc;
	v_i(k+1)=v_i(k)+dt*v_i_inc;
    h_i=h_i+dt*h_i_inc; 
    n_i=n_i+dt*n_i_inc; 
    q_i=q_i+dt*q_i_inc;
    s_i=s_i+dt*s_i_inc;
    
    if v_e(k+1)<-20 & v_e(k)>=-20,
        num_spikes_e=num_spikes_e+1;
        t_e_spikes(num_spikes_e)=t_old*(-20-v_e(k+1))+t_new*(v_e(k)+20);
        t_e_spikes(num_spikes_e)=t_e_spikes(num_spikes_e)/(v_e(k)-v_e(k+1));
    end;
    
    s_i_plot(k)=s_i;

end;
period=t_e_spikes(num_spikes_e)-t_e_spikes(num_spikes_e-1)

% plot the voltage traces

subplot(211);
t=[0:m_steps]*dt;
plot(t,v_e,'-r','Linewidth',2);
hold on;
plot(t,v_i,'-b','Linewidth',2);
hold off;
set(gca,'Fontsize',16);
xlabel('$t$','Fontsize',20);
title('voltage traces of E-cell (red) and I-cell (blue)','Fontsize',20);
shg;
