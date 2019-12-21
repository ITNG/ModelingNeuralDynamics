clear; clf; 
tic;

% Define network parameters: 

i_ext_i=1.5;
g_ii=0.5;        
v_rev_i=-75; 
tau_r_i=0.5; tau_peak_i=0.5; tau_d_i=9; 
t_final=200;   % Time (in ms) simulated. 
dt=0.01;       % Time step used in solving the differential equations.


% Process network parameters:

tau_dq_i=tau_d_q_function(tau_d_i,tau_r_i,tau_peak_i);
dt05=dt/2; m_steps=round(t_final/dt);



% initialize dynamic variables

v_i(1)=-75; h_i=0.1; n_i=0.1; q_i=0; s_i=0;

% solve the system of Hodgkin-Huxley-like equations using the midpoint method

num_spikes_i=0; 

for k=1:m_steps,
    t_old=(k-1)*dt; t_new=k*dt; t_mid=(t_old+t_new)/2;
	
	v_i_inc=0.1*(-65-v_i(k))+9*n_i.^4.*(-90-v_i(k)) ...
           +35*m_i_inf(v_i(k)).^3.*h_i.*(55-v_i(k)) ...
           +(g_ii*s_i).*(v_rev_i-v_i(k)) ...
           +i_ext_i;
	n_i_inc=(n_i_inf(v_i(k))-n_i)./tau_n_i(v_i(k));
    h_i_inc=(h_i_inf(v_i(k))-h_i)./tau_h_i(v_i(k));
    q_i_inc=(1+tanh(v_i(k)/10))/2.*(1-q_i)/0.1-q_i./tau_dq_i;
    s_i_inc=q_i.*(1-s_i)./tau_r_i-s_i./tau_d_i;

	 
	v_i_tmp=v_i(k)+dt05*v_i_inc;
	n_i_tmp=n_i+dt05*n_i_inc;
    h_i_tmp=h_i+dt05*h_i_inc;
    q_i_tmp=q_i+dt05*q_i_inc;   
    s_i_tmp=s_i+dt05*s_i_inc;    

	
	v_i_inc=0.1*(-65-v_i_tmp)+9*n_i_tmp.^4.*(-90-v_i_tmp) ...
           +35*m_i_inf(v_i_tmp).^3.*h_i_tmp.*(55-v_i_tmp) ...
           +(g_ii'*s_i_tmp).*(v_rev_i-v_i_tmp) ...
           +i_ext_i;
	n_i_inc=(n_i_inf(v_i_tmp)-n_i_tmp)./tau_n_i(v_i_tmp);
    h_i_inc=(h_i_inf(v_i_tmp)-h_i_tmp)./tau_h_i(v_i_tmp);
    q_i_inc=(1+tanh(v_i_tmp/10))/2.*(1-q_i_tmp)/0.1-q_i_tmp./tau_dq_i;
    s_i_inc=q_i_tmp.*(1-s_i_tmp)./tau_r_i-s_i_tmp./tau_d_i;
        

	v_i(k+1)=v_i(k)+dt*v_i_inc; 
    h_i=h_i+dt*h_i_inc; 
    n_i=n_i+dt*n_i_inc; 
    q_i=q_i+dt*q_i_inc;
    s_i=s_i+dt*s_i_inc;
    
    if v_i(k+1)<-20 & v_i(k)>=-20,
        num_spikes_i=num_spikes_i+1;
        t_i_spikes(num_spikes_i)=t_old*(-20-v_i(k+1))+t_new*(v_i(k)+20);
        t_i_spikes(num_spikes_i)=t_i_spikes(num_spikes_i)/(v_i(k)-v_i(k+1));
    end;

end;

base_period=t_i_spikes(num_spikes_i)-t_i_spikes(num_spikes_i-1);

% plot the voltage trace

subplot(211);
t=[0:m_steps]*dt;
plot(t,v_i,'-b','Linewidth',2);
set(gca,'Fontsize',16); 
xlabel('$t$ [ms]','Fontsize',20);
ylabel('$v$ [mV]','Fontsize',20);

shg;

i_ext_i=1.5; g_ii=0.5; tau_d_i=9;
i_ext_i=i_ext_i*0.99;


% initialize dynamic variables

v_i(1)=-75; h_i=0.1; n_i=0.1; q_i=0; s_i=0;

% solve the system of Hodgkin-Huxley-like equations using the midpoint method

num_spikes_i=0; 

for k=1:m_steps,
    t_old=(k-1)*dt; t_new=k*dt; t_mid=(t_old+t_new)/2;
	
	v_i_inc=0.1*(-65-v_i(k))+9*n_i.^4.*(-90-v_i(k)) ...
           +35*m_i_inf(v_i(k)).^3.*h_i.*(55-v_i(k)) ...
           +(g_ii*s_i).*(v_rev_i-v_i(k)) ...
           +i_ext_i;
	n_i_inc=(n_i_inf(v_i(k))-n_i)./tau_n_i(v_i(k));
    h_i_inc=(h_i_inf(v_i(k))-h_i)./tau_h_i(v_i(k));
    q_i_inc=(1+tanh(v_i(k)/10))/2.*(1-q_i)/0.1-q_i./tau_dq_i;
    s_i_inc=q_i.*(1-s_i)./tau_r_i-s_i./tau_d_i;

	 
	v_i_tmp=v_i(k)+dt05*v_i_inc;
	n_i_tmp=n_i+dt05*n_i_inc;
    h_i_tmp=h_i+dt05*h_i_inc;
    q_i_tmp=q_i+dt05*q_i_inc;   
    s_i_tmp=s_i+dt05*s_i_inc;    

	
	v_i_inc=0.1*(-65-v_i_tmp)+9*n_i_tmp.^4.*(-90-v_i_tmp) ...
           +35*m_i_inf(v_i_tmp).^3.*h_i_tmp.*(55-v_i_tmp) ...
           +(g_ii'*s_i_tmp).*(v_rev_i-v_i_tmp) ...
           +i_ext_i;
	n_i_inc=(n_i_inf(v_i_tmp)-n_i_tmp)./tau_n_i(v_i_tmp);
    h_i_inc=(h_i_inf(v_i_tmp)-h_i_tmp)./tau_h_i(v_i_tmp);
    q_i_inc=(1+tanh(v_i_tmp/10))/2.*(1-q_i_tmp)/0.1-q_i_tmp./tau_dq_i;
    s_i_inc=q_i_tmp.*(1-s_i_tmp)./tau_r_i-s_i_tmp./tau_d_i;
        

	v_i(k+1)=v_i(k)+dt*v_i_inc; 
    h_i=h_i+dt*h_i_inc; 
    n_i=n_i+dt*n_i_inc; 
    q_i=q_i+dt*q_i_inc;
    s_i=s_i+dt*s_i_inc;
    
    if v_i(k+1)<-20 & v_i(k)>=-20,
        num_spikes_i=num_spikes_i+1;
        t_i_spikes(num_spikes_i)=t_old*(-20-v_i(k+1))+t_new*(v_i(k)+20);
        t_i_spikes(num_spikes_i)=t_i_spikes(num_spikes_i)/(v_i(k)-v_i(k+1));
    end;

end;

period_reduced_I=t_i_spikes(num_spikes_i)-t_i_spikes(num_spikes_i-1);
percentage_change=(base_period-period_reduced_I)/base_period*100

i_ext_i=1.5; g_ii=0.5; tau_d_i=9;
g_ii=g_ii*1.01;

% initialize dynamic variables

v_i(1)=-75; h_i=0.1; n_i=0.1; q_i=0; s_i=0;

% solve the system of Hodgkin-Huxley-like equations using the midpoint method

num_spikes_i=0; 

for k=1:m_steps,
    t_old=(k-1)*dt; t_new=k*dt; t_mid=(t_old+t_new)/2;
	
	v_i_inc=0.1*(-65-v_i(k))+9*n_i.^4.*(-90-v_i(k)) ...
           +35*m_i_inf(v_i(k)).^3.*h_i.*(55-v_i(k)) ...
           +(g_ii*s_i).*(v_rev_i-v_i(k)) ...
           +i_ext_i;
	n_i_inc=(n_i_inf(v_i(k))-n_i)./tau_n_i(v_i(k));
    h_i_inc=(h_i_inf(v_i(k))-h_i)./tau_h_i(v_i(k));
    q_i_inc=(1+tanh(v_i(k)/10))/2.*(1-q_i)/0.1-q_i./tau_dq_i;
    s_i_inc=q_i.*(1-s_i)./tau_r_i-s_i./tau_d_i;

	 
	v_i_tmp=v_i(k)+dt05*v_i_inc;
	n_i_tmp=n_i+dt05*n_i_inc;
    h_i_tmp=h_i+dt05*h_i_inc;
    q_i_tmp=q_i+dt05*q_i_inc;   
    s_i_tmp=s_i+dt05*s_i_inc;    

	
	v_i_inc=0.1*(-65-v_i_tmp)+9*n_i_tmp.^4.*(-90-v_i_tmp) ...
           +35*m_i_inf(v_i_tmp).^3.*h_i_tmp.*(55-v_i_tmp) ...
           +(g_ii'*s_i_tmp).*(v_rev_i-v_i_tmp) ...
           +i_ext_i;
	n_i_inc=(n_i_inf(v_i_tmp)-n_i_tmp)./tau_n_i(v_i_tmp);
    h_i_inc=(h_i_inf(v_i_tmp)-h_i_tmp)./tau_h_i(v_i_tmp);
    q_i_inc=(1+tanh(v_i_tmp/10))/2.*(1-q_i_tmp)/0.1-q_i_tmp./tau_dq_i;
    s_i_inc=q_i_tmp.*(1-s_i_tmp)./tau_r_i-s_i_tmp./tau_d_i;
        

	v_i(k+1)=v_i(k)+dt*v_i_inc; 
    h_i=h_i+dt*h_i_inc; 
    n_i=n_i+dt*n_i_inc; 
    q_i=q_i+dt*q_i_inc;
    s_i=s_i+dt*s_i_inc;
    
    if v_i(k+1)<-20 & v_i(k)>=-20,
        num_spikes_i=num_spikes_i+1;
        t_i_spikes(num_spikes_i)=t_old*(-20-v_i(k+1))+t_new*(v_i(k)+20);
        t_i_spikes(num_spikes_i)=t_i_spikes(num_spikes_i)/(v_i(k)-v_i(k+1));
    end;

end;

period_raised_g_ii=t_i_spikes(num_spikes_i)-t_i_spikes(num_spikes_i-1);
percentage_change=(base_period-period_raised_g_ii)/base_period*100

i_ext_i=1.5; g_ii=0.5; tau_d_i=9;
tau_d_i=tau_d_i*1.01;

% Process network parameters:

tau_dq_i=tau_d_q_function(tau_d_i,tau_r_i,tau_peak_i);

% initialize dynamic variables

v_i(1)=-75; h_i=0.1; n_i=0.1; q_i=0; s_i=0;

% solve the system of Hodgkin-Huxley-like equations using the midpoint method

num_spikes_i=0; 

for k=1:m_steps,
    t_old=(k-1)*dt; t_new=k*dt; t_mid=(t_old+t_new)/2;
	
	v_i_inc=0.1*(-65-v_i(k))+9*n_i.^4.*(-90-v_i(k)) ...
           +35*m_i_inf(v_i(k)).^3.*h_i.*(55-v_i(k)) ...
           +(g_ii*s_i).*(v_rev_i-v_i(k)) ...
           +i_ext_i;
	n_i_inc=(n_i_inf(v_i(k))-n_i)./tau_n_i(v_i(k));
    h_i_inc=(h_i_inf(v_i(k))-h_i)./tau_h_i(v_i(k));
    q_i_inc=(1+tanh(v_i(k)/10))/2.*(1-q_i)/0.1-q_i./tau_dq_i;
    s_i_inc=q_i.*(1-s_i)./tau_r_i-s_i./tau_d_i;

	 
	v_i_tmp=v_i(k)+dt05*v_i_inc;
	n_i_tmp=n_i+dt05*n_i_inc;
    h_i_tmp=h_i+dt05*h_i_inc;
    q_i_tmp=q_i+dt05*q_i_inc;   
    s_i_tmp=s_i+dt05*s_i_inc;    

	
	v_i_inc=0.1*(-65-v_i_tmp)+9*n_i_tmp.^4.*(-90-v_i_tmp) ...
           +35*m_i_inf(v_i_tmp).^3.*h_i_tmp.*(55-v_i_tmp) ...
           +(g_ii'*s_i_tmp).*(v_rev_i-v_i_tmp) ...
           +i_ext_i;
	n_i_inc=(n_i_inf(v_i_tmp)-n_i_tmp)./tau_n_i(v_i_tmp);
    h_i_inc=(h_i_inf(v_i_tmp)-h_i_tmp)./tau_h_i(v_i_tmp);
    q_i_inc=(1+tanh(v_i_tmp/10))/2.*(1-q_i_tmp)/0.1-q_i_tmp./tau_dq_i;
    s_i_inc=q_i_tmp.*(1-s_i_tmp)./tau_r_i-s_i_tmp./tau_d_i;
        

	v_i(k+1)=v_i(k)+dt*v_i_inc; 
    h_i=h_i+dt*h_i_inc; 
    n_i=n_i+dt*n_i_inc; 
    q_i=q_i+dt*q_i_inc;
    s_i=s_i+dt*s_i_inc;
    
    if v_i(k+1)<-20 & v_i(k)>=-20,
        num_spikes_i=num_spikes_i+1;
        t_i_spikes(num_spikes_i)=t_old*(-20-v_i(k+1))+t_new*(v_i(k)+20);
        t_i_spikes(num_spikes_i)=t_i_spikes(num_spikes_i)/(v_i(k)-v_i(k+1));
    end;

end;

period_raised_tau_d_i=t_i_spikes(num_spikes_i)-t_i_spikes(num_spikes_i-1);
percentage_change=(base_period-period_raised_tau_d_i)/base_period*100



toc

