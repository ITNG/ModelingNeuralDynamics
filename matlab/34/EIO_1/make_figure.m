clear; clf; rng('default'); rng(63806);
tic;

% Define strengths of h- and A-currents in O-LM cells

g_h=12;
g_A=22;

% Define network parameters: 

num_e=200; num_i=50; num_o=50;
sigma_e=0.05; i_ext_e=1.8*ones(num_e,1).*(1+sigma_e*randn(num_e,1)); 
sigma_i=0.10; i_ext_i=1.0*ones(num_i,1).*(1+sigma_i*randn(num_i,1));
sigma_o=0.05; i_ext_o=-2.0*ones(num_o,1).*(1+sigma_o*randn(num_o,1));

g_hat_ee=0.00; g_hat_ei=0.25; g_hat_eo=0.00; 
g_hat_ie=0.25; g_hat_ii=0.25; g_hat_io=0.50; 
g_hat_oe=1.00; g_hat_oi=0.50; g_hat_oo=0.00; 

p_ee=1.0; p_ei=1.0; p_eo=1.0;
p_ie=1.0; p_ii=1.0; p_io=1.0;
p_oe=1.0; p_oi=1.0; p_oo=1.0;

        % See explanation below to understand what the preceding 
        % parameters mean.
        
v_rev_e=0; v_rev_i=-75; v_rev_o=-75;
tau_rise_e=0.5; tau_peak_e=0.5; tau_d_e=3; 
tau_rise_i=0.5; tau_peak_i=0.5; tau_d_i=9; 
tau_rise_o=0.5; tau_peak_o=0.5; tau_d_o=20; 

t_final=1000;    % Time (in ms) simulated. 
dt=0.01;        % Time step used in solving the differential equations.


% Process network parameters a bit:

u_ee=rand(num_e,num_e); u_ei=rand(num_e,num_i); u_eo=rand(num_e,num_o);
u_ie=rand(num_i,num_e); u_ii=rand(num_i,num_i); u_io=rand(num_i,num_o);
u_oe=rand(num_o,num_e); u_oi=rand(num_o,num_i); u_oo=rand(num_o,num_o);

g_ee=g_hat_ee*(u_ee<p_ee)/(num_e*p_ee); 
g_ei=g_hat_ei*(u_ei<p_ei)/(num_e*p_ei);
g_eo=g_hat_eo*(u_eo<p_eo)/(num_e*p_eo);

g_ie=g_hat_ie*(u_ie<p_ie)/(num_i*p_ie); 
g_ii=g_hat_ii*(u_ii<p_ii)/(num_i*p_ii);
g_io=g_hat_io*(u_io<p_io)/(num_i*p_io);

g_oe=g_hat_oe*(u_oe<p_oe)/(num_o*p_oe); 
g_oi=g_hat_oi*(u_oi<p_oi)/(num_o*p_oi);
g_oo=g_hat_oo*(u_oo<p_oo)/(num_o*p_oo);


        % Consider, for example, the i-th e-cell and the j-th i-cell. the
        % probability that there is a synaptic connection at all from the
        % i-th e-cell to the j-th i-cell is p_ei. if there is such a
        % connection, its strength is g_hat_ei/(num_e*p_ei). Note that
        % num_e*p_ei is the expected number of excitatory inputs into an
        % inhibitory cell. Therefore dividing by this quantity has the
        % effect that the expected value of the total excitatory
        % conductance affecting an inhibitory cell is g_hat_ei. 

tau_dq_e=tau_d_q_function(tau_d_e,tau_rise_e,tau_peak_e);
tau_dq_i=tau_d_q_function(tau_d_i,tau_rise_i,tau_peak_i);
tau_dq_o=tau_d_q_function(tau_d_o,tau_rise_o,tau_peak_o);

dt05=dt/2; m_steps=round(t_final/dt);

% initialize dynamic variables

iv=rtm_init(i_ext_e,rand(num_e,1));
v_e=iv(:,1); m_e=m_e_inf(v_e); h_e=iv(:,2); n_e=iv(:,3); 
q_e=zeros(num_e,1); s_e=zeros(num_e,1);

iv=wb_init(i_ext_i,rand(num_i,1));
v_i=iv(:,1); m_i=m_i_inf(v_i); h_i=iv(:,2); n_i=iv(:,3); 
q_i=zeros(num_i,1); s_i=zeros(num_i,1);

iv=olm_init(i_ext_o,rand(num_o,1));
v_o=iv(:,1); m_o=m_o_inf(v_o); h_o=iv(:,2); n_o=iv(:,3); 
r_o=iv(:,4); a_o=iv(:,5); b_o=iv(:,6);
q_o=zeros(num_o,1); s_o=zeros(num_o,1);

% solve the system of Hodgkin-Huxley-like equations using the midpoint method

num_spikes_e=0; t_e_spikes=[]; i_e_spikes=[];
num_spikes_i=0; t_i_spikes=[]; i_i_spikes=[];
num_spikes_o=0; t_o_spikes=[]; i_o_spikes=[];

lfp_v(1)=mean(v_e);
lfp_s(1)=mean(s_e);

for k=1:m_steps,
    
    t_old=(k-1)*dt; t_new=k*dt; t_mid=(t_old+t_new)/2;
	
	v_e_inc=0.1*(-67-v_e)+80*n_e.^4.*(-100-v_e)+100*m_e.^3.*h_e.*(50-v_e) ...
               +(g_ee'*s_e).*(v_rev_e-v_e) ...
               +(g_ie'*s_i).*(v_rev_i-v_e) ...
               +(g_oe'*s_o).*(v_rev_o-v_e) ...
               +i_ext_e;
	n_e_inc=(n_e_inf(v_e)-n_e)./tau_n_e(v_e);
    h_e_inc=(h_e_inf(v_e)-h_e)./tau_h_e(v_e);
    q_e_inc=(1+tanh(v_e/10))/2.*(1-q_e)/0.1-q_e./tau_dq_e;
    s_e_inc=q_e.*(1-s_e)./tau_rise_e-s_e./tau_d_e;
    
	v_i_inc=0.1*(-65-v_i)+9*n_i.^4.*(-90-v_i)+35*m_i.^3.*h_i.*(55-v_i) ...
               +(g_ei'*s_e).*(v_rev_e-v_i) ...
               +(g_ii'*s_i).*(v_rev_i-v_i) ...
               +(g_oi'*s_o).*(v_rev_o-v_i) ...
               +i_ext_i;
	n_i_inc=(n_i_inf(v_i)-n_i)./tau_n_i(v_i);
    h_i_inc=(h_i_inf(v_i)-h_i)./tau_h_i(v_i);
    q_i_inc=(1+tanh(v_i/10))/2.*(1-q_i)/0.1-q_i./tau_dq_i;
    s_i_inc=q_i.*(1-s_i)./tau_rise_i-s_i./tau_d_i;
    
    v_o_inc=0.05*(-70-v_o)+23*n_o.^4.*(-100-v_o)+30*m_o.^3.*h_o.*(90-v_o) ...
               +g_h*r_o.*(-32.9-v_o) ...
               +g_A*a_o.*b_o.*(-90-v_o) ...
               +(g_eo'*s_e).*(v_rev_e-v_o) ...
               +(g_io'*s_i).*(v_rev_i-v_o) ...
               +(g_oo'*s_o).*(v_rev_o-v_o) ...
               +i_ext_o;
    v_o_inc=v_o_inc/1.3;
	n_o_inc=(n_o_inf(v_o)-n_o)./tau_n_o(v_o);
    h_o_inc=(h_o_inf(v_o)-h_o)./tau_h_o(v_o);
    r_o_inc=(r_o_inf(v_o)-r_o)./tau_r_o(v_o);
    a_o_inc=(a_o_inf(v_o)-a_o)./tau_a_o(v_o);
    b_o_inc=(b_o_inf(v_o)-b_o)./tau_b_o(v_o);
    q_o_inc=(1+tanh(v_o/10))/2.*(1-q_o)/0.1-q_o./tau_dq_o;
    s_o_inc=q_o.*(1-s_o)./tau_rise_o-s_o./tau_d_o;


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
    
    v_o_tmp=v_o+dt05*v_o_inc;
	n_o_tmp=n_o+dt05*n_o_inc;
	m_o_tmp=m_o_inf(v_o_tmp);
    h_o_tmp=h_o+dt05*h_o_inc;
    r_o_tmp=r_o+dt05*r_o_inc;
    a_o_tmp=a_o+dt05*a_o_inc;
    b_o_tmp=b_o+dt05*b_o_inc;
    q_o_tmp=q_o+dt05*q_o_inc;   
    s_o_tmp=s_o+dt05*s_o_inc;   
    
    

	v_e_inc=0.1*(-67-v_e_tmp)+80*n_e_tmp.^4.*(-100-v_e_tmp) ...
               +100*m_e_tmp.^3.*h_e_tmp.*(50-v_e_tmp) ...
               +(g_ee'*s_e_tmp).*(v_rev_e-v_e_tmp) ...
               +(g_ie'*s_i_tmp).*(v_rev_i-v_e_tmp) ...
               +(g_oe'*s_o_tmp).*(v_rev_o-v_e_tmp) ...
               +i_ext_e;
	n_e_inc=(n_e_inf(v_e_tmp)-n_e_tmp)./tau_n_e(v_e_tmp);
    h_e_inc=(h_e_inf(v_e_tmp)-h_e_tmp)./tau_h_e(v_e_tmp);
    q_e_inc=(1+tanh(v_e_tmp/10))/2.*(1-q_e_tmp)/0.1-q_e_tmp./tau_dq_e;
    s_e_inc=q_e_tmp.*(1-s_e_tmp)./tau_rise_e-s_e_tmp./tau_d_e;
    
	v_i_inc=0.1*(-65-v_i_tmp)+9*n_i_tmp.^4.*(-90-v_i_tmp) ...
               +35*m_i_tmp.^3.*h_i_tmp.*(55-v_i_tmp) ...
               +(g_ei'*s_e_tmp).*(v_rev_e-v_i_tmp) ...
               +(g_ii'*s_i_tmp).*(v_rev_i-v_i_tmp) ...
               +(g_oi'*s_o_tmp).*(v_rev_o-v_i_tmp) ...
               +i_ext_i;
	n_i_inc=(n_i_inf(v_i_tmp)-n_i_tmp)./tau_n_i(v_i_tmp);
    h_i_inc=(h_i_inf(v_i_tmp)-h_i_tmp)./tau_h_i(v_i_tmp);
    q_i_inc=(1+tanh(v_i_tmp/10))/2.*(1-q_i_tmp)/0.1-q_i_tmp./tau_dq_i;
    s_i_inc=q_i_tmp.*(1-s_i_tmp)./tau_rise_i-s_i_tmp./tau_d_i;
    
    v_o_inc=0.05*(-70-v_o_tmp)+23*n_o_tmp.^4.*(-100-v_o_tmp) ...
               +30*m_o_tmp.^3.*h_o_tmp.*(90-v_o_tmp) ...
               +g_h*r_o_tmp.*(-32.9-v_o_tmp) ...
               +g_A*a_o_tmp.*b_o_tmp.*(-90-v_o_tmp) ...
               +(g_eo'*s_e_tmp).*(v_rev_e-v_o_tmp) ...
               +(g_io'*s_i_tmp).*(v_rev_i-v_o_tmp) ...
               +(g_oo'*s_o_tmp).*(v_rev_o-v_o_tmp) ...
               +i_ext_o;
    v_o_inc=v_o_inc/1.3;
	n_o_inc=(n_o_inf(v_o_tmp)-n_o_tmp)./tau_n_o(v_o_tmp);
    h_o_inc=(h_o_inf(v_o_tmp)-h_o_tmp)./tau_h_o(v_o_tmp);
    r_o_inc=(r_o_inf(v_o_tmp)-r_o_tmp)./tau_r_o(v_o_tmp);
    a_o_inc=(a_o_inf(v_o_tmp)-a_o_tmp)./tau_a_o(v_o_tmp);
    b_o_inc=(b_o_inf(v_o_tmp)-b_o_tmp)./tau_b_o(v_o_tmp);
    q_o_inc=(1+tanh(v_o_tmp/10))/2.*(1-q_o_tmp)/0.1-q_o_tmp./tau_dq_o;
    s_o_inc=q_o_tmp.*(1-s_o_tmp)./tau_rise_o-s_o_tmp./tau_d_o;
        
    v_e_old=v_e;
    v_i_old=v_i;
    v_o_old=v_o;

	v_e=v_e+dt*v_e_inc;
	m_e=m_e_inf(v_e); 
    h_e=h_e+dt*h_e_inc; 
    n_e=n_e+dt*n_e_inc; 
    q_e=q_e+dt*q_e_inc;
    s_e=s_e+dt*s_e_inc;
    
	v_i=v_i+dt*v_i_inc;
    m_i=m_i_inf(v_i); 
    h_i=h_i+dt*h_i_inc; 
    n_i=n_i+dt*n_i_inc; 
    q_i=q_i+dt*q_i_inc;
    s_i=s_i+dt*s_i_inc;
    
    v_o=v_o+dt*v_o_inc;
	n_o=n_o+dt*n_o_inc;
	m_o=m_o_inf(v_o);
    h_o=h_o+dt*h_o_inc;
    r_o=r_o+dt*r_o_inc;
    a_o=a_o+dt*a_o_inc;
    b_o=b_o+dt*b_o_inc;
    q_o=q_o+dt*q_o_inc;   
    s_o=s_o+dt*s_o_inc;   
    
    

	% Determine which and how many e-, i-, and o-cells spiked in the current 
    % time step:

    which_e=find(v_e_old>-20 & v_e <=-20); 
    which_i=find(v_i_old>-20 & v_i <=-20);
    which_o=find(v_o_old>-20 & v_o <=-20);
    
    l_e=length(which_e); 
    l_i=length(which_i);
    l_o=length(which_o);
    
    if l_e>0, 
        range=num_spikes_e+1:num_spikes_e+l_e; i_e_spikes(range)=which_e; 
        t_e_spikes(range)= ...
            ((-20-v_e(which_e))*(k-1)*dt+(v_e_old(which_e)+20)*k*dt)./ ...
                (-v_e(which_e)+v_e_old(which_e));
        num_spikes_e=num_spikes_e+l_e;
    end 
    if l_i>0, 
        range=num_spikes_i+1:num_spikes_i+l_i; i_i_spikes(range)=which_i; 
        t_i_spikes(range)= ...
            ((-20-v_i(which_i))*(k-1)*dt+(v_i_old(which_i)+20)*k*dt)./ ...
                (-v_i(which_i)+v_i_old(which_i));
        num_spikes_i=num_spikes_i+l_i;
    end 
    if l_o>0, 
        range=num_spikes_o+1:num_spikes_o+l_o; i_o_spikes(range)=which_o; 
        t_o_spikes(range)= ...
            ((-20-v_o(which_o))*(k-1)*dt+(v_o_old(which_o)+20)*k*dt)./ ...
                (-v_o(which_o)+v_o_old(which_o));
        num_spikes_o=num_spikes_o+l_o;
    end 
    
    lfp_v(k+1)=mean(v_e);
    lfp_s(k+1)=mean(s_e);
    
end;


% plot the spike rastergram

rastergram;
toc

