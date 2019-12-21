clear; clf; rng('default'); rng(63806);
tic;

% Define network parameters: 


num_e=200; 
sigma_e=0.05; i_ext_e=3.0*ones(num_e,1).*(1+sigma_e*randn(num_e,1)); 
g_m=1.0; % M-current conductance; the only parameter of the neuronal model,
         % other than external drive, not explicitly written into the code.
g_hat_gap=0.2;
p_gap=0.1;
        

t_final=500;   % Time (in ms) simulated. 
dt=0.01;       % Time step used in solving the differential equations.


% Process network parameters a bit:

g_gap=zeros(num_e,num_e);
for i=1:num_e-1,
    for j=i+1:num_e,
        u=rand(1,1);
        if u<p_gap,
            g_gap(i,j)=g_hat_gap/(num_e-1)/p_gap;
            g_gap(j,i)=g_gap(i,j);
        end;
    end;
end;
g_gap_column_sum=sum(g_gap)';
     

dt05=dt/2; m_steps=round(t_final/dt);

% initialize dynamic variables

iv=rtm_init_with_m_current(i_ext_e,rand(num_e,1),g_m);
v_e=iv(:,1); m_e=m_e_inf(v_e); h_e=iv(:,2); n_e=iv(:,3); w=iv(:,4);

% solve the system of Hodgkin-Huxley-like equations using the midpoint method

num_spikes_e=0; t_e_spikes=[]; i_e_spikes=[];

for k=1:m_steps,
    t_old=(k-1)*dt; t_new=k*dt; t_mid=(t_old+t_new)/2;
	
	v_e_inc=0.1*(-67-v_e)+80*n_e.^4.*(-100-v_e)+100*m_e.^3.*h_e.*(50-v_e) ...
               +g_gap*v_e-g_gap_column_sum.*v_e ...
               +g_m*w.*(-100-v_e)+i_ext_e;
	n_e_inc=(n_e_inf(v_e)-n_e)./tau_n_e(v_e);
    h_e_inc=(h_e_inf(v_e)-h_e)./tau_h_e(v_e);
    w_inc=(w_inf(v_e)-w)./tau_w(v_e);
	
	v_e_tmp=v_e+dt05*v_e_inc;
	n_e_tmp=n_e+dt05*n_e_inc;
	m_e_tmp=m_e_inf(v_e_tmp);
    h_e_tmp=h_e+dt05*h_e_inc;
    w_tmp=w+dt05*w_inc;

	v_e_inc=0.1*(-67-v_e_tmp)+80*n_e_tmp.^4.*(-100-v_e_tmp)+100*m_e_tmp.^3.*h_e_tmp.*(50-v_e_tmp) ...
               +g_gap*v_e_tmp-g_gap_column_sum.*v_e_tmp ...
               +g_m*w_tmp.*(-100-v_e_tmp)+i_ext_e;
	n_e_inc=(n_e_inf(v_e_tmp)-n_e_tmp)./tau_n_e(v_e_tmp);
    h_e_inc=(h_e_inf(v_e_tmp)-h_e_tmp)./tau_h_e(v_e_tmp);
        
    v_e_old=v_e;

	v_e=v_e+dt*v_e_inc;
	m_e=m_e_inf(v_e); h_e=h_e+dt*h_e_inc; n_e=n_e+dt*n_e_inc; 
    w=w+dt*w_inc;
    

	% Determine which and how many e- and i-cells spiked in the current 
    % time step:

    which_e=find(v_e_old>-20 & v_e <=-20); 
    l_e=length(which_e); 
    if l_e>0, 
        range=num_spikes_e+1:num_spikes_e+l_e; i_e_spikes(range)=which_e; 
        t_e_spikes(range)= ...
            ((-20-v_e(which_e))*(k-1)*dt+(v_e_old(which_e)+20)*k*dt)./ ...
                (-v_e(which_e)+v_e_old(which_e));
        num_spikes_e=num_spikes_e+l_e;
    end 

end;

f_hat_e=round(num_spikes_e/num_e/t_final*1000)

% plot the spike rastergram

rastergram;
toc

