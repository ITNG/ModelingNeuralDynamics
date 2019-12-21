clear; clf; rng('default'); rng(63806);
tic;

% Define network parameters: 

num_i=100;
sigma_i=0.05;  i_ext_i=1.5*(1+randn(num_i,1)*sigma_i);
g_hat_ii=0.5;   p_ii=0.5;
g_hat_gap=0.1;  p_gap=0.05;    % Individual gap junctions are of strength
                               % g_hat_gap/(p_gap*(num_i-1)).
                               % Any two different I-cells are gap-junctionally
                               % connected with probability p_gap.

v_rev_i=-75; 
tau_r_i=0.5; tau_peak_i=0.5; tau_d_i=9; 
t_final=500;   % Time (in ms) simulated. 
dt=0.01;       % Time step used in solving the differential equations.

% Process network parameters:

u_ii=rand(num_i,num_i); 
g_ii=g_hat_ii*(u_ii<p_ii)/(num_i*p_ii);
        
tau_dq_i=tau_d_q_function(tau_d_i,tau_r_i,tau_peak_i);

dt05=dt/2; m_steps=round(t_final/dt);

for i=1:num_i-1,
    for j=i+1:num_i,
        u=rand(1,1);
        G_gap(i,j)=(u<p_gap)*g_hat_gap/(p_gap*(num_i-1));
        G_gap(j,i)=G_gap(i,j);
    end;
end;
c=sum(G_gap); c=c';


% initialize dynamic variables

iv=wb_init(i_ext_i,rand(num_i,1));
v_i=iv(:,1); m_i=m_i_inf(v_i); h_i=iv(:,2); n_i=iv(:,3); 
q_i=zeros(num_i,1); s_i=zeros(num_i,1);

% solve the system of Hodgkin-Huxley-like equations using the midpoint method

num_spikes_i=0; t_i_spikes=[]; i_i_spikes=[];

for k=1:m_steps,
    t_old=(k-1)*dt; t_new=k*dt; t_mid=(t_old+t_new)/2;
	
	v_i_inc=0.1*(-65-v_i)+9*n_i.^4.*(-90-v_i)+35*m_i.^3.*h_i.*(55-v_i) ...
               +(g_ii'*s_i).*(v_rev_i-v_i) ...
               +i_ext_i+G_gap*v_i-c.*v_i;
	n_i_inc=(n_i_inf(v_i)-n_i)./tau_n_i(v_i);
    h_i_inc=(h_i_inf(v_i)-h_i)./tau_h_i(v_i);
    q_i_inc=(1+tanh(v_i/10))/2.*(1-q_i)/0.1-q_i./tau_dq_i;
    s_i_inc=q_i.*(1-s_i)./tau_r_i-s_i./tau_d_i;

	 
	v_i_tmp=v_i+dt05*v_i_inc;
	n_i_tmp=n_i+dt05*n_i_inc;
	m_i_tmp=m_i_inf(v_i_tmp);
    h_i_tmp=h_i+dt05*h_i_inc;
    q_i_tmp=q_i+dt05*q_i_inc;   
    s_i_tmp=s_i+dt05*s_i_inc;    

	
	v_i_inc=0.1*(-65-v_i_tmp)+9*n_i_tmp.^4.*(-90-v_i_tmp)+ ...
               +35*m_i_tmp.^3.*h_i_tmp.*(55-v_i_tmp) ...
               +(g_ii'*s_i_tmp).*(v_rev_i-v_i_tmp) ...
               +i_ext_i+G_gap*v_i_tmp-c.*v_i_tmp;
	n_i_inc=(n_i_inf(v_i_tmp)-n_i_tmp)./tau_n_i(v_i_tmp);
    h_i_inc=(h_i_inf(v_i_tmp)-h_i_tmp)./tau_h_i(v_i_tmp);
    q_i_inc=(1+tanh(v_i_tmp/10))/2.*(1-q_i_tmp)/0.1-q_i_tmp./tau_dq_i;
    s_i_inc=q_i_tmp.*(1-s_i_tmp)./tau_r_i-s_i_tmp./tau_d_i;
        
    v_i_old=v_i;

	v_i=v_i+dt*v_i_inc;
    m_i=m_i_inf(v_i); 
    h_i=h_i+dt*h_i_inc; 
    n_i=n_i+dt*n_i_inc; 
    q_i=q_i+dt*q_i_inc;
    s_i=s_i+dt*s_i_inc;
    

	% Determine which and how many i-cells spiked in the current 
    % time step:

    which_i=find(v_i_old>-20 & v_i <=-20);
    l_i=length(which_i);
    if l_i>0, 
        range=num_spikes_i+1:num_spikes_i+l_i; i_i_spikes(range)=which_i; 
        t_i_spikes(range)= ...
            ((-20-v_i(which_i))*(k-1)*dt+(v_i_old(which_i)+20)*k*dt)./ ...
                (-v_i(which_i)+v_i_old(which_i));
        num_spikes_i=num_spikes_i+l_i;
    end 
    

end;

% plot the spike rastergram

subplot(211);
rastergram;
set(gca,'Fontsize',16); 
set(gca,'Ytick',[1,num_i]);
xlabel('$t$ [ms]','Fontsize',20)
axis([0,t_final,0,num_i+1]); 

shg;

toc

