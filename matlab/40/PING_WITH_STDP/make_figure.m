clear; clf; rng('default'); rng(63806);
tic;

% Define network parameters: 

num_e=200; num_i=50;
i_ext_e=0.2+rand(num_e,1)*1.8;
i_ext_i=0.25*ones(num_i,1);
g_hat_ei=0.5; g_hat_ie=0.5; g_hat_ii=0.25; 
p_ei=0.5; p_ie=0.5; p_ii=0.5;
        % See explanation below to understand what the preceding eight
        % parameters mean.
v_rev_e=0; v_rev_i=-75; 
tau_r_e=0.5; tau_peak_e=0.5; tau_d_e=3; 
tau_r_i=0.5; tau_peak_i=0.5; tau_d_i=9; 
t_final=500;   % Time (in ms) simulated. 
dt=0.01;       % Time step used in solving the differential equations.


% Process network parameters a bit:

u_ei=rand(num_e,num_i);
u_ie=rand(num_i,num_e); u_ii=rand(num_i,num_i);
g_ei=g_hat_ei*(u_ei<p_ei)/(num_e*p_ei);
g_ie=g_hat_ie*(u_ie<p_ie)/(num_i*p_ie); g_ii=g_hat_ii*(u_ii<p_ii)/(num_i*p_ii);
        % Consider, for example, the i-th e-cell and the j-th i-cell. the
        % probability that there is a synaptic connection at all from the
        % i-th e-cell to the j-th i-cell is p_ei. if there is such a
        % connection, its strength is g_hat_ei/(num_e*p_ei). Note that
        % num_e*p_ei is the expected number of excitatory inputs into an
        % inhibitory cell. Therefore dividing by this quantity has the
        % effect that the expected value of the total excitatory
        % conductance affecting an inhibitory cell is g_hat_ei. 

tau_dq_e=tau_d_q_function(tau_d_e,tau_r_e,tau_peak_e);
tau_dq_i=tau_d_q_function(tau_d_i,tau_r_i,tau_peak_i);
dt05=dt/2; m_steps=round(t_final/dt);

% STDP parameters

g_ee=0.05*ones(num_e,num_e)/(num_e-1);   % initial strengths of 
for i=1:num_e,                          % recurrent excitatatory
    g_ee(i,i)=0;                        % connections.
end;

C=1.45;     % constant that appears in the approximate delta function
            % C*(1+tanh(v_e/10)) (with v_e=membrane potential of an E-cell)
K_plus=g_ee;         % maximum increase in E-to-E synaptic strength 
                     % resulting from a single spike pair
K_minus=g_ee*2/3;    % maximum decrease in E-to-E synaptic strength
                     % resulting from a single spike pair
tau_plus=10;         % If a post-synaptic spike follows a pre-synaptic
                     % one with with a delay delta, then the synaptic
                     % strength increases by
                     % K_plus*exp(-delta/tau_plus), for an E-to-E
                     % synapse.
tau_minus=10;        % analogous quantity for decreases in synaptic 
                     % strength resulting from a pre-synaptic spike
                     % following a post-synaptic one
B=8*g_ee;            % Upper bounds on the strengths of E-to-E synapses. 
                     % g_ee is time-dependent here, but B is fixed.

delta=max(g_ee/2,10^(-6));   % The changes resulting from STDP are constrained
                            % so the strength does not fall below zero, and
                            % does not rise above B. The function that 
                            % accomplishes this is R_delta, defined below.
                            % R_delta is almost piecewise linear, but the
                            % corners are smoothed, with transitions of
                            % width delta.
                            
% initialize dynamic variables

iv=rtm_init(i_ext_e,rand(num_e,1));
v_e=iv(:,1); m_e=m_e_inf(v_e); h_e=iv(:,2); n_e=iv(:,3); 
q_e=zeros(num_e,1); s_e=zeros(num_e,1);
a=100*ones(num_e,1);

v_i=-75*ones(num_i,1); m_i=m_i_inf(v_i); h_i=h_i_inf(v_i); n_i=n_i_inf(v_i); 
z=zeros(num_i,1); q_i=z; s_i=z; 

% solve the system of Hodgkin-Huxley-like equations using the midpoint method

num_spikes_e=0; t_e_spikes=[]; i_e_spikes=[];
num_spikes_i=0; t_i_spikes=[]; i_i_spikes=[];
lfp=zeros(m_steps+1,1); lfp(1)=mean(v_e);
for k=1:m_steps,
    t_old=(k-1)*dt; t_new=k*dt; t_mid=(t_old+t_new)/2;
	
	v_e_inc=0.1*(-67-v_e)+80*n_e.^4.*(-100-v_e)+100*m_e.^3.*h_e.*(50-v_e) ...
               +(g_ee'*s_e).*(v_rev_e-v_e)+(g_ie'*s_i).*(v_rev_i-v_e) ...
               +i_ext_e;
	n_e_inc=(n_e_inf(v_e)-n_e)./tau_n_e(v_e);
    h_e_inc=(h_e_inf(v_e)-h_e)./tau_h_e(v_e);
    q_e_inc=(1+tanh(v_e/10))/2.*(1-q_e)/0.1-q_e./tau_dq_e;
    s_e_inc=q_e.*(1-s_e)./tau_r_e-s_e./tau_d_e;
    a_inc=1-20*a.*(1+tanh(v_e/10));
    y=g_ee+diag(exp(-a/tau_plus).*(1-exp(-5*a/tau_plus)))*K_plus;
    y=min(B,y)-delta/2.*log(1+exp(-2*abs(B-y)./delta));
    y=y-g_ee;
    z=g_ee-K_minus*diag(exp(-a/tau_minus).*(1-exp(-5*a/tau_minus)));
    z=max(0,z)+delta/2.*log(1+exp(-2*abs(z)./delta));
    z=z-g_ee;
    g_ee_inc=C*(1+ones(num_e,1)*tanh(v_e'/10)).*y ...
            +C*(1+tanh(v_e/10)*ones(1,num_e)).*z;   
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
    a_tmp=a+dt05*a_inc;
    g_ee_tmp=g_ee+dt05*g_ee_inc;
	v_i_tmp=v_i+dt05*v_i_inc;
	n_i_tmp=n_i+dt05*n_i_inc;
	m_i_tmp=m_i_inf(v_i_tmp);
    h_i_tmp=h_i+dt05*h_i_inc;
    q_i_tmp=q_i+dt05*q_i_inc;   
    s_i_tmp=s_i+dt05*s_i_inc;    

	v_e_inc=0.1*(-67-v_e_tmp)+80*n_e_tmp.^4.*(-100-v_e_tmp)+100*m_e_tmp.^3.*h_e_tmp.*(50-v_e_tmp) ...
               +(g_ee'*s_e_tmp).*(v_rev_e-v_e_tmp)+(g_ie'*s_i_tmp).*(v_rev_i-v_e_tmp) ...
               +i_ext_e;
	n_e_inc=(n_e_inf(v_e_tmp)-n_e_tmp)./tau_n_e(v_e_tmp);
    h_e_inc=(h_e_inf(v_e_tmp)-h_e_tmp)./tau_h_e(v_e_tmp);
    q_e_inc=(1+tanh(v_e_tmp/10))/2.*(1-q_e_tmp)/0.1-q_e_tmp./tau_dq_e;
    s_e_inc=q_e_tmp.*(1-s_e_tmp)./tau_r_e-s_e_tmp./tau_d_e;
    a_inc=1-20*a_tmp.*(1+tanh(v_e_tmp/10));
    y=g_ee_tmp+ ...
        diag(exp(-a_tmp/tau_plus).*(1-exp(-5*a_tmp/tau_plus)))*K_plus;
    y=min(B,y)-delta/2.*log(1+exp(-2*abs(B-y)./delta));
    y=y-g_ee_tmp;
    z=g_ee_tmp- ...
        K_minus*diag(exp(-a_tmp/tau_minus).*(1-exp(-5*a_tmp/tau_minus)));
    z=max(0,z)+delta/2.*log(1+exp(-2*abs(z)./delta));
    z=z-g_ee_tmp;
    g_ee_inc=C*(1+ones(num_e,1)*tanh(v_e_tmp'/10)).*y ...
            +C*(1+tanh(v_e_tmp/10)*ones(1,num_e)).*z;     
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
    a=a+dt*a_inc;
    g_ee=g_ee+dt*g_ee_inc;
	v_i=v_i+dt*v_i_inc;
    m_i=m_i_inf(v_i); h_i=h_i+dt*h_i_inc; n_i=n_i+dt*n_i_inc; 
    q_i=q_i+dt*q_i_inc;
    s_i=s_i+dt*s_i_inc;
    

	% Determine which and how many e- and i-cells spiked in the current 
    % time step:

    which_e=find(v_e_old>-20 & v_e <=-20); which_i=find(v_i_old>-20 & v_i <=-20);
    l_e=length(which_e); l_i=length(which_i);
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
    
    lfp(k+1)=mean(v_e);

end;

% plot the spike rastergram

figure(1);
rastergram;

% turn the two-dimensional array g_ee into a one-dimensional one, 
% skipping the diagonal entries. 

vec_g=zeros(num_e*(num_e-1),1);
ij=0;
for i=1:num_e,
    for j=1:i-1,
        ij=ij+1;
        vec_g(ij)=g_ee(i,j);
    end;
    for j=i+1:num_e,
        ij=ij+1;
        vec_g(ij)=g_ee(i,j);
    end;
end;
R=max(vec_g);
sigma=10^(-4);
g=-0.5*R+[1:300]'/300*2*R-R/300;
density=zeros(length(g),1);
for ijk=1:num_e*(num_e-1),
    density=density+ ...
    exp(-(g-vec_g(ijk)).^2/(2*sigma^2))/sqrt(2*pi*sigma^2);
end;
density=density/num_e/(num_e-1);
figure(2);
subplot(111);
plot(g,density,'-k','Linewidth',5);
set(gca,'Fontsize',18);
xlabel('$\overline{g}_{EE}$  [mS/cm$^2$]','Fontsize',24);
ylabel('density [connections per mS/cm$^2$]', ...
    'Fontsize',24);
%set(gca,'Ytick',[]);
axis([0,2.5*10^(-3),0,max(density)*1.2]);

shg;

toc

