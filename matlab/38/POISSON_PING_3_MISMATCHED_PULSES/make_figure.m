clear; clf; rng('default'); rng(63806);
tic;

% Define network parameters: 

num_e=200; num_i=50;
sigma_e=0.00; i_ext_e=0.60*ones(num_e,1).*(1+sigma_e*randn(num_e,1)); 
sigma_i=0.00; i_ext_i=0.60*ones(num_i,1).*(1+sigma_i*randn(num_i,1));
g_hat_ee=0; g_hat_ei=1.25; g_hat_ie=1.25; g_hat_ii=0.4; 
p_ee=1.0; p_ei=1.0; p_ie=1.0; p_ii=1.0;
        % See explanation below to understand what the preceding eight
        % parameters mean.
v_rev_e=0; v_rev_i=-75; 
tau_r_e=0.3; tau_peak_e=0.3; tau_d_e=3; 
tau_r_i=0.3; tau_peak_i=0.3; tau_d_i=9; 
f_stoch=40;  g_stoch=0.1; % These are the parameters defining the Poisson
                          % streams of excitatory synaptic inputs into the
                          % E-cells: They have mean frequency f_stoch, 
                          % and maximum conductance g_stoch.
t_final=500;    % Time (in ms) simulated. 
dt=0.01;        % Time step used in solving the differential equations.

nn=5;    % number of E-cells receiving extra pulsatile input
ave=0.1; % temporal average of extra pulsatile input
P=29;    % period of extra pulsatile input
alpha=4; % parameter determining peakedness of pulsatile input
s=(1:1000)/100; 
mu=mean(exp(alpha*cos(pi*s).^2)-1);
phi_vec=(1:9)/10; 
         % vector of phases at which the pulse train is to be injected



% Process network parameters a bit:

u_ee=rand(num_e,num_e); u_ei=rand(num_e,num_i);
u_ie=rand(num_i,num_e); u_ii=rand(num_i,num_i);
g_ee=g_hat_ee*(u_ee<p_ee)/(num_e*p_ee); g_ei=g_hat_ei*(u_ei<p_ei)/(num_e*p_ei);
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


for ijk=1:length(phi_vec),
    phi=phi_vec(ijk);
    % initialize dynamic variables

    iv=rtm_init(i_ext_e,rand(num_e,1));
    v_e=iv(:,1); m_e=m_e_inf(v_e); h_e=iv(:,2); n_e=iv(:,3); 
    z=zeros(num_e,1); q_e=z; s_e=z; 
    q_stoch=zeros(num_e,1); s_stoch=zeros(num_e,1);

    iv=wb_init(i_ext_i,rand(num_i,1));
    v_i=iv(:,1); m_i=m_i_inf(v_i); h_i=iv(:,2); n_i=iv(:,3); 
    z=zeros(num_i,1); q_i=z; s_i=z; 


    % solve the system of Hodgkin-Huxley-like equations using the midpoint method

    num_spikes_e=0; t_e_spikes=[]; i_e_spikes=[];
    num_spikes_i=0; t_i_spikes=[]; i_i_spikes=[];

    for k=1:m_steps,
        t_old=(k-1)*dt; t_new=k*dt; t_mid=(t_old+t_new)/2;

        v_e_inc=0.1*(-67-v_e)+80*n_e.^4.*(-100-v_e)+100*m_e.^3.*h_e.*(50-v_e) ...
                   +(g_ee'*s_e).*(v_rev_e-v_e)+(g_ie'*s_i).*(v_rev_i-v_e) ...
                   +i_ext_e+g_stoch.*s_stoch.*(v_rev_e-v_e);
        v_e_inc(1:nn)=v_e_inc(1:nn)+ave*(exp(alpha*cos(pi*((k-1)*dt/P-phi)).^2)-1)/mu;
        n_e_inc=(n_e_inf(v_e)-n_e)./tau_n_e(v_e);
        h_e_inc=(h_e_inf(v_e)-h_e)./tau_h_e(v_e);
        q_e_inc=(1+tanh(v_e/10))/2.*(1-q_e)/0.1-q_e./tau_dq_e;
        s_e_inc=q_e.*(1-s_e)./tau_r_e-s_e./tau_d_e;
        q_stoch_inc=-q_stoch./tau_dq_e;
        s_stoch_inc=q_stoch.*(1-s_stoch)./tau_r_e-s_stoch./tau_d_e;
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
        q_stoch_tmp=q_stoch+dt05*q_stoch_inc;
        s_stoch_tmp=s_stoch+dt05*s_stoch_inc;
        v_i_tmp=v_i+dt05*v_i_inc;
        n_i_tmp=n_i+dt05*n_i_inc;
        m_i_tmp=m_i_inf(v_i_tmp);
        h_i_tmp=h_i+dt05*h_i_inc;
        q_i_tmp=q_i+dt05*q_i_inc;   
        s_i_tmp=s_i+dt05*s_i_inc;    

        v_e_inc=0.1*(-67-v_e_tmp)+80*n_e_tmp.^4.*(-100-v_e_tmp)+100*m_e_tmp.^3.*h_e_tmp.*(50-v_e_tmp) ...
                   +(g_ee'*s_e_tmp).*(v_rev_e-v_e_tmp)+(g_ie'*s_i_tmp).*(v_rev_i-v_e_tmp) ...
                   +i_ext_e+g_stoch.*s_stoch_tmp.*(v_rev_e-v_e_tmp);
        v_e_inc(1:nn)=v_e_inc(1:nn)+ave*(exp(alpha*cos(pi*((k-1/2)*dt/P-phi)).^2)-1)/mu;
        n_e_inc=(n_e_inf(v_e_tmp)-n_e_tmp)./tau_n_e(v_e_tmp);
        h_e_inc=(h_e_inf(v_e_tmp)-h_e_tmp)./tau_h_e(v_e_tmp);
        q_e_inc=(1+tanh(v_e_tmp/10))/2.*(1-q_e_tmp)/0.1-q_e_tmp./tau_dq_e;
        s_e_inc=q_e_tmp.*(1-s_e_tmp)./tau_r_e-s_e_tmp./tau_d_e;
        q_stoch_inc=-q_stoch_tmp./tau_dq_e;
        s_stoch_inc=q_stoch_tmp.*(1-s_stoch_tmp)./tau_r_e-s_stoch_tmp./tau_d_e;
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
        q_stoch=q_stoch+dt*q_stoch_inc;
        s_stoch=s_stoch+dt*s_stoch_inc;
        v_i=v_i+dt*v_i_inc;
        m_i=m_i_inf(v_i); h_i=h_i+dt*h_i_inc; n_i=n_i+dt*n_i_inc; 
        q_i=q_i+dt*q_i_inc;
        s_i=s_i+dt*s_i_inc;


        % Determine which and how many e- and i-cells spiked in the current 
        % time step:

        which_e=find(v_e_old>-20 & v_e <=-20); 
        which_i=find(v_i_old>-20 & v_i <=-20);
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

        % Random input arrivals: 

        u=rand(num_e,1);
        ind=(u<f_stoch/1000*dt);
        % ind(i)=1 with probability f_stoch/1000*dt
        q_stoch(ind)=1;

        spike_count(k)=l_e;

    end;
    ind=find(i_e_spikes<=nn);
    spike_count=length(ind)
    spike_count_vec(ijk)=spike_count;
    ijk
    phi
end;

subplot(111);
plot(phi_vec,spike_count_vec,'ok','Markersize',10);
set(gca,'Fontsize',24);
xlabel('$\varphi$','Fontsize',28);
nn_string=num2str(nn);
ylabel(['\# spikes in E-cells 1$-$', nn_string],'Fontsize',28);

P=polyfit(phi_vec,spike_count_vec,1);
A=P(1); 
B=P(2);
hold on;
plot(phi_vec,A*phi_vec+B,'-g','Linewidth',4);
hold off;
axis([0,1,10,40]); axis('square');
shg;

toc

