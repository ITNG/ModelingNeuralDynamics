clear; clf; rng('default'); rng(63806);
tic;

% Define network parameters: 

num_e=2; num_i=1;
i_ext_e=zeros(2,1); i_ext_e(1)=0.4; i_ext_e(2)=0.8; 
i_ext_i=zeros(1,1);
g_ee=zeros(num_e,num_e); g_ei=zeros(num_e,num_i);
g_ie=zeros(num_i,num_e); g_ii=zeros(num_i,num_i);

g_ei(:,1)=0.125*ones(2,1);
g_ie(1,:)=0.25*ones(1,2);
g_ii(1,1)=0.25;

v_rev_e=0; v_rev_i=-75; 
tau_r_e=0.5; tau_peak_e=0.5; tau_d_e=3; 
tau_r_i=0.5; tau_peak_i=0.5; tau_d_i=9; 
tau_dq_e=tau_d_q_function(tau_d_e,tau_r_e,tau_peak_e);
tau_dq_i=tau_d_q_function(tau_d_i,tau_r_i,tau_peak_i);

t_final=500;    % Time (in ms) simulated. 
dt=0.01;        % Time step used in solving the differential equations.
dt05=dt/2; 
m_steps=round(t_final/dt);

% STDP parameters

g_ee(1,2)=0.05; % initial strengths of recurrent excitatory
g_ee(2,1)=0.05; % connections.
C=1.45;     % constant that appears in the approximate delta function
            % C*(1+tanh(v_e/10)) (with v_e=membrane potential of an E-cell)
K_plus=g_ee;       % maximum increase in E-to-E synaptic strength 
                     % resulting from a single spike pair
K_minus=g_ee*2/3;     % maximum decrease in E-to-E synaptic strength
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

delta=max(g_ee/2,10^(-6));  % The changes resulting from STDP are constrained
                            % so the strength does not fall below zero, and
                            % does not rise above B. The "max" and "min"
                            % functions appearing when this constraint 
                            % is implemented are smoothed, with a smoothing
                            % parameter of delta.
                                     

% initialize dynamic variables

v_e=ones(num_e,1);
v_e(1)=-70; v_e(2)=-70;
m_e=m_e_inf(v_e); h_e=h_e_inf(v_e); n_e=n_e_inf(v_e); 
z=zeros(num_e,1); q_e=z; s_e=z;
a=100*ones(num_e,1);

v_i=-75*ones(num_i,1); m_i=m_i_inf(v_i); h_i=h_i_inf(v_i); n_i=n_i_inf(v_i); 
z=zeros(num_i,1); q_i=z; s_i=z; 

% solve the system of Hodgkin-Huxley-like equations using the midpoint method

num_spikes_e=0; t_e_spikes=[]; i_e_spikes=[];
num_spikes_i=0; t_i_spikes=[]; i_i_spikes=[];

for k=1:m_steps,
    t_old=(k-1)*dt; t_new=k*dt; t_mid=(t_old+t_new)/2;
	
	v_e_inc=0.1*(-67-v_e)+80*n_e.^4.*(-100-v_e)+100*m_e.^3.*h_e.*(50-v_e) ...
               +(g_ee'*s_e).*(v_rev_e-v_e)+(g_ie'*s_i).*(v_rev_i-v_e) ...
               +i_ext_e;
	n_e_inc=(n_e_inf(v_e)-n_e)./tau_n_e(v_e);
    h_e_inc=(h_e_inf(v_e)-h_e)./tau_h_e(v_e);
    q_e_inc=(1+tanh(v_e/10))/2.*(1-q_e)/0.1-q_e./tau_dq_e;
    s_e_inc=q_e.*(1-s_e)./tau_r_e-s_e./tau_d_e;
    a_inc=1-5*a.*(1+tanh(v_e/10));
    
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
    a_inc=1-5*a_tmp.*(1+tanh(v_e_tmp/10));
    y=g_ee_tmp+ ...
        diag(exp(-a_tmp/tau_plus).*(1-exp(-5*a_tmp/tau_plus)))*K_plus;
    y=min(B,y)-delta/2.*log(1+exp(-2*abs(B-y)./delta));
    y=y-g_ee_tmp;
    z=g_ee_tmp-K_minus* ...
        diag(exp(-a_tmp/tau_minus).*(1-exp(-5*a_tmp/tau_minus)));
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

    which_e=find(v_e_old<-40 & v_e >=-40); which_i=find(v_i_old<-40 & v_i >=-40);
    l_e=length(which_e); l_i=length(which_i);
    if l_e>0, 
        range=num_spikes_e+1:num_spikes_e+l_e; i_e_spikes(range)=which_e; 
        t_e_spikes(range)= ...
            ((v_e(which_e)+40)*(k-1)*dt+(-v_e_old(which_e)-40)*k*dt)./ ...
                (v_e(which_e)-v_e_old(which_e));
        num_spikes_e=num_spikes_e+l_e;
    end;
    if l_i>0, 
        range=num_spikes_i+1:num_spikes_i+l_i; i_i_spikes(range)=which_i; 
        t_i_spikes(range)= ...
            ((40+v_i(which_i))*(k-1)*dt+(-v_i_old(which_i)-40)*k*dt)./ ...
                (v_i(which_i)-v_i_old(which_i));
        num_spikes_i=num_spikes_i+l_i;
    end 
    
    g_12(k)=g_ee(1,2);
    g_21(k)=g_ee(2,1);

end;

% plot the spike rastergram

subplot(411);

if num_spikes_i>0, plot(t_i_spikes,i_i_spikes,'.b','Markersize',20); hold on; end;
if num_spikes_e>0, plot(t_e_spikes,i_e_spikes+num_i,'.r','Markersize',20);  hold on; end; 
plot([0,t_final],[num_i+1/2,num_i+1/2],'--k','Linewidth',1);
hold off;
set(gca,'Fontsize',12);
set(gca,'Ytick',[]);
axis([0,t_final,0,4]);

subplot(412);
t=[1:m_steps]*dt;
plot(t,g_12,'-k','Linewidth',2);
axis([0,t_final,0,B(1,2)]);
set(gca,'Fontsize',12);
ylabel('$\overline{g}_{{\rm syn},EE,12}$','Fontsize',18);

subplot(413);
plot(t,g_21,'-k','Linewidth',2);
axis([0,t_final,0,B(2,1)]);
set(gca,'Fontsize',12);
ylabel('$\overline{g}_{{\rm syn},EE,21}$','Fontsize',18);

subplot(414);
ind1=find(i_e_spikes==1);
t1=t_e_spikes(ind1);
ind2=find(i_e_spikes==2);
t2=t_e_spikes(ind2);
done=0;
for j=1:length(t2),
    ind=find(t1>t2(j));
    if length(ind)>0,
        done=done+1;
        d(j)=min(t1(ind))-t2(j);
    end;
end;
plot(t2(1:length(d)),d,'.k','Markersize',20);
set(gca,'Fontsize',12);
axis([0,t_final,0,12]);
xlabel('$t$ [ms]','Fontsize',18);
ylabel('lag','Fontsize',18);
    


shg;



toc

