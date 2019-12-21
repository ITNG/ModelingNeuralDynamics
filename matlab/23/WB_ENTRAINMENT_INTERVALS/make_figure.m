clear; clf;

% Warning: This takes quite a while to run. 
% To reassure you that the computer is working, 
% it will print out the values of g_syn as it
% processes them.

tic;
c=1;
g_k=9; 
g_na=35; 
g_l=0.1;
v_k=-90;
v_na=55;
v_l=-65;

i_ext=0;

T=50;                 % period of inputs 
N=200;
a=0.142; b=0.195;
g_syn_vec=a+[0:N]/N*(b-a);
tau_d=2; tau_r=0.5; tau_peak=tau_r;     % input kinetics
tau_d_q=tau_d_q_function(tau_d,tau_r,tau_peak);

t_final=10000;
dt=0.01;
dt05=dt/2;
m_steps=round(t_final/dt);

z=zeros(m_steps+1,1);
v=z; m=z; h=z; n=z;

for ijk=1:length(g_syn_vec),
    g_syn=g_syn_vec(ijk);
    g_syn
    v(1)=-65;
    m(1)=m_inf(v(1));
    h(1)=h_inf(v(1));
    n(1)=n_inf(v(1));
    q(1)=0;
    s(1)=0;

    for k=1:m_steps,

        t=(k-1)*dt;
        if abs(round(t/T)-t/T)<10^(-12) & k>1,
            q(k)=1;
        end;

        v_inc=(g_k*n(k)^4*(v_k-v(k))+g_na*m(k)^3*h(k)*(v_na-v(k))+ ...
               g_l*(v_l-v(k))-g_syn*s(k)*v(k)+i_ext)/c;
        h_inc=alpha_h(v(k))*(1-h(k))-beta_h(v(k))*h(k); 
        n_inc=alpha_n(v(k))*(1-n(k))-beta_n(v(k))*n(k); 
        s_inc=q(k)*(1-s(k))/tau_r-s(k)/tau_d;
        q_inc=-q(k)/tau_d_q;

        v_tmp=v(k)+dt05*v_inc;
        m_tmp=m_inf(v_tmp);
        h_tmp=h(k)+dt05*h_inc;
        n_tmp=n(k)+dt05*n_inc;
        s_tmp=s(k)+dt05*s_inc;
        q_tmp=q(k)+dt05*q_inc;

        v_inc=(g_k*n_tmp^4*(v_k-v_tmp)+g_na*m_tmp^3*h_tmp*(v_na-v_tmp)+ ...
               g_l*(v_l-v_tmp)-g_syn*s_tmp*v_tmp+i_ext)/c;
        m_inc=alpha_m(v_tmp)*(1-m_tmp)-beta_m(v_tmp)*m_tmp;
        h_inc=alpha_h(v_tmp)*(1-h_tmp)-beta_h(v_tmp)*h_tmp; 
        n_inc=alpha_n(v_tmp)*(1-n_tmp)-beta_n(v_tmp)*n_tmp; 
        s_inc=q_tmp*(1-s_tmp)/tau_r-s_tmp/tau_d;
        q_inc=-q_tmp/tau_d_q;

        v(k+1)=v(k)+dt*v_inc;
        m(k+1)=m_inf(v(k+1));
        h(k+1)=h(k)+dt*h_inc;
        n(k+1)=n(k)+dt*n_inc;
        s(k+1)=s(k)+dt*s_inc;
        q(k+1)=q(k)+dt*q_inc;

    end;


    t=[0:m_steps]'*dt;
    ind=find(v(1:m_steps)>=-20&v(2:m_steps+1)<-20&t(1:m_steps)>=t_final/2);
    t_spikes=(t(ind).*(-v(ind+1)-20)+t(ind+1).*(20+v(ind)))./(v(ind)-v(ind+1));
    num_spikes=length(t_spikes);
    periods=t_spikes(2:num_spikes)-t_spikes(1:num_spikes-1);
    sigma=(max(periods)-min(periods))/mean(periods);
    n=mean(periods)/T;
    n_vec(ijk)=n;
    sigma_vec(ijk)=sigma;
    
end;


subplot(111);
plot(g_syn_vec,n_vec,'.k','Markersize',5);
hold on;
ind=find(sigma_vec<10^(-3));
plot(g_syn_vec(ind),n_vec(ind),'.r','Markersize',15);
hold off;
axis([0.14,0.19,0,18]);
set(gca,'Fontsize',24);
xlabel('$\overline{g}_{\rm syn}$','Fontsize',32);
ylabel('$n$','Fontsize',32);
axis('square');
set(gca,'Xtick',[0.14:0.01:0.19]);
set(gca,'Ytick',[5:5:15])
shg;
toc;