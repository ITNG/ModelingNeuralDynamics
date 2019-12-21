

clear; clf;

c=1;
g_k=9; 
g_na=35; 
g_l=0.1;
v_k=-90;
v_na=55;
v_l=-65;

U=0.5;
C=1.25;

tau_rec=500;
tau_d_q=5;
tau_r=3;
tau_d=9;

i_ext=0.5; 

t_final=200;
dt=0.01;
dt05=dt/2;
m_steps=round(t_final/dt);

z=zeros(m_steps+1,1);
v=z; m=z; h=z; n=z; p=z; q=z; s=z;

v(1)=-70;
m(1)=m_inf(v(1));
h(1)=h_inf(v(1));
n(1)=n_inf(v(1));
p(1)=1;
q(1)=0;
s(1)=0;


num_spikes=0;


for k=1:m_steps,
    
    v_inc=(g_k*n(k)^4*(v_k-v(k))+g_na*m(k)^3*h(k)*(v_na-v(k))+g_l*(v_l-v(k))+i_ext)/c;
    n_inc=alpha_n(v(k))*(1-n(k))-beta_n(v(k))*n(k);
    h_inc=alpha_h(v(k))*(1-h(k))-beta_h(v(k))*h(k);
    p_inc=-C*(1+tanh(v(k)/10))*p(k)*log(1/(1-U))+(1-p(k)-q(k))/tau_rec;
    q_inc= C*(1+tanh(v(k)/10))*p(k)*log(1/(1-U))-q(k)/tau_d_q;
    s_inc=q(k)*(1-s(k))/tau_r-s(k)/tau_d;
    
    v_tmp=v(k)+dt05*v_inc;
    h_tmp=h(k)+dt05*h_inc;
    n_tmp=n(k)+dt05*n_inc;
    m_tmp=m_inf(v_tmp);
    p_tmp=p(k)+dt05*p_inc;
    q_tmp=q(k)+dt05*q_inc;
    s_tmp=s(k)+dt05*s_inc;
    
    v_inc=(g_k*n_tmp^4*(v_k-v_tmp)+g_na*m_tmp^3*h_tmp*(v_na-v_tmp)+g_l*(v_l-v_tmp)+i_ext)/c;
    h_inc=alpha_h(v_tmp)*(1-h_tmp)-beta_h(v_tmp)*h_tmp;
    n_inc=alpha_n(v_tmp)*(1-n_tmp)-beta_n(v_tmp)*n_tmp;
    p_inc=-C*(1+tanh(v_tmp/10))*p_tmp*log(1/(1-U))+(1-p_tmp-q_tmp)/tau_rec;
    q_inc= C*(1+tanh(v_tmp/10))*p_tmp*log(1/(1-U))-q_tmp/tau_d_q;
    s_inc=q_tmp*(1-s_tmp)/tau_r-s_tmp/tau_d;
    
    v(k+1)=v(k)+dt*v_inc;
    h(k+1)=h(k)+dt*h_inc;
    n(k+1)=n(k)+dt*n_inc;
    m(k+1)=m_inf(v(k+1));
    p(k+1)=p(k)+dt*p_inc;
    q(k+1)=q(k)+dt*q_inc;
    s(k+1)=s(k)+dt*s_inc;
    
    if v(k+1)<-20 && v(k)>=-20,
        num_spikes=num_spikes+1;
    end;
    
end;

gamma=C*(1+tanh(v/10));
integral_of_delta_function_per_period=sum(gamma)*dt/num_spikes

t=[0:m_steps]*dt;
subplot(221);
set(gca,'Fontsize',14);
plot(t,v,'-k','Linewidth',2);
shg;
hold off;
%xlabel('$t$ [ms]','Fontsize',18); 
ylabel('$v$ [mV]','Fontsize',18);
axis([0,t_final,-100,50]);


subplot(222);
set(gca,'Fontsize',14);
plot(t,p,'-k','Linewidth',2);
shg;
hold off;
%xlabel('$t$ [ms]','Fontsize',18); 
ylabel('$p$','Fontsize',18);
axis([0,t_final,0,1]);


subplot(223);
set(gca,'Fontsize',14);
plot(t,q,'-k','Linewidth',2);
shg;
hold off;
xlabel('$t$ [ms]','Fontsize',18); 
ylabel('$q$','Fontsize',18);
axis([0,t_final,0,max(q)*1.1]);

subplot(224);
set(gca,'Fontsize',14);
plot(t,s,'-k','Linewidth',2);
shg;
hold off;
xlabel('$t$ [ms]','Fontsize',18); 
ylabel('$s$','Fontsize',18);
axis([0,t_final,0,max(s)*1.1]);

shg;

    
    
