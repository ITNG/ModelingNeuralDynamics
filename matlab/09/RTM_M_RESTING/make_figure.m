clear; clf;
c=1;
g_k=80; 
g_na=100;
g_l=0.1;
v_k=-100;
v_na=50;
v_l=-67;
g_m=0.25;

i_ext=0; 

t_final=600; 
dt=0.01;
dt05=dt/2;
m_steps=round(t_final/dt);

z=zeros(m_steps+1,1);
v=z; m=z; h=z; n=z; w=z;

v(1)=-70;
m(1)=m_inf(v(1));
h(1)=h_inf(v(1));
n(1)=n_inf(v(1));
w(1)=0;
    
for k=1:m_steps,

    v_inc=(g_k*n(k)^4*(v_k-v(k))+g_na*m(k)^3*h(k)*(v_na-v(k))+g_l*(v_l-v(k)) ...
          +g_m*w(k)*(v_k-v(k))+i_ext)/c;
    n_inc=alpha_n(v(k))*(1-n(k))-beta_n(v(k))*n(k);
    h_inc=alpha_h(v(k))*(1-h(k))-beta_h(v(k))*h(k);
    w_inc=(w_inf(v(k))-w(k))/tau_w(v(k));

    v_tmp=v(k)+dt05*v_inc;
    m_tmp=m_inf(v_tmp);
    h_tmp=h(k)+dt05*h_inc;
    n_tmp=n(k)+dt05*n_inc;
    w_tmp=w(k)+dt05*w_inc;

    v_inc=(g_k*n_tmp^4*(v_k-v_tmp)+g_na*m_tmp^3*h_tmp*(v_na-v_tmp)+g_l*(v_l-v_tmp) ...
         +g_m*w_tmp*(v_k-v_tmp)+i_ext)/c;
    m_inc=alpha_m(v_tmp)*(1-m_tmp)-beta_m(v_tmp)*m_tmp;
    h_inc=alpha_h(v_tmp)*(1-h_tmp)-beta_h(v_tmp)*h_tmp;
    n_inc=alpha_n(v_tmp)*(1-n_tmp)-beta_n(v_tmp)*n_tmp;
    w_inc=(w_inf(v_tmp)-w_tmp)/tau_w(v_tmp);

    v(k+1)=v(k)+dt*v_inc;
    m(k+1)=m_inf(v(k+1));
    h(k+1)=h(k)+dt*h_inc;
    n(k+1)=n(k)+dt*n_inc;
    w(k+1)=w(k)+dt*w_inc;

end;

t=[0:m_steps]*dt;
subplot(211);
plot(t,v,'-k','Linewidth',2);
shg;
set(gca,'Fontsize',16);
ylabel('$v$ [mV]','Fontsize',20);
subplot(212);
plot(t,w,'-k','Linewidth',2);
set(gca,'Fontsize',16);
xlabel('$t$ [ms]','Fontsize',20); 
ylabel('$w$','Fontsize',20);
axis([0,t_final,0,0.05]);

