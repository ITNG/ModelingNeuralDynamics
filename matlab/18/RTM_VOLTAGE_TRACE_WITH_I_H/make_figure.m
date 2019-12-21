clear; clf;
c=1;
g_k=80; 
g_na=100;
g_l=0.1;
v_k=-100;
v_na=50;
v_l=-67;

g_h=1.0;
v_h=-32.9;

t_final=300; 
dt=0.01;
dt05=dt/2;
m_steps=round(t_final/dt);

z=zeros(m_steps+1,1);
v=z; m=z; h=z; n=z; r=z;

v(1)=-70;
m(1)=0;
h(1)=1;
n(1)=0;
r(1)=0.5;
    
for k=1:m_steps,

    i_ext=-4*((k-1)*dt>=t_final/3&(k-1)*dt<2*t_final/3);
    
    v_inc=(g_k*n(k)^4*(v_k-v(k))+g_na*m(k)^3*h(k)*(v_na-v(k))+g_l*(v_l-v(k)) ...
          +g_h*r(k)*(v_h-v(k))+i_ext)/c;
    n_inc=alpha_n(v(k))*(1-n(k))-beta_n(v(k))*n(k);
    h_inc=alpha_h(v(k))*(1-h(k))-beta_h(v(k))*h(k);
    r_inc=(r_inf(v(k))-r(k))/tau_r(v(k));

    v_tmp=v(k)+dt05*v_inc;
    m_tmp=m_inf(v_tmp);
    h_tmp=h(k)+dt05*h_inc;
    n_tmp=n(k)+dt05*n_inc;
    r_tmp=r(k)+dt05*r_inc;

    v_inc=(g_k*n_tmp^4*(v_k-v_tmp)+g_na*m_tmp^3*h_tmp*(v_na-v_tmp)+g_l*(v_l-v_tmp) ...
         +g_h*r_tmp*(v_h-v_tmp)+i_ext)/c;
    h_inc=alpha_h(v_tmp)*(1-h_tmp)-beta_h(v_tmp)*h_tmp;
    n_inc=alpha_n(v_tmp)*(1-n_tmp)-beta_n(v_tmp)*n_tmp;
    r_inc=(r_inf(v_tmp)-r_tmp)/tau_r(v_tmp);

    v(k+1)=v(k)+dt*v_inc;
    m(k+1)=m_inf(v(k+1));
    h(k+1)=h(k)+dt*h_inc;
    n(k+1)=n(k)+dt*n_inc;
    r(k+1)=r(k)+dt*r_inc;

end;

t=[0:m_steps]*dt;
subplot(211);
plot(t,v,'-k','Linewidth',2);
set(gca,'Fontsize',16);
ylabel('$v$ [mV]','Fontsize',20);

shg;

