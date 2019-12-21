clear; clf;
c=1.3;
g_k=23; 
g_na=30; 
g_h=12; 
g_A=22;
g_l=0.05;
v_k=-100;
v_na=90;
v_l=-70;
v_h=-32.9; 
v_A=-90; 


i_ext=0;


t_final=500;
dt=0.01;
dt05=dt/2;
m_steps=round(t_final/dt);

z=zeros(m_steps+1,1);
v=z; m=z; h=z; n=z; r=z; a=z; b=z;

v(1)=-63;
m(1)=m_o_inf(v(1));
h(1)=h_o_inf(v(1));
n(1)=n_o_inf(v(1));
r(1)=r_o_inf(v(1));
a(1)=a_o_inf(v(1));
b(1)=b_o_inf(v(1));


for k=1:m_steps,
    
    
    v_inc=(g_k*n(k)^4*(v_k-v(k))+g_na*m(k)^3*h(k)*(v_na-v(k)) ...
         +g_l*(v_l-v(k))+g_h*r(k)*(v_h-v(k))+g_A*a(k)*b(k)*(v_A-v(k)) ...
         + i_ext)/c;
    h_inc=alpha_h_o(v(k))*(1-h(k))-beta_h_o(v(k))*h(k); 
    n_inc=alpha_n_o(v(k))*(1-n(k))-beta_n_o(v(k))*n(k); 
    r_inc=(r_o_inf(v(k))-r(k))/tau_r_o(v(k));
    a_inc=(a_o_inf(v(k))-a(k))/tau_a_o(v(k));
    b_inc=(b_o_inf(v(k))-b(k))/tau_b_o(v(k));

    v_tmp=v(k)+dt05*v_inc;
    m_tmp=m_o_inf(v_tmp);
    h_tmp=h(k)+dt05*h_inc;
    n_tmp=n(k)+dt05*n_inc;
    r_tmp=r(k)+dt05*r_inc;
    a_tmp=a(k)+dt05*a_inc;
    b_tmp=b(k)+dt05*b_inc;
    
    
    v_inc=(g_k*n_tmp^4*(v_k-v_tmp)+g_na*m_tmp^3*h_tmp*(v_na-v_tmp) ...
         +g_l*(v_l-v_tmp)+g_h*r_tmp*(v_h-v_tmp) ...
         +g_A*a_tmp*b_tmp*(v_A-v_tmp)+i_ext)/c;
    m_inc=alpha_m_o(v_tmp)*(1-m_tmp)-beta_m_o(v_tmp)*m_tmp;
    h_inc=alpha_h_o(v_tmp)*(1-h_tmp)-beta_h_o(v_tmp)*h_tmp; 
    n_inc=alpha_n_o(v_tmp)*(1-n_tmp)-beta_n_o(v_tmp)*n_tmp; 
    r_inc=(r_o_inf(v_tmp)-r_tmp)/tau_r_o(v_tmp);
    a_inc=(a_o_inf(v_tmp)-a_tmp)/tau_a_o(v_tmp);
    b_inc=(b_o_inf(v_tmp)-b_tmp)/tau_b_o(v_tmp);
    
    v(k+1)=v(k)+dt*v_inc;
    m(k+1)=m_o_inf(v(k+1));
    h(k+1)=h(k)+dt*h_inc;
    n(k+1)=n(k)+dt*n_inc;
    r(k+1)=r(k)+dt*r_inc;
    a(k+1)=a(k)+dt*a_inc;
    b(k+1)=b(k)+dt*b_inc;
    
    
end;


t=(0:m_steps)*dt;
subplot(211);
plot(t,v,'-k','Linewidth',2);
shg;
hold off;
 set(gca,'Fontsize',16);
ylabel('$v$','Fontsize',20);
axis([0,t_final,-100,90]);

subplot(212);
plot(t,a.*b,'-k','Linewidth',2);
shg;
set(gca,'Fontsize',16);
ylabel('$ab$','Fontsize',20);
axis([0,t_final,0.,0.05]);

    
