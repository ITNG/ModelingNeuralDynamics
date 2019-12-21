clear;
c=1;
g_k=36; 
g_na=120;
g_l=0.3;
v_k=-82;
v_na=45;
v_l=-59;

i_ext=10;

t_final=50;
dt=0.01; dt05=dt/2;
m_steps=round(t_final/dt);

v=zeros(m_steps+1,1);
m=v; h=v; n=v;

v(1)=-50;
m(1)=alpha_m(v(1))/(alpha_m(v(1))+beta_m(v(1)));
h(1)=0.6;
n(1)=0.4;

subplot(111);
for k=1:m_steps,
    
    v_inc=(g_na*m(k)^3*h(k)*(v_na-v(k))+g_k*n(k)^4*(v_k-v(k))+g_l*(v_l-v(k))+i_ext)/c;
    m_inc=alpha_m(v(k))*(1-m(k))-beta_m(v(k))*m(k);
    h_inc=alpha_h(v(k))*(1-h(k))-beta_h(v(k))*h(k);
    n_inc=alpha_n(v(k))*(1-n(k))-beta_n(v(k))*n(k);
    
    v_tmp=v(k)+dt05*v_inc;
    m_tmp=m(k)+dt05*m_inc;
    h_tmp=h(k)+dt05*h_inc;
    n_tmp=n(k)+dt05*n_inc;
    
    v_inc=(g_na*m_tmp^3*h_tmp*(v_na-v_tmp)+g_k*n_tmp^4*(v_k-v_tmp)+g_l*(v_l-v_tmp)+i_ext)/c;
    m_inc=alpha_m(v_tmp)*(1-m_tmp)-beta_m(v_tmp)*m_tmp;
    h_inc=alpha_h(v_tmp)*(1-h_tmp)-beta_h(v_tmp)*h_tmp;
    n_inc=alpha_n(v_tmp)*(1-n_tmp)-beta_n(v_tmp)*n_tmp;
    
    v(k+1)=v(k)+dt*v_inc;
    m(k+1)=m(k)+dt*m_inc;
    h(k+1)=h(k)+dt*h_inc;
    n(k+1)=n(k)+dt*n_inc;
    
    if abs(k*dt-0.1)<10^(-8),
        vec=[v(k+1)-v(k); n(k+1)-n(k)];
        arrow(-100,50,0,1,v(k+1),n(k+1),vec,0.05,2,'-k')
    end;
    
    if abs(k*dt-1)<10^(-8),
        vec=[v(k+1)-v(k); n(k+1)-n(k)];
        arrow(-100,50,0,1,v(k+1),n(k+1),vec,0.05,2,'-k')
    end;
    
    if abs(k*dt-20)<10^(-8),
        vec=[v(k+1)-v(k-100); n(k+1)-n(k-100)];
        arrow(-100,50,0,1,v(k+1),n(k+1),vec,0.05,2,'-k')
    end;
end;

t=[0:m_steps]*dt;


plot(v,n,'-k','Linewidth',2)
hold off;
set(gca,'Fontsize',24);
xlabel('$v$ [mV]','Fontsize',32); ylabel('$n$','fontsize',32); 
axis([-100,50,0,1]); axis('square');
shg;
    
    
