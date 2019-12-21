

clear; clf;
c=1;
g_k=224; 
g_na=112;
g_l=0.5;
v_k=-90;
v_na=60;
v_l=-70;

i_ext=6.9;

v_vec=-100+[0:500000]/500000*150;
f=@(v) (g_k*n_inf(v).^2.*(v_k-v)+g_na*m_inf(v).^3.*h_inf(v).*(v_na-v)+ ...
    g_l*(v_l-v)+i_ext)/c;

k=[1:length(v_vec)-1]; k_p=[2:length(v_vec)];
ind=find(f(v_vec(k)).*f(v_vec(k_p))<0);
ind=min(ind);
v_star=(v_vec(ind)+v_vec(ind+1))/2;
m_star=m_inf(v_star);
h_star=h_inf(v_star);
n_star=n_inf(v_star);





t_final=40;
dt=0.01;
dt05=dt/2;
m_steps=round(t_final/dt);

z=zeros(m_steps+1,1);
v=z; m=z; h=z; n=z;


t=[0:m_steps]*dt;

 
v(1)=v_star+3;
n(1)=n_star;
m(1)=m_star;
h(1)=h_star;


for k=1:m_steps,
    
    v_inc=(g_k*n(k)^2*(v_k-v(k))+g_na*m(k)^3*h(k)*(v_na-v(k))+g_l*(v_l-v(k))+i_ext)/c;
    n_inc=alpha_n(v(k))*(1-n(k))-beta_n(v(k))*n(k);
    h_inc=alpha_h(v(k))*(1-h(k))-beta_h(v(k))*h(k);
    
    v_tmp=v(k)+dt05*v_inc;
    m_tmp=m_inf(v_tmp);
    h_tmp=h(k)+dt05*h_inc;
    n_tmp=n(k)+dt05*n_inc;
    
    v_inc=(g_k*n_tmp^2*(v_k-v_tmp)+g_na*m_tmp^3*h_tmp*(v_na-v_tmp)+g_l*(v_l-v_tmp)+i_ext)/c;
    m_inc=alpha_m(v_tmp)*(1-m_tmp)-beta_m(v_tmp)*m_tmp;
    h_inc=alpha_h(v_tmp)*(1-h_tmp)-beta_h(v_tmp)*h_tmp;
    n_inc=alpha_n(v_tmp)*(1-n_tmp)-beta_n(v_tmp)*n_tmp;
    
    v(k+1)=v(k)+dt*v_inc;
    m(k+1)=m_inf(v(k+1));
    h(k+1)=min(h_star,h(k)+dt*h_inc);
    n(k+1)=n(k)+dt*n_inc;
    
end;

t=[0:m_steps]*dt;
subplot(211);
plot(t,v,'-k','Linewidth',2);
axis([0,t_final,-95,55]);
shg;
hold off;
 set(gca,'Fontsize',16);
xlabel('$t$ [ms]','Fontsize',20); ylabel('$v$ [mV]','Fontsize',20);
    

