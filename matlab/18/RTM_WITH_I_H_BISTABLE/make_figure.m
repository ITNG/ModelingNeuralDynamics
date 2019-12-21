clear; clf;
c=1;
g_k=80; 
g_na=100;
g_l=0.1;
v_k=-100;
v_na=50;
v_l=-67;

g_h=1;
v_h=-32.9;

i_ext=-3.19; 




f=@(v) g_na*m_inf(v).^3.*h_inf(v).*(v_na-v)+g_k*n_inf(v).^4.*(v_k-v)+ ...
       g_l*(v_l-v)+g_h*r_inf(v).*(v_h-v)+i_ext;
   
v=-100+[0:100000]/100000*150;
kl=[1:length(v)-1]; kr=[2:length(v)];
ind=find(f(v(kl)).*f(v(kr))<=0);
v_left=min(v(ind));
v_right=min(v(ind+1));
while v_right-v_left>10^(-12),
    v_c=(v_left+v_right)/2;
    if f(v_c)*f(v_left)<=0,
        v_right=v_c;
    else
        v_left=v_c;
    end;
end;
v_star=(v_left+v_right)/2;
m_star=m_inf(v_star);
h_star=h_inf(v_star);
n_star=n_inf(v_star);
r_star=r_inf(v_star);


t_final=3000; 
dt=0.01;
dt05=dt/2;
m_steps=round(t_final/dt);

z=zeros(m_steps+1,1);
v=z; m=z; h=z; n=z; r=z;

v(1)=v_star+0.01;
m(1)=m_star;
h(1)=h_star;
n(1)=n_star;
r(1)=r_star;
    
for k=1:m_steps,

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
    m_inc=alpha_m(v_tmp)*(1-m_tmp)-beta_m(v_tmp)*m_tmp;
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

v(1)=v_star+1;
m(1)=m_star;
h(1)=h_star;
n(1)=n_star;
r(1)=r_star;
    
for k=1:m_steps,

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
    m_inc=alpha_m(v_tmp)*(1-m_tmp)-beta_m(v_tmp)*m_tmp;
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
subplot(212);
plot(t,v,'-k','Linewidth',2);
set(gca,'Fontsize',16);
xlabel('$t$ [ms]','Fontsize',20);
ylabel('$v$ [mV]','Fontsize',20);