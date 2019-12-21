clear; clf;

c=1;
g_k=9; 
g_na=35;
g_l=0.1;
v_k=-90;
v_na=55;
v_l=-65;

num=2;
gap_strength=0.01;
g_gap=zeros(num,num);
g_gap(1,2)=gap_strength;
g_gap(2,1)=gap_strength;
d_gap=sum(g_gap)';

i_ext=[1; 0];

t_final=100;
dt=0.01;
dt05=dt/2;
m_steps=round(t_final/dt);

z=zeros(num,m_steps+1);
v=z; m=z; h=z; n=z;

v(:,1)=[-63,-63];
m(:,1)=m_inf(v(1));
h(:,1)=h_inf(v(1));
n(:,1)=n_inf(v(1));


for k=1:m_steps,
    
    v_inc=(g_k*n(:,k).^4.*(v_k-v(:,k))+g_na*m(:,k).^3.*h(:,k).*(v_na-v(:,k))+ ...
           g_l*(v_l-v(:,k))+ ...
           g_gap*v(:,k)-d_gap.*v(:,k)+i_ext)/c;
    n_inc=alpha_n(v(:,k)).*(1-n(:,k))-beta_n(v(:,k)).*n(:,k);
    h_inc=alpha_h(v(:,k)).*(1-h(:,k))-beta_h(v(:,k)).*h(:,k);
    
    v_tmp=v(:,k)+dt05*v_inc;
    h_tmp=h(:,k)+dt05*h_inc;
    n_tmp=n(:,k)+dt05*n_inc;
    m_tmp=m_inf(v_tmp);
    
    v_inc=(g_k*n_tmp.^4.*(v_k-v_tmp)+g_na*m_tmp.^3.*h_tmp.*(v_na-v_tmp)+ ...
           g_l*(v_l-v_tmp)+ ...
           g_gap*v_tmp-d_gap.*v_tmp+i_ext)/c;
    h_inc=alpha_h(v_tmp).*(1-h_tmp)-beta_h(v_tmp).*h_tmp;
    n_inc=alpha_n(v_tmp).*(1-n_tmp)-beta_n(v_tmp).*n_tmp;
    
    v(:,k+1)=v(:,k)+dt*v_inc;
    h(:,k+1)=h(:,k)+dt*h_inc;
    n(:,k+1)=n(:,k)+dt*n_inc;
    m(:,k+1)=m_inf(v(:,k+1));
    
end;

t=[0:m_steps]*dt;
subplot(211);
set(gca,'Fontsize',16);
plot(t,v(1,:),'-k','Linewidth',3);
ylabel('$v_1$ [mV]','Fontsize',20);
axis([0,t_final,-100,50]);
hold on;

subplot(212);
set(gca,'Fontsize',16);
plot(t,v(2,:),'-k','Linewidth',3);
xlabel('$t$ [ms]','Fontsize',20); ylabel('$v_2$ [mV]','Fontsize',20);
axis([0,t_final,min(v(2,round(m_steps/2):m_steps))-.75,max(v(2,round(m_steps/2):m_steps))+.75]);
hold on;

v(:,1)=[-63,-63];
m(:,1)=m_inf(v(1));
h(:,1)=h_inf(v(1));
n(:,1)=n_inf(v(1));
    
for k=1:m_steps,
    
    thr=-50;
    g_gap=[0 gap_strength; gap_strength 0]*(v(1,k)<=thr)+zeros(2,2)*(v(1,k)>=thr);
    d_gap=sum(g_gap)';
    
    v_inc=(g_k*n(:,k).^4.*(v_k-v(:,k))+g_na*m(:,k).^3.*h(:,k).*(v_na-v(:,k))+ ...
           g_l*(v_l-v(:,k))+ ...
           g_gap*v(:,k)-d_gap.*v(:,k)+i_ext)/c;
    n_inc=alpha_n(v(:,k)).*(1-n(:,k))-beta_n(v(:,k)).*n(:,k);
    h_inc=alpha_h(v(:,k)).*(1-h(:,k))-beta_h(v(:,k)).*h(:,k);
    
    v_tmp=v(:,k)+dt05*v_inc;
    h_tmp=h(:,k)+dt05*h_inc;
    n_tmp=n(:,k)+dt05*n_inc;
    m_tmp=m_inf(v_tmp);
    
    v_inc=(g_k*n_tmp.^4.*(v_k-v_tmp)+g_na*m_tmp.^3.*h_tmp.*(v_na-v_tmp)+ ...
           g_l*(v_l-v_tmp)+ ...
           g_gap*v_tmp-d_gap.*v_tmp+i_ext)/c;
    h_inc=alpha_h(v_tmp).*(1-h_tmp)-beta_h(v_tmp).*h_tmp;
    n_inc=alpha_n(v_tmp).*(1-n_tmp)-beta_n(v_tmp).*n_tmp;
    
    v(:,k+1)=v(:,k)+dt*v_inc;
    h(:,k+1)=h(:,k)+dt*h_inc;
    n(:,k+1)=n(:,k)+dt*n_inc;
    m(:,k+1)=m_inf(v(:,k+1));
    
end;


subplot(211);
plot(t,v(1,:),'-r','Linewidth',1);
hold off;

subplot(212);
plot(t,v(2,:),'-r','Linewidth',1);
hold off;

shg;
    
    
