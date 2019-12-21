clear; clf;

% This is very similar to the code in HH_F_I_CURVE. See there for 
% more extensive comments. 

tic 
    
c=1;
g_k=80; 
g_na=100;
g_l=0.1;
v_k=-100;
v_na=50;
v_l=-67;
g_m=0.2;

i_ext=0.506;




v=-100+[0:100000]/100000*150;
f=@(v) g_na*m_inf(v).^3.*h_inf(v).*(v_na-v)+...
       g_k*n_inf(v).^4.*(v_k-v)+...
       g_l*(v_l-v)+...
       g_m*w_inf(v).*(v_k-v)+i_ext;
   
l=length(v);
kl=[1:l-1]; kr=[2:l];
ind=find(f(v(kl)).*f(v(kr))<=0);
ind=min(ind);
v_left=v(ind); 
v_right=v(ind+1);
while (v_right-v_left>10^(-14)),
    v_c=(v_right+v_left)/2;
    if f(v_c)*f(v_right)<=0,
        v_left=v_c;
    else
        v_right=v_c;
    end;
end;
v_star=(v_left+v_right)/2
m_star=m_inf(v_star)
h_star=h_inf(v_star)
n_star=n_inf(v_star)
w_star=w_inf(v_star)


t_final=500;
dt=0.01; dt05=dt/2;
m_steps=round(t_final/dt);

z=zeros(m_steps+1,1);
v=z; m=z; h=z; n=z;


v(1)=v_star+1; 
m(1)=m_star;
h(1)=h_star;
n(1)=n_star;
w(1)=w_star;

for k=1:m_steps,
    
    v_inc=(g_na*m(k)^3*h(k)*(v_na-v(k))+ ...
        g_k*n(k)^4*(v_k-v(k))+g_l*(v_l-v(k))+ ...
        g_m*w(k)*(v_k-v(k))+i_ext)/c;
    h_inc=alpha_h(v(k))*(1-h(k))-beta_h(v(k))*h(k);
    n_inc=alpha_n(v(k))*(1-n(k))-beta_n(v(k))*n(k);
    w_inc=(w_inf(v(k))-w(k))/tau_w(v(k));
    
    v_tmp=v(k)+dt05*v_inc;
    m_tmp=m_inf(v_tmp);
    h_tmp=h(k)+dt05*h_inc;
    n_tmp=n(k)+dt05*n_inc;
    w_tmp=w(k)+dt05*w_inc;
    
    v_inc=(g_na*m_tmp^3*h_tmp*(v_na-v_tmp)+ ...
        g_k*n_tmp^4*(v_k-v_tmp)+g_l*(v_l-v_tmp)+ ...
        g_m*w_tmp*(v_k-v_tmp)+i_ext)/c;
    h_inc=alpha_h(v_tmp)*(1-h_tmp)-beta_h(v_tmp)*h_tmp;
    n_inc=alpha_n(v_tmp)*(1-n_tmp)-beta_n(v_tmp)*n_tmp;
    w_inc=(w_inf(v_tmp)-w(k))/tau_w(v_tmp);
    
    v(k+1)=v(k)+dt*v_inc;
    m(k+1)=m_inf(v(k+1));
    h(k+1)=h(k)+dt*h_inc;
    n(k+1)=n(k)+dt*n_inc;
    w(k+1)=w(k)+dt*w_inc;
    
end;


subplot(311);
t=[0:m_steps]*dt;
plot(t,h,'-k','Linewidth',2);
hold on;
plot(t,h_inf(v),'--b','Linewidth',2);
plot([0,t_final],[h_star,h_star],'-r','Linewidth',2);
hold off;
set(gca,'Fontsize',14);
ylabel('$h$','Fontsize',18);
axis([0,t_final,0.95,1]);

subplot(312);
t=[0:m_steps]*dt;
plot(t,n,'-k','Linewidth',2);
hold on;
plot(t,n_inf(v),'--b','Linewidth',2);
plot([0,t_final],[n_star,n_star],'-r','Linewidth',2);
hold off;
set(gca,'Fontsize',14);
ylabel('$n$','Fontsize',18);
axis([0,t_final,0,0.2]);

subplot(313);
t=[0:m_steps]*dt;
plot(t,w,'-k','Linewidth',2);
hold on;
plot(t,w_inf(v),'--b','Linewidth',2);
plot([0,t_final],[w_star,w_star],'-r','Linewidth',2);
ind=find(w<w_star);
plot(t(ind),0.045,'.m','Markersize',20);
hold off;
set(gca,'Fontsize',14);
xlabel('$t$ [ms]','Fontsize',18);
ylabel('$w$','Fontsize',18);
axis([0,t_final,0.045,0.06]);
shg;

