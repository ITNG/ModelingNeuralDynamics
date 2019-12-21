clf; clear;
tau_m=10; 
I=0.13;
w_k=0.05;
tilde_w_k=0.08;
tau_w=40;

dt=0.01;
dt05=dt/2;
t_final=100;
m_steps=round(t_final/dt);

v=zeros(m_steps+1,1);

for k=1:m_steps,
    t=(k-1)*dt;
    v_inc=-v(k)/tau_m+I-w_k*exp(-t/tau_w)*v(k);
    v_tmp=v(k)+dt05*v_inc;
    v_inc=-v_tmp/tau_m+I-w_k*exp(-(t+dt05)/tau_w)*v_tmp;
    v(k+1)=v(k)+dt*v_inc;
end;
subplot(211);
plot([0:m_steps]*dt,v,'-k','Linewidth',2);
set(gca,'Fontsize',16);
xlabel('$t$','Fontsize',20); 
title('$v$ (solid) and $\tilde{v}$ (dashes)','Fontsize',20);

for k=1:m_steps,
    t=(k-1)*dt;
    v_inc=-v(k)/tau_m+I-tilde_w_k*exp(-t/tau_w)*v(k);
    v_tmp=v(k)+dt05*v_inc;
    v_inc=-v_tmp/tau_m+I-tilde_w_k*exp(-(t+dt05)/tau_w)*v_tmp;
    v(k+1)=v(k)+dt*v_inc;
end;
hold on;
plot([0:m_steps]*dt,v,'--k','Linewidth',2);
axis([0,60,0,1]);
hold off;
shg;
