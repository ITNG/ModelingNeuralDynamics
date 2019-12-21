

clear; clf; rng('default'); rng(63806);

global alpha Period g_bar m;

tau=10;
alpha=1; 
Period=25;
g_bar=0.1;
N=2000; m=mean(exp(alpha*cos(pi*[0:N-1]/N).^2)-1);



t_final=100;
dt=0.01; dt05=dt/2;
m_steps=round(t_final/dt);
t=[0:m_steps]'*dt;

subplot(311);
plot(t,g(t),'-b','Linewidth',2);
hold on;
plot([0,t_final],[g_bar,g_bar],'-r','Linewidth',2);
hold off;
set(gca,'Fontsize',14);
title('$g$ (blue) and $\overline{g}$ (red)','Fontsize',18);
axis([0,t_final,0,max(g(t))*1.2]);
text(-12,max(g(t))*1.2*0.5,'A','Fontsize',24);


v=zeros(m_steps+1,1);

subplot(312);
I=0.15;
k_old=1;
num_spikes=0;
for k=1:m_steps,
    v_inc=-v(k)/tau+I-g((k-1)*dt)*v(k);
    v_tmp=v(k)+dt05*v_inc;
    v_inc=-v_tmp/tau+I-g((k-1/2)*dt)*v_tmp;
    v(k+1)=v(k)+dt*v_inc;
    if v(k+1)>1,
        plot([k_old:k+1]*dt,v(k_old:k+1),'-b','Linewidth',2);
        k_old=k+1;
        v(k+1)=0;
        num_spikes=num_spikes+1;
        hold on;
        plot([k*dt,k*dt],[0,5],'-b','Linewidth',2);
    end;       
end;
plot([k_old:m_steps+1]*dt,v(k_old:m_steps+1),'-b','Linewidth',2);
set(gca,'Fontsize',14);
axis([0,t_final,0,6]);
text(-12,3,'B','Fontsize',24);
freq=num_spikes/t_final*1000


k_old=1;
num_spikes=0;
for k=1:m_steps,
    v_inc=-v(k)/tau+I-g_bar*v(k);
    v_tmp=v(k)+dt05*v_inc;
    v_inc=-v_tmp/tau+I-g_bar*v_tmp;
    v(k+1)=v(k)+dt*v_inc;
    if v(k+1)>1,
        plot([k_old:k+1]*dt,v(k_old:k+1),'-r','Linewidth',2);
        k_old=k+1;
        v(k+1)=0;
        num_spikes=num_spikes+1;
        hold on;
        plot([k*dt,k*dt],[0,1],'-r','Linewidth',2);
    end;       
end;
plot([k_old:m_steps+1]*dt,v(k_old:m_steps+1),'-r','Linewidth',2);
hold off;
%xlabel('$t$ [ms]','Fontsize',18);
title('$v$ (blue) and $\overline{v}$ (red), $I=0.15$','Fontsize',18);
freq_bar=num_spikes/t_final*1000

subplot(313);
I=0.2;
k_old=1;
num_spikes=0;
for k=1:m_steps,
    v_inc=-v(k)/tau+I-g((k-1)*dt)*v(k);
    v_tmp=v(k)+dt05*v_inc;
    v_inc=-v_tmp/tau+I-g((k-1/2)*dt)*v_tmp;
    v(k+1)=v(k)+dt*v_inc;
    if v(k+1)>1,
        plot([k_old:k+1]*dt,v(k_old:k+1),'-b','Linewidth',2);
        k_old=k+1;
        v(k+1)=0;
        num_spikes=num_spikes+1;
        hold on;
        plot([k*dt,k*dt],[0,5],'-b','Linewidth',2);
    end;       
end;
plot([k_old:m_steps+1]*dt,v(k_old:m_steps+1),'-b','Linewidth',2);
set(gca,'Fontsize',14);
axis([0,t_final,0,6]);
text(-12,3,'C','Fontsize',24);
freq=num_spikes/t_final*1000


k_old=1;
num_spikes=0;
for k=1:m_steps,
    v_inc=-v(k)/tau+I-g_bar*v(k);
    v_tmp=v(k)+dt05*v_inc;
    v_inc=-v_tmp/tau+I-g_bar*v_tmp;
    v(k+1)=v(k)+dt*v_inc;
    if v(k+1)>1,
        plot([k_old:k+1]*dt,v(k_old:k+1),'-r','Linewidth',2);
        k_old=k+1;
        v(k+1)=0;
        num_spikes=num_spikes+1;
        hold on;
        plot([k*dt,k*dt],[0,1],'-r','Linewidth',2);
    end;       
end;
plot([k_old:m_steps+1]*dt,v(k_old:m_steps+1),'-r','Linewidth',2);
hold off;
xlabel('$t$ [ms]','Fontsize',18);
title('$v$ (blue) and $\overline{v}$ (red), $I=0.2$','Fontsize',18);
freq_bar=num_spikes/t_final*1000


shg;