clear; clf; rng('default'); rng(63806);

global alpha Period g_bar m;

alpha=5;
Period=25;
g_bar=0.1;
N=2000; m=mean(exp(alpha*cos(pi*[0:N-1]/N).^2)-1);
tau=10;
I=0.1;
tau_noise=3;
sigma_noise=0.08;

t_final=500;
dt=0.01; dt05=dt/2;
m_steps=round(t_final/dt);
gamma=sigma_noise*sqrt(1-exp(-2*dt/tau_noise));

t=[0:m_steps]'*dt;

subplot(311);
plot(t,g(t),'-b','Linewidth',2);
set(gca,'Fontsize',14);
ylabel('$g$','Fontsize',18);
axis([0,t_final,0,0.5]);
text(-75,0.25,'A','Fontsize',24);
hold on;
plot(t,g_bar*ones(m_steps+1,1),'-r','Linewidth',2);
hold off;


v=zeros(m_steps+1,1);
s_noise=zeros(m_steps+1,1);

s_noise(1)=sigma_noise*randn(1,1);
for k=1:m_steps,
    s_noise(k+1)=s_noise(k)*exp(-dt/tau_noise)+gamma*randn(1,1);
end;

subplot(312);
k_old=1;
num_spikes=0;
for k=1:m_steps,
    v_inc=-v(k)/tau+I+s_noise(k) ...
        -g((k-1)*dt)*v(k);
    v_tmp=v(k)+dt05*v_inc;
    v_inc=-v_tmp/tau+I+(s_noise(k)+s_noise(k+1))/2 ...
        -g((k-1/2)*dt)*v_tmp;
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
hold off;
set(gca,'Fontsize',14);
axis([0,t_final,-1,6]);
text(-75,2.5,'B','Fontsize',24);
ylabel('$v$','Fontsize',18);
freq=num_spikes/t_final*1000


subplot(313);
k_old=1;
num_spikes=0;
for k=1:m_steps,
    v_inc=-v(k)/tau+I+s_noise(k)-g_bar*v(k);
    v_tmp=v(k)+dt05*v_inc;
    v_inc=-v_tmp/tau+I+(s_noise(k)+s_noise(k+1))/2-g_bar*v_tmp;
    v(k+1)=v(k)+dt*v_inc;
    if v(k+1)>1,
        plot([k_old:k+1]*dt,v(k_old:k+1),'-r','Linewidth',2);
        k_old=k+1;
        v(k+1)=0;
        num_spikes=num_spikes+1;
        hold on;
        plot([k*dt,k*dt],[0,5],'-r','Linewidth',2);
    end;       
end;
plot([k_old:m_steps+1]*dt,v(k_old:m_steps+1),'-r','Linewidth',2);
hold off;
set(gca,'Fontsize',14);
axis([0,t_final,-1,6]);
text(-75,2.5,'C','Fontsize',24);
ylabel('$\overline{v}$','Fontsize',18);
xlabel('$t$ [ms]','Fontsize',18);
freq_bar=num_spikes/t_final*1000


shg;