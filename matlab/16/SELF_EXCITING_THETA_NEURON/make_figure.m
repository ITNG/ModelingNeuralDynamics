clear; clf; 
I=-0.05;
w_max=0.20;
tau_w=20;

t_final=50;

dt=0.01; dt05=dt/2; m_steps=round(t_final/dt);

theta=zeros(m_steps+1,1); 
w=zeros(m_steps+1,1);
theta(1)=pi/2;
num_spikes=0;
for k=1:m_steps
    theta_inc=1-cos(theta(k))+(I+w(k))*(1+cos(theta(k)));
    w_inc=-w(k)/tau_w;
    theta_tmp=theta(k)+dt05*theta_inc;
    w_tmp=w(k)+dt05*w_inc;
    theta_inc=1-cos(theta_tmp)+(I+w_tmp)*(1+cos(theta_tmp));
    w_inc=-w_tmp/tau_w;
    theta(k+1)=theta(k)+dt*theta_inc;
    w(k+1)=w(k)+dt*w_inc;
    if theta(k+1)>pi,
        theta(k+1)=-pi;
        w(k+1)=w_max;
        num_spikes=num_spikes+1;
        k_spikes(num_spikes)=k;
    end;
end;

subplot(211);
t=[0:m_steps]*dt;
plot(t,1-cos(theta),'-k','Linewidth',2);
set(gca,'Fontsize',16);
axis([0,t_final,0,2]);
ylabel('$1-\cos \theta$','Fontsize',20); 
subplot(212);
plot(t(1:k_spikes(1)),w(1:k_spikes(1)),'-k','Linewidth',2);
hold on;
t_mid=(t(k_spikes(1))+t(k_spikes(1)+1))/2;
plot([t_mid,t_mid],[w(k_spikes(1)),w_max],':k','Linewidth',2);
for i=1:num_spikes-1,
    plot(t(k_spikes(i)+1:k_spikes(i+1)), ...
        w(k_spikes(i)+1:k_spikes(i+1)),'-k','Linewidth',2);
    t_mid=(t(k_spikes(i+1))+t(k_spikes(i+1)+1))/2;
    plot([t_mid,t_mid],[w(k_spikes(i+1)),w_max],':k','Linewidth',2);
end;
plot(t(k_spikes(num_spikes)+1:length(w)), ...
    w(k_spikes(num_spikes)+1:length(w)),'-k','Linewidth',2);
hold off;
set(gca,'Fontsize',16);
if w_max>0,
    axis([0,t_final,0,w_max]);
end;
xlabel('$t$ [ms]','Fontsize',20);
ylabel('$z$','Fontsize',20); 
shg;
