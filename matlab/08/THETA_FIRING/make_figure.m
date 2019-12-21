tau_m=1/2;
I=0.505;



t_final=150;

dt=0.001; dt05=dt/2;

m_steps=round(t_final/dt);

theta=zeros(m_steps+1,1);

for k=1:m_steps,
    theta_inc=-cos(theta(k))/tau_m+2*I*(1+cos(theta(k)));
    theta_tmp=theta(k)+dt05*theta_inc;
    theta_inc=-cos(theta_tmp)/tau_m+2*I*(1+cos(theta_tmp));
    theta(k+1)=theta(k)+dt*theta_inc;
end;

subplot(211);
plot([0:m_steps]*dt,1-cos(theta),'-k','Linewidth',2);
set(gca,'Fontsize',16);
ylabel('$1-cos(\theta)$','Fontsize',20); xlabel('$t$','Fontsize',20);

shg;

