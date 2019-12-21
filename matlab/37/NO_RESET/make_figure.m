clear; clf;

tau_m=10;
I=0.12;


T=25;
epsilon=10;
g_bar=1/7;
tau_m_hat=1/(1/tau_m+g_bar*T/epsilon);

v0=1.0;

subplot(211);

for ijk=1:5;
    t=[0:100]/100*epsilon;
    v=v0*exp(-t/tau_m_hat)+tau_m_hat*I*(1-exp(-t/tau_m_hat));
    plot(t+(ijk-1)*T,v,'-k','Linewidth',2);
    hold on;
    v0=v(101);
    t=[0:100]/100*(T-epsilon);
    v=v0*exp(-t/tau_m)+tau_m*I*(1-exp(-t/tau_m));
    plot(t+(ijk-1)*T+epsilon,v,'-k','Linewidth',2);
    v0=v(101);
end;
set(gca,'Fontsize',16);
xlabel('$t$','Fontsize',20);
ylabel('$v$','Fontsize',20);
axis([0,5*T,0,1]);
hold off;
shg;