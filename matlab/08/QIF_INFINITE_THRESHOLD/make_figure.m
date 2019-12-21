clear;
tau_m=2;
I=0.15;

T=2*tau_m/sqrt(tau_m*I-1/4)*atan(1/(2*sqrt(tau_m*I-1/4)));
t_ast=tau_m/sqrt(tau_m*I-1/4)*(pi/2-atan(1/(2*sqrt(tau_m*I-1/4))));

t_period=[0:100]/100*T; t=t_period;
v_0_to_1=0.5+sqrt(tau_m*I-1/4)*tan(sqrt(tau_m*I-1/4)/tau_m*t- ...
    atan(1/(2*sqrt(tau_m*I-1/4))));

t_blowup=[0:99]/100*t_ast; t=t_blowup;
v_1_to_inf=0.5+sqrt(tau_m*I-1/4)*tan(sqrt(tau_m*I-1/4)/tau_m*t+ ...
    atan(1/(2*sqrt(tau_m*I-1/4))));

v_minus_inf_to_0=1-v_1_to_inf(100:-1:1);

subplot(211);

plot(t_period,v_0_to_1,'-k','Linewidth',3);
hold on;
plot(T+t_blowup,v_1_to_inf,'-k','Linewidth',1);
plot([T+t_ast,T+t_ast],[-1,2],'--k','Linewidth',1);

for ijk=1:5,
    plot(ijk*(T+2*t_ast)-t_ast+t_blowup,v_minus_inf_to_0,'-k','Linewidth',1);
    plot(ijk*(T+2*t_ast)+t_period,v_0_to_1,'-k','Linewidth',3);
    plot(ijk*(T+2*t_ast)+T+t_blowup,v_1_to_inf,'-k','Linewidth',1);
    plot([ijk*(T+2*t_ast)+T+t_ast,ijk*(T+2*t_ast)+T+t_ast],[-1,2],'--k','Linewidth',1);
end;

axis([0,150,-1,2]);
hold off;
set(gca,'Fontsize',16);
axis([0,150,-1,2]);
xlabel('$t$','Fontsize',20); ylabel('$v$','Fontsize',20);

hold off;
shg;

