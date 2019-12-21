clear; clf; 

T=20;
tau=30;
n=10;
r=exp(-T/tau);
% left=(1-r)/(1-r^n); 
% right=(1-r)/(1-r^(n-1));
% epsilon=(left+right)/2;
epsilon=0.4875;
t=[0:100]/100*T;
t_final=1000;

subplot(211);

plot([0,T],[0,0],'-k','Linewidth',4);
hold on;
plot([T,T],[0,epsilon],':k','Linewidth',1);
v_post=epsilon;

for k=1:round(t_final/T),
    plot(k*T+t,v_post*exp(-t/tau),'-k','Linewidth',2);
    v_pre=v_post*r;
    v_post=v_pre+epsilon;
    v_post=v_post*(v_post<1);
    if v_post==0,
        plot([(k+1)*T,(k+1)*T],[0,5],'-k','Linewidth',2);
    end;
    plot([(k+1)*T,(k+1)*T],[v_pre,v_post],':k','Linewidth',1);
end;
hold off;
set(gca,'Fontsize',16);
xlabel('$t$ [ms]','Fontsize',20);
ylabel('$v$','Fontsize',20);

axis([0,t_final,0,6]);

shg;

