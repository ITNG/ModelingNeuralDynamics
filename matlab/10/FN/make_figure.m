clear;
v=[-100:100]/100*3;
a=1.25; tau_n=15.625; I_ext=-0.5;


subplot(221);
plot(v,a*v,'-k','Linewidth',2);
set(gca,'Fontsize',16);
xlabel('$v$','Fontsize',20); ylabel('$n$','Fontsize',16);
axis([-3,3,-3,3]); 
hold on;
	

plot(v,v-v.^3/3+I_ext,'-r','Linewidth',2);

t_final=400;
dt=0.01; dt05=dt/2;
m_steps=round(t_final/dt);
v=zeros(m_steps+1,1); n=v;

v(1)=-1;
n(1)=-2;


for k=1:m_steps,
    v_inc=v(k)-v(k)^3/3-n(k)+I_ext;
    n_inc=(a*v(k)-n(k))/tau_n;
    v_tmp=v(k)+dt05*v_inc;
    n_tmp=n(k)+dt05*n_inc;
    v_inc=v_tmp-v_tmp^3/3-n_tmp+I_ext;
    n_inc=(a*v_tmp-n_tmp)/tau_n;
    v(k+1)=v(k)+dt*v_inc;
    n(k+1)=n(k)+dt*n_inc;
end;

plot(v,n,'-b','Linewidth',2);
hold off;

hold off;

subplot(222);
plot([0:m_steps]*dt,v,'-k','Linewidth',2);
set(gca,'Fontsize',16);
xlabel('$t$','Fontsize',16); ylabel('$v$','Fontsize',16);
axis([0,t_final,-3,3]);

shg;


