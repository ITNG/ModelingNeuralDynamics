clear; clf;

tic
%a=100/21; tau_n=100/sqrt(3); % parameters in the FN-equation
a=5; tau_n=60;               % parameters in the FN-equation
 % the canard explosion happens in this interval;



t_final=10000;               % we find the amplitude of the limit cycle
                            % by solving from time 0 to time t_final ...
dt=0.01; dt05=dt/2;         % ... with time step dt
m_steps=round(t_final/dt);
v=zeros(m_steps+1,1); 
n=zeros(m_steps+1,1);
subplot(111);
I=-4.256889
v(1)=0; n(1)=0;
for k=1:m_steps,
    v_inc=v(k)-v(k)^3/3-n(k)+I;
    n_inc=(a*v(k)-n(k))/tau_n;
    v_tmp=v(k)+dt05*v_inc;
    n_tmp=n(k)+dt05*n_inc;
    v_inc=v_tmp-v_tmp^3/3-n_tmp+I;
    n_inc=(a*v_tmp-n_tmp)/tau_n;
    v(k+1)=v(k)+dt*v_inc;
    n(k+1)=n(k)+dt*n_inc;
end;
vv=v(4*m_steps/5:m_steps+1); 
nn=n(4*m_steps/5:m_steps+1); 
plot(vv,nn,'-b','Linewidth',5);   
hold on;

set(gca,'Fontsize',24);
xlabel('$v$','Fontsize',32);
ylabel('$n$','Fontsize',32);
A=-2; B=0.5; C=-5; D=-4;
axis([A,B,C,D]); axis('square');

v_plot=A+(0:100)/100*(B-A);
n_plot=v_plot-v_plot.^3/3+I;
plot(v_plot,n_plot,'-g','Linewidth',2);
n_plot=a*v_plot;
plot(v_plot,n_plot,'-r','Linewidth',2);
hold off;


shg;


