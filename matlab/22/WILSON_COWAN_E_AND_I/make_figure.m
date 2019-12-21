clear; clf;

t_final=300;
I_E=20; I_I=0;
w_EE=1.5; w_IE=1;
w_EI=1; w_II=0;
tau_E=5; tau_I=10; 

dt=0.01; dt05=dt/2; m_steps=round(t_final/dt);

z=zeros(m_steps+1,1); E=z; I=z;
E(1)=50; I(1)=10;


for k=1:m_steps,
    E_inc=(f(w_EE*E(k)-w_IE*I(k)+I_E)-E(k))/tau_E;
    I_inc=(g(w_EI*E(k)-w_II*I(k)+I_I)-I(k))/tau_I;
    E_tmp=E(k)+dt05*E_inc;
    I_tmp=I(k)+dt05*I_inc;
    E_inc=(f(w_EE*E_tmp-w_IE*I_tmp+I_E)-E_tmp)/tau_E;
    I_inc=(g(w_EI*E_tmp-w_II*I_tmp+I_I)-I_tmp)/tau_I;
    E(k+1)=E(k)+dt*E_inc;
    I(k+1)=I(k)+dt*I_inc;
end;

t=[0:m_steps]*dt;
subplot(211);
plot(t,E,'-r','Linewidth',2)
hold on;
plot(t,I,'-b','Linewidth',2);
set(gca,'Fontsize',16);
xlabel('$t$ [ms]','Fontsize',20);
title('$E$ (red) and $I$ (blue)','Fontsize',20);
axis([0,t_final,0,100]);
hold off;

shg;

