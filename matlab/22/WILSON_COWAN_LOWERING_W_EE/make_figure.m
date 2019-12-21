clear; clf;

w_EE=1.5; w_IE=1;
w_EI=1; w_II=0;
tau_E=5; tau_I=10; 
I_E=20; I_I=0;

t_final=1000;

dt=0.01; dt05=dt/2; m_steps=round(t_final/dt);

z=zeros(m_steps+1,1); E=z; I=z; E(1)=10; I(1)=10;

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

subplot(331);
A=0; B=100; C=0; D=100;
ind=round(4*m_steps/5):m_steps+1;
E=E(ind); I=I(ind); 
plot(E,I,'-k','Linewidth',2);
axis([A,B,C,D]);
axis('square');
set(gca,'Fontsize',12);
xlabel('$E$','Fontsize',16);
ylabel('$I$','Fontsize',16);
str_w_EE=num2str(w_EE);
title(['$w_{EE}=$', str_w_EE],'Fontsize',16);

w_EE=1.25; 

E(1)=10; I(1)=10;

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

subplot(332);
A=0; B=100; C=0; D=100;
E=E(ind); I=I(ind); 
plot(E,I,'-k','Linewidth',2);
axis([A,B,C,D]);
axis('square');
set(gca,'Fontsize',12);
xlabel('$E$','Fontsize',16);
str_w_EE=num2str(w_EE);
title(['$w_{EE}=$', str_w_EE],'Fontsize',16);


w_EE=1.0; 

E(1)=10; I(1)=10;

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

subplot(333);
A=0; B=100; C=0; D=100;
E=E(ind); I=I(ind); 
plot(E,I,'-k','Linewidth',2);
axis([A,B,C,D]);
axis('square');
set(gca,'Fontsize',12);
xlabel('$E$','Fontsize',16);
str_w_EE=num2str(w_EE);
title(['$w_{EE}=$', str_w_EE],'Fontsize',16);




shg;