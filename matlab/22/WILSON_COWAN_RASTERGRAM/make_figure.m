clear; clf; rng('default'); rng(63806);


w_EE=1.5; w_IE=1;
w_EI=1; w_II=0;
tau_E=5; tau_I=10; 
I_E=20; I_I=0;

t_final=300;

dt=0.01; dt05=dt/2; m_steps=round(t_final/dt);

z=zeros(m_steps+1,1); E=z; I=z; E(1)=50; I(1)=10;

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



t=[0:m_steps]*dt; t_left=t(1:m_steps); t_right=t(2:m_steps+1);
num_e=80; num_i=20;
U=rand(num_e+num_i,m_steps);

subplot(211);
for k=1:m_steps,
    ind=find(U(num_i+1:num_i+num_e,k)<=dt*(E(k)+E(k+1))/2/1000);
    plot((k-1/2)*dt*ones(length(ind),1),ind+num_i,'.r','Markersize',6);
    hold on;
end;
for k=1:m_steps,
    ind=find(U(1:num_i,k)<=dt*(I(k)+I(k+1))/2/1000);
    plot((k-1/2)*dt*ones(length(ind),1),ind,'.b','Markersize',6);
    hold on;
end;
set(gca,'Fontsize',16);
xlabel('$t$ [ms]','Fontsize',20);
ylabel('neuronal index','Fontsize',20);
axis([0,t_final,0,num_e+num_i]);
hold off;

box on;

shg;

