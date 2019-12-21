clear; clf; rng('default'); rng(63806);

global alpha Period g_bar m;

alpha=1;
Period=25;
g_bar=0.1;
N=2000; m=mean(exp(alpha*cos(pi*[0:N-1]/N).^2)-1);
tau=10;
I_vec=[1:100]/100*0.2+0.1;

tic
t_final=2000;
dt=0.01; dt05=dt/2;
m_steps=round(t_final/dt);
t=[0:m_steps]'*dt;

v=zeros(m_steps+1,1);

for ij=1:length(I_vec),
    I=I_vec(ij);
    v(1)=0;
    k_old=1;
    num_spikes=0;
    for k=1:m_steps,
        v_inc=-v(k)/tau+I-g((k-1)*dt)*v(k);
        v_tmp=v(k)+dt05*v_inc;
        v_inc=-v_tmp/tau+I-g((k-1/2)*dt)*v_tmp;
        v(k+1)=v(k)+dt*v_inc;
        if v(k+1)>1,
           
            k_old=k+1;
            v(k+1)=0;
            num_spikes=num_spikes+1;
        end;       
    end;
    f_vec(ij)=num_spikes/t_final*1000;
end;

subplot(211);
L=length(I_vec);
for k=1:L-1,
    if f_vec(k+1)-f_vec(k)<=1,
        plot([I_vec(k),I_vec(k+1)],[f_vec(k),f_vec(k+1)],'-b','Linewidth',2);
        hold on;
    else
        plot([(I_vec(k)+I_vec(k+1))/2,(I_vec(k)+I_vec(k+1))/2], ...
            [f_vec(k),f_vec(k+1)],'--b','Linewidth',1);
        hold on;
    end;
end;


f=[1:200];
% Solve for I in terms of f:

A=f/1000/(1+tau*g_bar);
A=1./A;
A=A/tau;
A=exp(A);
I=A*(1+tau*g_bar)./tau./(A-1);
plot(I,f,'-r','Linewidth',2);
axis([min(I_vec),max(I_vec),0,200]);


set(gca,'Fontsize',16);
xlabel('$I$','Fontsize',20);
ylabel('$f$ [Hz]','Fontsize',20);

plot([min(I_vec),1/tau+g_bar],[0,0],'-r','Linewidth',4);

I_onset=1/tau+g_bar;
plot(I_onset,0,'.r','Markersize',40);
for k=1:L-1,
    if f_vec(k)==0 & f_vec(k+1)>0, 
        I_onset=(I_vec(k)+I_vec(k+1))/2;
        plot(I_onset,0,'.b','Markersize',40);
    end;
end;
hold off;





shg;
toc
