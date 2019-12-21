clear; clf;

tau=10;
num=2;
g_gap=0.01;
beta=3;
epsilon=beta*g_gap;
G=[0 g_gap;
   g_gap 0];
c=sum(G)';


i_ext=[0.125; 0.09];

t_final=100;
dt=0.002;
dt05=dt/2;
m_steps=round(t_final/dt);

z=zeros(num,m_steps+1);
v=z; m=z; h=z; n=z;

v(:,1)=[0.4;0.9];



for k=1:m_steps,
    
    v_inc=-v(:,k)/tau+i_ext+G*v(:,k)-c.*v(:,k);
    v(:,k+1)=v(:,k)+dt*v_inc;
    w=v(:,k+1);
    ind=find(w==max(w)); ind=min(ind);
    if w(ind)>1,
        v(ind,k+1)=0;
        v(3-ind,k+1)=v(3-ind,k+1)+epsilon;
        if v(3-ind,k+1)>1,
            v(3-ind,k+1)=0;
        end;
        tt=k*dt;
        subplot(2,1,ind);
        plot([tt,tt],[0,6],'-k','Linewidth',3);
        hold on;
    end;
end;

t=[0:m_steps]*dt;
set(gca,'Fontsize',16);
plot(t,v(1,:),'-k','Linewidth',3);
ylabel('$v_1$ [mV]','Fontsize',20);
hold on;
axis([0,t_final,0,6]);

subplot(212);
set(gca,'Fontsize',16);
plot(t,v(2,:),'-k','Linewidth',3);
xlabel('$t$ [ms]','Fontsize',20); ylabel('$v_2$ [mV]','Fontsize',20);
hold on;
axis([0,t_final,0.85,0.95]);
shg;

epsilon=0;
v(:,1)=[0.4;0.9];



for k=1:m_steps,
    
    v_inc=-v(:,k)/tau+i_ext+G*v(:,k)-c.*v(:,k);
    v(:,k+1)=v(:,k)+dt*v_inc;
    w=v(:,k+1);
    ind=find(w==max(w)); ind=min(ind);
    if w(ind)>1,
        v(ind,k+1)=0;
        v(3-ind,k+1)=v(3-ind,k+1)+epsilon;
        if v(3-ind,k+1)>1,
            v(3-ind,k+1)=0;
        end;
        tt=k*dt;
        subplot(2,1,ind);
        plot([tt,tt],[0,6],'-r','Linewidth',1);
        hold on;
    end;
end;

subplot(2,1,1);
plot(t,v(1,:),'-r','Linewidth',1);
hold off;

subplot(2,1,2);
plot(t,v(2,:),'-r','Linewidth',1);
xlabel('$t$ [ms]','Fontsize',20); ylabel('$v_2$ [mV]','Fontsize',20);
hold off;

shg;
    
    
