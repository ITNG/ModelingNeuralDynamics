
clear; clf;

I_E=20; I_I=0;
w_EE=1.5; w_IE=1;
w_EI=1; w_II=0;
tau_E=5; tau_I=10; 
t_final=300;

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

subplot(111);
t=[0:m_steps]*dt;
plot(E,I,'-k','Linewidth',4);
A=0; B=100; C=0; D=100; axis([A,B,C,D]); axis('square');
set(gca,'Fontsize',24);
xlabel('$E$','Fontsize',32); ylabel('$I$','Fontsize',32);
hold on;
epsilon=0.05; width=3; col='-k';
k=round(m_steps/4);
arrow(A,B,C,D,E(k),I(k),[E(k)-E(k-1); I(k)-I(k-1)],epsilon,width,col)

N=1000;
I_left=[0:N-1]/N*100;
I_right=[1:N]/N*100;
i_red=0;
for i=0:N,
    E=i/N*100;
    R_left=f(w_EE*E-w_IE*I_left+I_E)-E;
    R_right=f(w_EE*E-w_IE*I_right+I_E)-E;
    ind=find(R_left.*R_right<0);
    if length(ind)>0,
    for j=1:length(ind);
        I_l=I_left(ind(j)); 
        I_r=I_right(ind(j));
        while I_r-I_l>10^(-7),
            I_c=(I_r+I_l)/2;
            R_l=f(w_EE*E-w_IE*I_l+I_E)-E;
            R_c=f(w_EE*E-w_IE*I_c+I_E)-E;
        
            if R_l*R_c<0,
                I_r=I_c;
            else
                I_l=I_c;
            end;
        end;
        I_c=(I_r+I_l)/2;
        i_red=i_red+1;
        E_red(i_red)=E; I_red(i_red)=I_c;
    end;
    end;
end;

i_blue=0;
for i=0:N,
    E=i/N*100;
    R_left=g(w_EI*E-w_II*I_left+I_I)-I_left;
    R_right=g(w_EI*E-w_II*I_right+I_I)-I_right;
    ind=find(R_left.*R_right<0);
    if length(ind)>0,
    for j=1:length(ind);
        I_l=I_left(ind(j)); I_r=I_right(ind(j));
        while I_r-I_l>10^(-7),
            I_c=(I_r+I_l)/2;
            R_l=g(w_EI*E-w_II*I_l+I_I)-I_l;
            R_c=g(w_EI*E-w_II*I_c+I_I)-I_l;
            if R_l*R_c<0,
                I_r=I_c;
            else
                I_l=I_c;
            end;
        end;
        I_c=(I_r+I_l)/2;
        i_blue=i_blue+1;
        E_blue(i_blue)=E; I_blue(i_blue)=I_c;
    end;
    end;
end;
plot(E_red,I_red,'-r','Linewidth',2);
plot(E_blue,I_blue,'-b','Linewidth',2);
hold off;
    

shg;
