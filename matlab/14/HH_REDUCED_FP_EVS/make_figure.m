clear; clf;
c=1;
g_k=36; 
g_na=120;
g_l=0.3;
v_k=-82;
v_na=45;
v_l=-59;

f=@(v) g_na*m_inf(v).^3.*(0.83-n_inf(v)).*(v_na-v)+ ...
       g_k*n_inf(v).^4.*(v_k-v)+ ...
       g_l*(v_l-v);

figure(3);
subplot(211);
i_ext_vec=[0:1000]/1000*15;
for ijk=1:length(i_ext_vec),
    i_ext=i_ext_vec(ijk);
    v_left=-100;
    v_right=50;
    while v_right-v_left>10^(-10),
        v_c=(v_left+v_right)/2;
        if (f(v_c)+i_ext)*(f(v_left)+i_ext)>0,
            v_left=v_c;
        else
            v_right=v_c;
        end;
    end;
    v_c=(v_left+v_right)/2;
    fp_vec(ijk)=v_c;
    n_c=n_inf(v_c);
    J(1,1,ijk)=g_na*3*m_inf(v_c)^2*m_inf_p(v_c)*...
        (0.83-n_c)*(v_na-v_c)- ...
        g_na*m_inf(v_c)^3*(0.83-n_c)- ...
        g_k*n_c^4-g_l;
    J(1,2,ijk)=-g_na*m_inf(v_c)^3*(v_na-v_c)+...
        4*g_k*n_c^3*(v_k-v_c);
    J(2,1,ijk)=alpha_n_p(v_c)*(1-n_c)-beta_n_p(v_c)*n_c;
    J(2,2,ijk)=-alpha_n(v_c)-beta_n(v_c);
    E=eig(J(:,:,ijk));
    if real(E(1))<0 & real(E(2))<0 & abs(imag(E(1)))>10^(-4),
        plot(i_ext,v_c,'.b','Markersize',10);
        typ(ijk)=1;
    end;
    hold on;
    if real(E(1))>0 & real(E(2))>0 & abs(imag(E(1)))>10^(-4),
        plot(i_ext,v_c,'.r','Markersize',10);
        typ(ijk)=2;
    end;
    imag_part(ijk)=imag(E(1));
    real_part(ijk)=real(E(1));
    v_c_vec(ijk)=v_c;
    hold on;
end;
hold off;
set(gca,'Fontsize',16);
xlabel('$I_{\rm ext}$','Fontsize',20);
ylabel('$v_\ast$','Fontsize',20);
shg;

figure(2);
subplot(211);
ind=find(typ==1);
plot(i_ext_vec(ind),v_c_vec(ind),'-k','Linewidth',2);
hold on;
ind=find(typ==2);
plot(i_ext_vec(ind),v_c_vec(ind),'--k','Linewidth',2);
hold off;
set(gca,'Fontsize',16);
xlabel('$I_{\rm ext}$','Fontsize',20);
ylabel('$v_\ast$','Fontsize',20);

figure(1);
subplot(111);
A=-1; B=1; C=-1; D=1;
plot(real_part,imag_part,'-k','Linewidth',2);
hold on;
plot(real_part,-imag_part,'-k','Linewidth',2);
N=length(real_part);
ind=find(real_part(1:N-1)<0.1&real_part(2:N)>0.1);
ind=ind(1);
x=real_part(ind);
y=imag_part(ind);
v=zeros(2,1);
v(1)=real_part(ind)-real_part(ind-200);
v(2)=imag_part(ind)-imag_part(ind-200);
epsilon=0.1;
width=2;
col='-k';
arrow(A,B,C,D,x,y,v,epsilon,width,col);
y=-y;
v(2)=-v(2);
arrow(A,B,C,D,x,y,v,epsilon,width,col);
hold off
set(gca,'Fontsize',24);
xlabel('${\rm Re}~\!(\lambda)$','Fontsize',32);
ylabel('${\rm Im}~\!(\lambda)$','Fontsize',32);
axis([A,B,C,D]);
axis('square');
hold on;
plot([0,0],[-1,1],'-b','Linewidth',1)
hold off;




shg;
    
