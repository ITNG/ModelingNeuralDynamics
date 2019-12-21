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

subplot(211);
i_ext_vec=[0:1000]/1000*15;
green_i=[]; green_v=[];
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
        plot(i_ext,v_c,'.r','Markersize',10);
    end;
    hold on;
    if real(E(1))>0 & real(E(2))>0 & abs(imag(E(1)))>10^(-4),
        green_v=[green_v,v_c];
        green_i=[green_i,i_ext];
    end;
    imag_part(ijk)=imag(E(1));
    real_part(ijk)=real(E(1));
    v_c_vec(ijk)=v_c;
    hold on;
end;
plot(green_i,green_v,'--g','Linewidth',3)
hold off;
set(gca,'Fontsize',16);
xlabel('$I$','Fontsize',20);
ylabel('$v_\ast$','Fontsize',20);
shg;
