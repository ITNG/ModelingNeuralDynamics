clear; clf;

global c g_k g_na g_l v_k v_na v_l i_ext

c=1;
g_k=224; 
g_na=112;
g_l=0.5;
v_k=-90;
v_na=60;
v_l=-70;

N_i=1000; N_v=1000;
i_ext_vec=[0:N_i]/N_i*7;

subplot(111);
i_magenta=[];
i_blue=[];
v_magenta=[];
v_blue=[];
for ijk=1:length(i_ext_vec),
    i_ext=i_ext_vec(ijk);
    v_min=min(v_k,v_l+i_ext/g_l);
    v_max=max(v_na,v_l+i_ext/g_l);

    v_vec=v_min+[0:N_v]/N_v*(v_max-v_min);

    v_vec_left=v_vec(1:N_v); v_vec_right=v_vec(2:N_v+1);
    ind=find(F(v_vec_left).*F(v_vec_right)<=0);
    
    
    for i=1:length(ind),
        A=v_vec_left(ind(i));
        B=v_vec_right(ind(i));
        while B-A>10^(-12),
            C=(A+B)/2;
            if F(C)*F(A)<=0,
                B=C;
            else
                A=C;
            end;
        end;
        v=(A+B)/2;
        n=n_inf(v);

        J(1,1)=g_na*3*m_inf(v)^2*m_inf_p(v)*(0.36-n)*(v_na-v) ...
            -g_na*m_inf(v)^3*(0.36-n) ...
            -g_k*n^2-g_l;
        J(1,2)=-g_na*m_inf(v)^3*(v_na-v)+g_k*2*n*(v_k-v);
        J(2,1)=alpha_n_p(v)*(1-n)-beta_n_p(v)*n;
        J(2,2)=-alpha_n(v)-beta_n(v);
        E=eig(J);
        if abs(imag(E(1)))>10^(-6),
            if real(E(1))>0,
                plot(i_ext,v,'.g','Markersize',10);
                hold on;
            end;
            if real(E(1))<0,
                plot(i_ext,v,'.r','Markersize',10);
                hold on;
            end;
        end;
        
        if abs(imag(E(1)))<10^(-6),
            if E(1)<0 & E(2)<0,
                plot(i_ext,v,'.k','Markersize',10);
                hold on;
            end;
            if E(1)*E(2)<=0,
                i_magenta=[i_magenta,i_ext];
                v_magenta=[v_magenta,v];
            end;
            if E(1)>0 & E(2)>0,
                i_blue=[i_blue,i_ext];
                v_blue=[v_blue,v];
            end;
        end;
  
    end;
   
end;
hold on;
plot(i_magenta,v_magenta,'--m','Linewidth',2);
plot(i_blue,v_blue,'--b','Linewidth',2);


hold off;
set(gca,'Fontsize',24);
xlabel('$I$ [$\mu$A/cm$^2$]','Fontsize',32);
ylabel('$v_\ast$ [mV]','Fontsize',32);
axis([min(i_ext_vec),max(i_ext_vec),-90,0]);


shg;






    
    
