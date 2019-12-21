clear; clf;

global c g_k g_na g_l v_k v_na v_l i_ext

c=1;
g_k=80; 
g_na=100;
g_l=0.1;
v_k=-100;
v_na=50;
v_l=-67;

N_i=1000; N_v=1000;
i_ext_vec=[0:N_i]/N_i*0.2;


subplot(111);


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
        while B-A>10^(-10),
            C=(A+B)/2;
            if F(C)*F(A)<=0,
                B=C;
            else
                A=C;
            end;
        end;
        v=(A+B)/2;
        n=n_inf(v);

        J(1,1)=g_na*3*m_inf(v)^2*m_inf_p(v)*(1-n)*(v_na-v) ...
            -g_na*m_inf(v)^3*(1-n) ...
            -g_k*n^4-g_l;
        J(1,2)=-g_na*m_inf(v)^3*(v_na-v)+g_k*4*n^3*(v_k-v);
        J(2,1)=alpha_n_p(v)*(1-n)-beta_n_p(v)*n;
        J(2,2)=-alpha_n(v)-beta_n(v);
        E=eig(J);
        if abs(imag(E(1)))>10^(-4),
            if real(E(1))>0,
                plot(i_ext,v,'.g','Markersize',10);
                hold on;
            end;
            if real(E(1))<0,
                plot(i_ext,v,'.r','Markersize',10);
                hold on;
            end;
        end;
            
        if E(1)<0 & E(2)<0,
            plot(i_ext,v,'.k','Markersize',10);
            hold on;
        end;
        if E(1)*E(2)<0,
            plot(i_ext,v,'.m','Markersize',10);
            hold on;
        end;
        if E(1)>0 & E(2)>0,
            plot(i_ext,v,'.b','Markersize',10);
            hold on;
        end;
       
            
    end;
   
end;


hold off;
set(gca,'Fontsize',24);
xlabel('$I$ [$\mu$A/cm$^2$]','Fontsize',32);
ylabel('$v_\ast$ [mV]','Fontsize',32);
axis([min(i_ext_vec),max(i_ext_vec),-70,-30]);

shg;






    
    
