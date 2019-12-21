clear; clf;

c=1;
g_na=20;
g_k=10; 
g_l=8;
v_na=60;
v_k=-90;
v_l=-80;
tau_n=0.15;

f=@(v) g_na*m_inf(v).*(v_na-v)+g_k*n_inf(v).*(v_k-v)+g_l*(v_l-v);

        % f is the function of which we have to determine the zeros
        % to find the fixed points.

i_ext_vec=-4+(0:1000)/1000*12; % external drives to be considered

v=-100+(0:3000)/3000*150;    % For a given external drive I, we check
                             % the sign of the function f in these 
                             % values of v. So if two zeros lie closer 
                             % than 150/3000 mV to each other, we may miss
                             % both of them.

maxv_c=-100;    % maxv_c and minv_c will be the maximum and minimum
minv_c=50;      % voltages in all fixed points found. (This is useful 
                % to keep track of for the "axis" statement at the end.)
subplot(111);
i_green=[]; v_green=[];
i_magenta=[]; v_magenta=[];
for ijk=1:length(i_ext_vec),
    I=i_ext_vec(ijk);
    k=(1:length(v)-1);
    ind=find((f(v(k))+I).*(f(v(k+1))+I)<=0);
    for klm=1:length(ind),
        j=ind(klm);
        v_low=v(j); v_high=v(j+1);
        while v_high-v_low>10^(-12),
            v_c=(v_low+v_high)/2;
            if (f(v_c)+I)*(f(v_high)+I)<=0,
                v_low=v_c;
            else
                v_high=v_c;
            end;
        end;
        v_c=(v_low+v_high)/2;
        maxv_c=max(maxv_c,v_c);
        minv_c=min(minv_c,v_c);
        n_c=n_inf(v_c);
        J=zeros(2,2);
        J(1,1)=g_na*m_inf_p(v_c)*(v_na-v_c)-g_na*m_inf(v_c)-g_k*n_c-g_l;
        J(1,2)=g_k*(v_k-v_c);
        J(2,1)=n_inf_p(v_c)/tau_n;
        J(2,2)=-1/tau_n;
        E=eig(J);
        if abs(imag(E(1)))<10^(-12) & real(E(1))<0 & real(E(2))<0,
            
            % we test abs(imag(E(1))<10^(-12), not
            % abs(imag(E(1))==0, because the imaginary part of 
            % the computed eigenvalue may be non-zero due to rounding
            % even when the imaginary part is in reality equal to zero.
            
            % stable node
            plot(I,v_c,'.k','Markersize',15);
            hold on;
        end;
        if abs(imag(E(1)))<10^(-12) & real(E(1))>0 & real(E(2))>0,
            % unstable node
            plot(I,v_c,'.b','Markersize',15);
            hold on;
        end;
        if abs(imag(E(1)))>10^(-12) & real(E(1))<0,
            % stable spiral
            plot(I,v_c,'.r','Markersize',15);
            hold on;
        end;
        if abs(imag(E(1)))>10^(-12) & real(E(1))>0,
            % unstable spiral
            i_green=[i_green,I];
            v_green=[v_green,v_c];
        end;
        if abs(imag(E(1)))<10^(-12) & real(E(1))*real(E(2))<0,
            % saddle
            i_magenta=[i_magenta,I];
            v_magenta=[v_magenta,v_c];
        end;
        
    end;
end;
plot(i_green,v_green,'--g','Linewidth',4);
plot(i_magenta,v_magenta,'--m','Linewidth',4);
hold off;
set(gca,'Fontsize',24);
xlabel('$I$','Fontsize',32);
ylabel('$v_\ast$','Fontsize',32);
axis([-4,8,minv_c-5,maxv_c+5]); 

shg;
    