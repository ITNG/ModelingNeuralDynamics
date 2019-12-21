function psif=psif(tau_m,I,g_I,tau_I,dt,x);

% x can be a column vector, in which case so is psi.

v_1=0;
v_2=x;
s=1;
N=length(x);
psif=zeros(N,1);
done=zeros(N,1);

dt05=dt/2;

while min(done)==0,
    v_1_old=v_1;
    v_2_old=v_2;
    v_1_inc=-v_1/tau_m+I-g_I*s*v_1;
    v_2_inc=-v_2/tau_m+I-g_I*s*v_2;
    v_1_tmp=v_1+dt05*v_1_inc;
    v_2_tmp=v_2+dt05*v_2_inc;
    s_tmp=s*exp(-dt05/tau_I);
    v_1_inc=-v_1_tmp/tau_m+I-g_I*s_tmp*v_1_tmp;
    v_2_inc=-v_2_tmp/tau_m+I-g_I*s_tmp*v_2_tmp;
    v_1=v_1+dt*v_1_inc;
    v_2=v_2+dt*v_2_inc;
    ind=find(v_2>1 & done==0);
    psif(ind)=(v_1_old*(v_2(ind)-1)+v_1*(1-v_2_old(ind)))./ ...
        (v_2(ind)-v_2_old(ind));
    done(ind)=1;
    s=s*exp(-dt/tau_I);
end;
    