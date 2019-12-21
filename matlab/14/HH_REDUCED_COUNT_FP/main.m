clear; clf;
c=1;
g_k=36; 
g_na=120;
g_l=0.3;
v_k=-82;
v_na=45;
v_l=-59;

v=-100+[0:10000]/10000*150;
n=length(v);
i_left=[1:n-1];
i_right=[2:n];


i_ext_vec=[0:1000]/1000*15;

for ijk=1:length(i_ext_vec),
    i_ext=i_ext_vec(ijk);
    f=g_na*m_inf(v).^3.*(0.83-n_inf(v)).*(v_na-v)+ ...
       g_k*n_inf(v).^4.*(v_k-v)+ ...
       g_l*(v_l-v)+i_ext;
    ind=find(f(i_left).*f(i_right)<=0);
    num_fp(ijk)=length(ind);
end;
max(num_fp)
min(num_fp)
    
