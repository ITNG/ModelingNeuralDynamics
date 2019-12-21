function alpha_n=alpha_n(v)
ind_one=find(abs(v-95)>0.1);
ind_two=find(abs(v-95)<=0.1);
one=(95-v)./(exp((-95+v)/-11.8)-1);
two=11.8+0.4999999999*(-95+v)+1033078/146283845*(-95+v).^2-3.54*10^-13*(-95+v).^3-80931/9574*10^-7*(-95+v).^4;
alpha_n(ind_one)=one(ind_one);
alpha_n(ind_two)=two(ind_two);