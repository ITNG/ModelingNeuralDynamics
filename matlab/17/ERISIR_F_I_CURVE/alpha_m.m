function alpha_m=alpha_m(v) 
one=20*(2*v-151)./(1-exp((151-2*v)/27));
two=540+20*(v-151/2)+20/81*(v-75.5).^2-4/177147*(v-75.5).^4;
ind_one=find(abs(v-75.5)>0.1);
ind_two=find(abs(v-75.5)<=0.1);
alpha_m(ind_one)=one(ind_one);
alpha_m(ind_two)=two(ind_two);