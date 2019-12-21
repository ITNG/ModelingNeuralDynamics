function alpha_m=alpha_m(v);
if abs(v+45)>1.0e-8,
    alpha_m=(v+45)/10./(1-exp(-(v+45)/10));
else
    alpha_m=1;
end;