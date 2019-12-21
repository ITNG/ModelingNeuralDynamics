function F=F(v);

global c g_k g_na g_l v_k v_na v_l i_ext

F=g_na*(m_inf(v)).^3.*(1-n_inf(v)).*(v_na-v)+g_k*n_inf(v).^4.*(v_k-v)+ ...
    g_l*(v_l-v)+i_ext;