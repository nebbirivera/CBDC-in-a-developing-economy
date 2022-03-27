function F = num_eqs(x, fixed)%, price_rule)

	A = fixed(1);
	B = fixed(2);
	tau_c_t = fixed(3);
	delta = fixed(4);
	theta = fixed(5);
	i_t = fixed(6);
	% if price_rule
	i_dc_t = fixed(7);
	% else
	
	% end
	ups_p_t = fixed(8);
	betta = fixed(9);
	alpha_k = fixed(10);
	mc_t = fixed(11);
	a_t = fixed(12);
	y_L_t = fixed(13);
	tau_r_t = fixed(14);
	un_t = fixed(15);
	g_t = fixed(16);
	xi_F = fixed(17);
	eta_F = fixed(18);
	xi_I = fixed(19);
	eta_I = fixed(20);
	omega = fixed(21);
	T = fixed(22);
	L_F_t = fixed(23);
	L_I_t = fixed(24);

	if all(x>=0)
		c_t = x(1);
		k_t = x(2);
		m_dc_t = x(3);
	else
		F = 1e10*abs(rand(3,1));
		return
	end

	% LGF = (m_dc_t^(theta)*(T^(theta)+((i_t-i_dc_t)/(T^(theta)*i_t))^(theta/(1-theta))))
	% LGF = (m_c_t^(theta)*(1+(T*i_t/(i_t-i_dc_t))^(theta/(1-theta))))
	% F(1) = (c_t*(1+(1+tau_c_t)*(A*((1+tau_c_t)*c_t/(m_c_t^(theta)*(1+(T*i_t/(i_t-i_dc_t))^(theta/(1-theta)))))+B*((m_c_t^(theta)*(1+(T*i_t/(i_t-i_dc_t))^(theta/(1-theta))))/((1+tau_c_t)*c_t))-2*sqrt(A*B))) + (delta*k_t)*(1+A*((delta*k_t)/(m_c_t^(theta)*(1 + (T*i_t/((i_t-i_dc_t)))^(theta/(1-theta))))) + B*((m_c_t^(theta)*(1 + (T*i_t/((i_t-i_dc_t)))^(theta/(1-theta))))/(delta*k_t)) - 2*sqrt(A*B)) + un_t^(omega/(omega-1))*(xi_F*(L_F_t*eta_F)^(1/(1-omega))+xi_I*(L_I_t*eta_I)^(1/(1-omega)))) - (a_t/ups_p_t)*(k_t)^(alpha_k)*(y_L_t)^(1-alpha_k)*(1-g_t);
	F(1) = ((c_t + (delta*k_t))/(1-g_t) + c_t*((1+tau_c_t)*(A*((1+tau_c_t)*c_t/(m_dc_t^(theta)*(T^(theta)+((i_t-i_dc_t)/(T^(theta)*i_t))^(theta/(1-theta)))))+B*((m_dc_t^(theta)*(T^(theta)+((i_t-i_dc_t)/(T^(theta)*i_t))^(theta/(1-theta))))/((1+tau_c_t)*c_t))-2*sqrt(A*B))) + (delta*k_t)*(A*((delta*k_t)/(m_dc_t^(theta)*(T^(theta)+((i_t-i_dc_t)/(T^(theta)*i_t))^(theta/(1-theta))))) + B*((m_dc_t^(theta)*(T^(theta)+((i_t-i_dc_t)/(T^(theta)*i_t))^(theta/(1-theta))))/(delta*k_t)) - 2*sqrt(A*B)) + un_t^(omega/(omega-1))*(xi_F*(L_F_t*eta_F)^(1/(1-omega))+xi_I*(L_I_t*eta_I)^(1/(1-omega)))) - (a_t/ups_p_t)*(k_t)^(alpha_k)*(y_L_t)^(1-alpha_k);
	F(2) = (A*(((1+tau_c_t)*c_t/(m_dc_t^(theta)*(T^(theta)+((i_t-i_dc_t)/(T^(theta)*i_t))^(theta/(1-theta)))))^2+((delta*k_t)/(m_dc_t^(theta)*(T^(theta)+((i_t-i_dc_t)/(T^(theta)*i_t))^(theta/(1-theta)))))^2) - 2*B)*theta*(m_dc_t*((i_t-i_dc_t)/(T^(theta)*i_t))^(1/(1-theta)))^(theta-1) - i_t/(1+i_t);
	F(3) = k_t - ((a_t*alpha_k*mc_t*y_L_t^(1-alpha_k))/(ups_p_t*((1+2*(A*((delta*k_t)/(m_dc_t^(theta)*(T^(theta)+((i_t-i_dc_t)/(T^(theta)*i_t))^(theta/(1-theta)))))-sqrt(A*B)))*((1/betta)-1+delta)/(1-tau_r_t))))^(1/(1-alpha_k));

	% F

	% norm(F)