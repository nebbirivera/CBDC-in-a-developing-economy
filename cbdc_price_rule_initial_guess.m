function [c_0 k_0 m_dc_0] = cbdc_price_rule_initial_guess(fixed, point)

	A = fixed(1);
	B = fixed(2);
	tau_c_t = fixed(3);
	delta = fixed(4);
	theta = fixed(5);
	i_t = fixed(6);
	i_dc_t = fixed(7);
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

	% Approximate around
	
	a = point .* [1; (1/0.075); 0.05*sqrt(A/B)];

	[aa bb] = cbdc_price_rule_get_coeffs(fixed,a);

	options = optimset('Display', 'off', 'TolX', 1e-8, 'TolFun', 1e-8,...
		'MaxFunEvals',1e4,'MaxIter',1e4);

	x0 = a;
	[x1 Fval] = fsolve(@(x) bb + aa*[x(1)-a(1);x(2)-a(2);x(3)-a(3)], x0, options);

	c_0 = x1(1);
	k_0 = x1(2);
	m_dc_0 = x1(3);

	if all([c_0 k_0 m_dc_0]<=0)
		display("WARNING: One initial guess value is non-positive.")
		return
	end