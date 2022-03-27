function [c_0 k_0 m_c_0] = nocbdc_initial_guess(fixed, point)

	A = fixed(1);
	B = fixed(2);
	tau_c_t = fixed(3);
	delta = fixed(4);
	theta = fixed(5);
	i_t = fixed(6);
	ups_p_t = fixed(7);
	betta = fixed(8);
	alpha_k = fixed(9);
	mc_t = fixed(10);
	a_t = fixed(11);
	y_L_t = fixed(12);
	tau_r_t = fixed(13);
	un_t = fixed(14);
	g_t = fixed(15);
	xi_F = fixed(16);
	eta_F = fixed(17);
	xi_I = fixed(16);
	eta_I = fixed(19);
	omega = fixed(20);
	T = fixed(21);
	L_F_t = fixed(22);
	L_I_t = fixed(23);

	% Approximate around
	a = point .* [1; (1/0.075); 0.5*sqrt(A/B)];

	[aa bb] = nocbdc_get_coeffs(fixed,a);
	options = optimset('Display', 'off', 'TolX', 1e-8, 'TolFun', 1e-8,...
		'MaxFunEvals',1e4,'MaxIter',1e4);
	x0 = a;
	[x1 Fval] = fsolve(@(x) bb + aa*[x(1)-a(1);x(2)-a(2);x(3)-a(3)], x0, options);

	c_0 = x1(1);
	k_0 = x1(2);
	m_c_0 = x1(3);

	if any([c_0 k_0 m_c_0]<=0)
		display("WARNING: One initial guess value is non-positive.")
		if c_0<0;c_0=1;end;if k_0<0;k_0=1;end;if m_c_0<0;m_c_0=1;end;
		return
	end