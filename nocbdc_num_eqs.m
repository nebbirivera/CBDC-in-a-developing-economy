% Copyright (C) 2021-22 Pablo Nebbi Rivera Moreno and Karol Lorena Triana Montaño
%
% This is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% It is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% For a copy of the GNU General Public License, see <http://www.gnu.org/licenses/>.

function F = nocbdc_num_eqs(x, fixed)%, price_rule)

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
	xi_I = fixed(18);
	eta_I = fixed(19);
	omega = fixed(20);
	T = fixed(21);
	L_F_t = fixed(22);
	L_I_t = fixed(23);
	
	c_t = x(1);
	k_t = x(2);
	m_c_t = x(3);

	F(1) = ((c_t + (delta*k_t))/(1-g_t) + c_t*((1+tau_c_t)*(A*((1+tau_c_t)*c_t/(m_c_t^(theta)))+B*((m_c_t^(theta))/((1+tau_c_t)*c_t))-2*sqrt(A*B))) + (delta*k_t)*(A*((delta*k_t)/(m_c_t^(theta))) + B*((m_c_t^(theta))/(delta*k_t)) - 2*sqrt(A*B)) + un_t^(omega/(omega-1))*(xi_F*(L_F_t*eta_F)^(1/(1-omega))+xi_I*(L_I_t*eta_I)^(1/(1-omega)))) - (a_t/ups_p_t)*(k_t)^(alpha_k)*(y_L_t)^(1-alpha_k);
	F(2) = (A*(((1+tau_c_t)*c_t/(m_c_t^(theta)))^2+((delta*k_t)/(m_c_t^(theta)))^2) - 2*B)*theta*m_c_t^(theta-1) - i_t/(1+i_t);
	F(3) = k_t - ((a_t*alpha_k*mc_t*y_L_t^(1-alpha_k))/(ups_p_t*((1+2*(A*((delta*k_t)/(m_c_t^(theta)))-sqrt(A*B)))*((1/betta)-1+delta)/(1-tau_r_t))))^(1/(1-alpha_k));
  