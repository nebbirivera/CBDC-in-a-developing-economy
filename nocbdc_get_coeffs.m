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

function [aa bb] = nocbdc_get_coeffs(fixed,a)

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

 
	c_t_s = a(1);
	k_t_s = a(2);
	m_c_t_s = a(3);
	
	bb = [...
	un_t^(omega/(omega-1))*(xi_F/(L_F_t*eta_F)^(1/(omega-1))+xi_I/(L_I_t*eta_I)^(1/(omega-1)))-(c_t_s+delta*k_t_s)/(g_t-1)+delta*k_t_s*((A*delta*k_t_s)/m_c_t_s^theta-2*(A*B)^(1/2)+(B*m_c_t_s^theta)/(delta*k_t_s))+c_t_s*(tau_c_t+1)*((A*c_t_s*(tau_c_t+1))/m_c_t_s^theta-2*(A*B)^(1/2)+(B*m_c_t_s^theta)/(c_t_s*(tau_c_t+1)))-(a_t*k_t_s^alpha_k*y_L_t^(1-alpha_k))/ups_p_t;...
	-i_t/(i_t+1)-m_c_t_s^(theta-1)*theta*(2*B-A*(delta^2*k_t_s^2/m_c_t_s^(2*theta)+c_t_s^2/m_c_t_s^(2*theta)*(tau_c_t+1)^2));...
	k_t_s-1/(-(a_t*alpha_k*mc_t*y_L_t^(1-alpha_k)*(tau_r_t-1))/(ups_p_t*((2*A*delta*k_t_s)/m_c_t_s^theta-2*(A*B)^(1/2)+1)*(delta+1/betta-1)))^(1/(alpha_k-1))];


	aa = [...
	(tau_c_t+1)*((A*c_t_s*(tau_c_t+1))/m_c_t_s^theta-2*(A*B)^(1/2)+(B*m_c_t_s^theta)/(c_t_s*(tau_c_t+1)))-1/(g_t-1)+c_t_s*((A*(tau_c_t+1))/m_c_t_s^theta-(B*m_c_t_s^theta)/(c_t_s^2*(tau_c_t+1)))*(tau_c_t+1),delta*((A*delta*k_t_s)/m_c_t_s^theta-2*(A*B)^(1/2)+(B*m_c_t_s^theta)/(delta*k_t_s))-delta/(g_t-1)+delta*k_t_s*((A*delta)/m_c_t_s^theta-(B*m_c_t_s^theta)/(delta*k_t_s^2))-(a_t*alpha_k*k_t_s^(alpha_k-1)*y_L_t^(1-alpha_k))/ups_p_t,delta*k_t_s*((B*m_c_t_s^(theta-1)*theta)/(delta*k_t_s)-(A*delta*k_t_s*theta)/m_c_t_s^(theta+1))-c_t_s*((A*c_t_s*theta*(tau_c_t+1))/m_c_t_s^(theta+1)-(B*m_c_t_s^(theta-1)*theta)/(c_t_s*(tau_c_t+1)))*(tau_c_t+1);...
	2*A*c_t_s/m_c_t_s^(2*theta)*m_c_t_s^(theta-1)*theta*(tau_c_t+1)^2,2*A*delta^2*k_t_s/m_c_t_s^(2*theta)*m_c_t_s^(theta-1)*theta,-A*m_c_t_s^(theta-1)*theta*((2*delta^2*k_t_s^2*theta)/m_c_t_s^(2*theta+1)+(2*c_t_s^2*theta*(tau_c_t+1)^2)/m_c_t_s^(2*theta+1))-m_c_t_s^(theta-2)*theta*(2*B-A*(delta^2*k_t_s^2/m_c_t_s^(2*theta)+c_t_s^2/m_c_t_s^(2*theta)*(tau_c_t+1)^2))*(theta-1);...
	0,(2*A*a_t*alpha_k*delta*mc_t*y_L_t^(1-alpha_k)*(tau_r_t-1))/(m_c_t_s^theta*ups_p_t*(alpha_k-1)*((2*A*delta*k_t_s)/m_c_t_s^theta-2*(A*B)^(1/2)+1)^2*(delta+1/betta-1)*(-(a_t*alpha_k*mc_t*y_L_t^(1-alpha_k)*(tau_r_t-1))/(ups_p_t*((2*A*delta*k_t_s)/m_c_t_s^theta-2*(A*B)^(1/2)+1)*(delta+1/betta-1)))^(1/(alpha_k-1)+1))+1,-(2*A*a_t*alpha_k*delta*k_t_s*mc_t*theta*y_L_t^(1-alpha_k)*(tau_r_t-1))/(m_c_t_s^(theta+1)*ups_p_t*(alpha_k-1)*((2*A*delta*k_t_s)/m_c_t_s^theta-2*(A*B)^(1/2)+1)^2*(delta+1/betta-1)*(-(a_t*alpha_k*mc_t*y_L_t^(1-alpha_k)*(tau_r_t-1))/(ups_p_t*((2*A*delta*k_t_s)/m_c_t_s^theta-2*(A*B)^(1/2)+1)*(delta+1/betta-1)))^(1/(alpha_k-1)+1))];
