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

function [c_0 k_0 i_dc_0] = cbdc_quant_rule_initial_guess(fixed, point, m_dc_to_gdp)
 	

% syms A B tau_c_t delta theta i_t i_dc_t ups_p_t betta alpha_k mc_t a_t...
% y_L_t tau_r_t un_t g_t xi_F eta_F xi_I eta_I omega T L_F_t L_I_t
A = fixed(1);
B = fixed(2);
tau_c_t = fixed(3);
delta = fixed(4);
theta = fixed(5);
i_t = fixed(6);
m_dc_to_gdp = fixed(7);
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

syms c_t_s k_t_s i_dc_t_s

assume(c_t_s,'positive');assume(k_t_s,'positive');assume(i_dc_t_s,'positive');

f1(c_t_s,k_t_s,i_dc_t_s) = ((c_t_s + (delta*k_t_s))/(1-g_t) + c_t_s*((1+tau_c_t)*(A*((1+tau_c_t)*c_t_s/((m_dc_to_gdp*((c_t_s+delta*k_t_s)/(1-g_t)))^(theta)*(T^(theta)+((i_t-i_dc_t_s)/(T^(theta)*i_t))^(theta/(1-theta)))))+B*(((m_dc_to_gdp*((c_t_s+delta*k_t_s)/(1-g_t)))^(theta)*(T^(theta)+((i_t-i_dc_t_s)/(T^(theta)*i_t))^(theta/(1-theta))))/((1+tau_c_t)*c_t_s))-2*sqrt(A*B))) + (delta*k_t_s)*(A*((delta*k_t_s)/((m_dc_to_gdp*((c_t_s+delta*k_t_s)/(1-g_t)))^(theta)*(T^(theta)+((i_t-i_dc_t_s)/(T^(theta)*i_t))^(theta/(1-theta))))) + B*(((m_dc_to_gdp*((c_t_s+delta*k_t_s)/(1-g_t)))^(theta)*(T^(theta)+((i_t-i_dc_t_s)/(T^(theta)*i_t))^(theta/(1-theta))))/(delta*k_t_s)) - 2*sqrt(A*B)) + un_t^(omega/(omega-1))*(xi_F*(L_F_t*eta_F)^(1/(1-omega))+xi_I*(L_I_t*eta_I)^(1/(1-omega)))) - (a_t/ups_p_t)*(k_t_s)^(alpha_k)*(y_L_t)^(1-alpha_k);
f2(c_t_s,k_t_s,i_dc_t_s) = (A*(((1+tau_c_t)*c_t_s/((m_dc_to_gdp*((c_t_s+delta*k_t_s)/(1-g_t)))^(theta)*(T^(theta)+((i_t-i_dc_t_s)/(T^(theta)*i_t))^(theta/(1-theta)))))^2+((delta*k_t_s)/((m_dc_to_gdp*((c_t_s+delta*k_t_s)/(1-g_t)))^(theta)*(T^(theta)+((i_t-i_dc_t_s)/(T^(theta)*i_t))^(theta/(1-theta)))))^2) - 2*B)*theta*((m_dc_to_gdp*((c_t_s+delta*k_t_s)/(1-g_t)))*((i_t-i_dc_t_s)/(T^(theta)*i_t))^(1/(1-theta)))^(theta-1) - i_t/(1+i_t);
f3(c_t_s,k_t_s,i_dc_t_s) = k_t_s - ((a_t*alpha_k*mc_t*y_L_t^(1-alpha_k))/(ups_p_t*((1+2*(A*((delta*k_t_s)/((m_dc_to_gdp*((c_t_s+delta*k_t_s)/(1-g_t)))^(theta)*(T^(theta)+((i_t-i_dc_t_s)/(T^(theta)*i_t))^(theta/(1-theta)))))-sqrt(A*B)))*((1/betta)-1+delta)/(1-tau_r_t))))^(1/(1-alpha_k));

df1_c(c_t_s,k_t_s,i_dc_t_s) = diff(f1, c_t_s);
df1_m_dc(c_t_s,k_t_s,i_dc_t_s) = diff(f1, i_dc_t_s);
df1_k(c_t_s,k_t_s,i_dc_t_s) = diff(f1, k_t_s);

df2_c(c_t_s,k_t_s,i_dc_t_s) = diff(f2, c_t_s);
df2_m_dc(c_t_s,k_t_s,i_dc_t_s) = diff(f2, i_dc_t_s);
df2_k(c_t_s,k_t_s,i_dc_t_s) = diff(f2, k_t_s);

df3_c(c_t_s,k_t_s,i_dc_t_s) = diff(f3, c_t_s);
df3_m_dc(c_t_s,k_t_s,i_dc_t_s) = diff(f3, i_dc_t_s);
df3_k(c_t_s,k_t_s,i_dc_t_s) = diff(f3, k_t_s);

a = point .* [1; (1/0.075); i_t/10];

f1l(c_t_s,k_t_s,i_dc_t_s) = subs(f1(a(1),a(2),a(3)) + df1_c(a(1),a(2),a(3))*(c_t_s-a(1)) + df1_k(a(1),a(2),a(3))*(k_t_s-a(2)) + df1_m_dc(a(1),a(2),a(3))*(i_dc_t_s-a(3)));
f2l(c_t_s,k_t_s,i_dc_t_s) = subs(f2(a(1),a(2),a(3)) + df2_c(a(1),a(2),a(3))*(c_t_s-a(1)) + df2_k(a(1),a(2),a(3))*(k_t_s-a(2)) + df2_m_dc(a(1),a(2),a(3))*(i_dc_t_s-a(3)));
f3l(c_t_s,k_t_s,i_dc_t_s) = subs(f3(a(1),a(2),a(3)) + df3_c(a(1),a(2),a(3))*(c_t_s-a(1)) + df3_k(a(1),a(2),a(3))*(k_t_s-a(2)) + df3_m_dc(a(1),a(2),a(3))*(i_dc_t_s-a(3)));

[c_0 k_0 i_dc_0] = vpasolve([f1l==0,f2l==0,f3l==0],[c_t_s k_t_s i_dc_t_s],a);

if any([c_0 k_0 i_dc_0]<=0)
	% display("WARNING: One initial guess value is non-positive.");
	if c_0<0;c_0=1;end;if k_0<0;k_0=1;end;if i_dc_0<0;i_dc_0=0;end;
	return
end
