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

var
	gdp_t (long_name='Real Gross Domestic Product')
	pi_t (long_name='General price inflation rate')
	L_t (long_name='Employed labor')
	m_c_t (long_name='Traditional Central Bank currency')
	G_t (long_name='Government expenditure')
	c_t (long_name='Consumption')
	infor_rate_t (long_name='Informality rate')
	Inv_t (long_name='Investment flow')
	unemp_rate_t (long_name='Unemployment rate')
	lab_force_part_rate_t (long_name='Labor force participation rate')
	vac_t (long_name='Aggregate average vacancies')
	X_t (long_name='Jaimovich-Rebelo dynamic shifter')
	k_t (long_name='Capital stock')
	L_F_t (long_name='Formal labor')
	L_I_t (long_name='Informal labor')
	un_t (long_name='Unemployment')
	O_t (long_name='Out of the labor force')
	w_F_t (long_name='Formal sector real wage')
	w_I_t (long_name='Informal sector real wage')
	r_t (long_name='Real return on capital')
	lambda_c_t (long_name='Budget constraint lagrange multiplier')
	lambda_Inv_t (long_name='Capital accumulation lagrange multiplier')
	lambda_F_t (long_name='Formal labor law of motion lag. multiplier')
	lambda_I_t (long_name='Informal labor law of motion lag. multiplier')
	lambda_X_t (long_name='Jaimovich-Rebelo shifter lag. multiplier')
	u_c_t (long_name='Consumption marginal utility')
	i_t (long_name='Main monetary policy nominal net interest rate')
	pi_opt_t (long_name='Optimal price inflation rate')
	sc_t (long_name='Consumption transaction costs')
	sInv_t (long_name='Investment trasaction costs')
	vc_t (long_name='Consumption based transaction velocity')
	vInv_t (long_name='Investment based transaction velocity')
	Ac_t (long_name='Consumption based velocity coefficient')
	AInv_t (long_name='Investment based velocity coefficient')
	S_TC_t (long_name='Shock to demand for transaction costs')
	LGF_t (long_name='Liquidity generating function')
	LGF_m_c_t (long_name='LGF derivative w.r.t. m_c_t')
	a_t (long_name='General TFP shock')
	a_F_t (long_name='Formal sector productivity shock')
	V_F_h_t (long_name='Value accrued to HH in the formal sector')
	V_F_f_t (long_name='Value accrued to firm in the formal sector')
	V_I_h_t (long_name='Value accrued to HH in the informal sector')
	V_I_f_t (long_name='Value accrued to firm in the informal sector')
	vac_F_t (long_name='Formal sector vacancies')
	vac_I_t (long_name='Informal sector vacancies')
	q_F_t (long_name='Prob. a formal vacancy finds a worker')
	q_I_t (long_name='Prob. an informal vacancy finds a worker')
	pr_F_t (long_name='Prob. a worker finds a vacancy in the formal sector')
	pr_I_t (long_name='Prob. a worker finds a vacancy in the informal sector')
	Y_t (long_name='Output')
	p_L_t (long_name='Labor final price *')
	p_F_t (long_name='Formal product price')
	p_I_t (long_name='Informal product price')
	y_L_t (long_name='Labor output *')
	mc_t (long_name='Real aggregate marginal cost')
	SDF_t (long_name='Stochastic discount factor')
	y_F_t (long_name='Formal sector output')
	y_I_t (long_name='Informal sector output')
	ups_p_t (long_name='General price dispersion')
	aux_p_1_t (long_name='Optimal price aux. variable no. 1')
	aux_p_2_t (long_name='Optimal price aux. variable no. 2')
	z_u_t (long_name='Preferences shock')
	tau_w_t (long_name='Payroll taxes to formal firms')
	tau_c_t (long_name='VAT rate')
	tau_r_t (long_name='Capital return rate')
	g_t (long_name='Fiscal expenditure')
	d_t (long_name='Government debt')
	tau_ls_t (long_name='Lump-sum taxes')
	match_t (long_name='Job creation')
	TFP_t (long_name='Total factor productivity')
	tc_avg_t (long_name='Average transaction costs distortion')
	w_avg_t (long_name='Average real wage of the economy')
;

parameters
	A (long_name='Transaction costs A parameter')
	B (long_name='Transaction costs B parameter')
	gama (long_name='Jaimovich-Rebelo shift parameter')
	betta (long_name='Household deterministic discount factor')
	pssi (long_name='Disutility of labor parameter')
	chi (long_name='Inverse of Frisch elasticity')
	phi_Inv (long_name='Ivestment adjustment cost parameter')
	delta (long_name='Capital depreciation rate')
	eta_I (long_name='Informal sector separation rate')
	eta_F (long_name='Formal sector separation rate')
	muu (long_name='Firms relative bargaining power')
	xi_I (long_name='Informal sector cost of a vacancy')
	xi_F (long_name='Formal sector cost of a vacancy')
	phi_p (long_name='Prob. of not being able to change price')
	zzeta (long_name='Disutility of unemployment parameter') 
	alpha_k (long_name='YEl w.r.t. capital')
	eps_p (long_name='EoS between intermediate varieties')
	eps_L (long_name='EoS between types of labor firms')
	theta (long_name='LGF el. w.r.t. each money')
	pi_obj (long_name='Long term objective inflation rate')
	i_ss (long_name='Steady state policy interest rate')
	phi_pi_i (long_name='Inflation reaction param. main monetary policy tool')
	phi_y_i (long_name='Output gap reaction param. main monetary policy tool')
	phi_pi_2 (long_name='Inflation reaction param. 2nd monetary policy tool')
	phi_y_2 (long_name='Output gap reaction param. 2nd monetary policy tool')
	T (long_name='Relative liquidity services of digital currency')
	omega (long_name='El. of matching w.r.t. unemployment')
	a_I (long_name='Average productivity of the informal sector')
	spr (long_name='Spread between CDBD and regular MP interest')
	g_ss (long_name='Steady government expenditure as a share of GDP')
	tau_c_ss (long_name='Steady state VAT')
	tau_r_ss (long_name='Steady state return of capital tax rate')
	tau_w_ss (long_name='Steady state payroll tax rate')
	rho_a (long_name='Persistance of technology shocks')
	rho_a_F (long_name='Persistance of informal technology shocks')
	rho_S_TC (long_name='Persistance of TC demand shocks')
	rho_g (long_name='Persistance of government expenditure shock')
	rho_u (long_name='Persistance of demand (preferences) shock')
	rho_tau_w (long_name='Persistance of payroll tax disturbances')
	rho_tau_c (long_name='Persistance of consumption tax disturbances')
	rho_tau_r (long_name='Persistance of capital returns disturbances')
	sig_a (long_name='TFP shock standard deviation')
	sig_g (long_name='Government expenditure shock standard deviation')
	sig_u (long_name='Demand (preferences) shock standard deviation')
	sig_a_F (long_name='Formal productivity shock standard deviation')
	sig_i (long_name='Monetary policy shock standard deviation')
	c_ss (long_name='Levels steady state consumption')
	k_ss (long_name='Levels steady state capital')
	r_ss (long_name='Steady state real return of capital')
	w_F_ss (long_name='Levels steady state formal sector wage')
	L_F_ss (long_name='Levels steady state formal labor')
	tau_ls_ss (long_name='Levels steady state lump-sum taxes')
	d_ss (long_name='Levels steady state government debt')
	gdp_ss (long_name='GDP steady state')
	m_c_ss (long_name='Levels steady state traditional money')
	gamma_g (long_name='Fiscal rule parameter')
;

varexo e_i_t e_S_TC_t e_a_t e_a_F_t e_z_u_t e_g_t;

load fixed_params.mat % Assigned in main.m
alpha_k=x_alpha_k;
eps_p=x_eps_p;
delta=x_delta;
betta=x_betta;
chi=x_chi;
pi_obj=x_pi_obj;
i_ss=x_i_ss;
eta_I=x_eta_I;
sig_S_TC=x_sig_S_TC;
rho_S_TC=x_rho_S_TC;
rho_a_F=x_rho_a_F;
rho_u=x_rho_u;
rho_a=x_rho_a;
sig_a_F=x_sig_a_F;
rho_g=x_rho_g;
sig_u=x_sig_u;
rho_tau_w=x_rho_tau_w;
rho_tau_c=x_rho_tau_c;
rho_tau_r=x_rho_tau_r;
omega=x_omega;
muu=x_muu;
eps_L=x_eps_L;
theta=x_theta;
sig_i=x_sig_i;
B=x_B;
A=x_A;
phi_p=x_phi_p;
a_I=x_a_I;
g_ss=x_g_ss;
tau_w_ss=x_tau_w_ss;
tau_c_ss=x_tau_c_ss;
tau_r_ss=x_tau_r_ss;
phi_pi_i=x_phi_pi_i;
phi_y_i=x_phi_y_i;
phi_pi_2=x_phi_pi_2;
phi_y_2=x_phi_y_2;
gamma_g=x_gamma_g;

% Calibrated to second moments
load calibrated_parameters_2ndm.mat
sig_g = x1(1);
eta_F = x1(2);
pssi = x1(3);
phi_Inv = x1(4);
gama = x1(5);
sig_a = x1(6);

% Calibrated to steady state moments
load calibrated_parameters_ssm
xi_F = xi_F_c;
xi_I = xi_I_c;
zzeta = zzeta_c;

model;
	[name='1. Consumption FOC']
	exp(lambda_c_t)*(1+tau_c_t)*(1+2*(Ac_t*exp(vc_t)-sqrt(A*B))) = exp(u_c_t) 
	+ gama*lambda_X_t*(exp(X_t(-1))/exp(c_t))^(1-gama);

	[name='2. Jaimovich-Rebelo shifter FOC']
	lambda_X_t = betta*(1-gama)*lambda_X_t(+1)*(exp(c_t(+1))/exp(X_t))^(gama)
	- pssi*((exp(L_t))^(1+chi)/(1+chi))*exp(u_c_t);

	[name='3. Investment FOC']
	exp(lambda_c_t)*(1+2*(AInv_t*exp(vInv_t)-sqrt(A*B))) = 
	exp(lambda_Inv_t)
	*(1-phi_Inv*(exp(Inv_t)/exp(Inv_t(-1))*((exp(Inv_t)/exp(Inv_t(-1)))-1)
		+0.5*(((exp(Inv_t)/exp(Inv_t(-1)))-1)^2))) 
	+ phi_Inv*betta*exp(lambda_Inv_t(+1))*((exp(Inv_t(+1))/exp(Inv_t))^2)
	*((exp(Inv_t(+1))/exp(Inv_t))-1);

    [name='4. Capital FOC']
    exp(lambda_Inv_t) = betta*(exp(lambda_c_t(+1))*(1-tau_r_t(+1))*r_t(+1) 
    + exp(lambda_Inv_t(+1))*(1-delta));

    [name='5. Bonds FOC']
    betta*exp(lambda_c_t(+1))*(1+i_t) = exp(lambda_c_t)*(1+pi_t(+1));

    [name='6. Value accrued to HH for working in the informal sector']
    exp(V_I_h_t) = exp(w_I_t) 
    - pssi*exp(X_t)*(exp(L_t))^chi*exp(u_c_t)/exp(lambda_c_t) 
    + (1-eta_I)*exp(SDF_t)*exp(V_I_h_t(+1));

    [name='7. Value accrued to firm from an emp. relation informal sector']
    exp(V_I_f_t) = exp(p_I_t)*a_I - exp(w_I_t) 
    + (1-eta_I)*exp(SDF_t)*exp(V_I_f_t(+1));

    [name='8. Nash bargaining in the informal sector']
	muu*exp(V_I_h_t) = (1-muu)*exp(V_I_f_t);

	[name='9. Household value in the informal sector definition']
	exp(lambda_I_t) = exp(V_I_h_t)*exp(lambda_c_t);

	[name='10. Informal firm value definition']
	xi_I = exp(V_I_f_t)*exp(q_I_t);

	[name='11. Value accrued to HH for working in the formal sector']
	exp(V_F_h_t) = exp(w_F_t) 
	- pssi*exp(X_t)*exp(L_t)^(chi)*exp(u_c_t)/exp(lambda_c_t) 
	+ (1-eta_F)*exp(SDF_t)*exp(V_F_h_t(+1));

	[name='12. Value accrued to a formal firm from employment relation']
	exp(V_F_f_t) = exp(p_F_t)*a_F_t - exp(w_F_t)*(1+tau_w_t) 
	+ (1-eta_F)*exp(SDF_t)*exp(V_F_f_t(+1));

	[name='13. Nash bargaining in the formal sector']
	muu*exp(V_F_h_t) = (1-muu)*exp(V_F_f_t);

	[name='14. Household value in the formal sector definition']
	exp(lambda_F_t) = exp(V_F_h_t)*exp(lambda_c_t);

	[name='15. Formal firm value definition']
	xi_F = exp(V_F_f_t)*exp(q_F_t);

	[name='16. Unemployment FOC']
	zzeta*z_u_t*exp(un_t) = exp(pr_F_t)*exp(lambda_F_t) 
	+ exp(pr_I_t)*exp(lambda_I_t);

	[name='17. m_c_t FOC']
	(Ac_t*exp(vc_t)^2 + AInv_t*exp(vInv_t)^2 - 2*B)*exp(LGF_m_c_t) = 
	(i_t)/(1+i_t);

	[name='18. Aggregate vacancies']
	exp(vac_t) = exp(vac_F_t) + exp(vac_I_t);

	[name='19. Return on capital/firms demand for capital']
	r_t*exp(k_t(-1)) = alpha_k*exp(mc_t)*exp(Y_t);

	[name='20. Firms demand for aggregate labor']
	exp(p_L_t)*exp(y_L_t) = (1-alpha_k)*exp(mc_t)*exp(Y_t);

	[name='21. Demand for informal labor']
	exp(y_L_t) = exp(y_I_t)*(exp(p_L_t)/exp(p_I_t))^(-eps_L);

	[name='22. Demand for formal labor']
	exp(y_L_t) = exp(y_F_t)*(exp(p_L_t)/exp(p_F_t))^(-eps_L);

	[name='23. Optimal price setting']
	1+pi_opt_t = (eps_p/(eps_p-1))*(1+pi_t)*exp(aux_p_1_t)/exp(aux_p_2_t);

	[name='24. aux_p_1_t']
	exp(aux_p_1_t) = exp(mc_t)*exp(Y_t) + phi_p*exp(SDF_t)
	*(1+pi_t(+1))^(eps_p)*exp(aux_p_1_t(+1));

	[name='25. aux_p_2_t']
	exp(aux_p_2_t) = exp(Y_t) + phi_p*exp(SDF_t)*(1+pi_t(+1))^(eps_p-1)
	*exp(aux_p_2_t(+1));

	[name='26. Production function']
	exp(Y_t) = (a_t/exp(ups_p_t))*(exp(k_t(-1)))^(alpha_k)
	*(exp(y_L_t))^(1-alpha_k);

	[name='27. Price dispersion evolution']
	exp(ups_p_t) = (1+pi_t)^(eps_p)*((1-phi_p)*(1+pi_opt_t)^(-eps_p) 
	+ phi_p*exp(ups_p_t(-1)));

	[name='28. Resource constraint/market clearing']
	exp(Y_t) = (exp(c_t)+ exp(Inv_t))/(1-g_t) 
	+ exp(c_t)*(1+tau_c_t)*exp(sc_t)+exp(Inv_t)*exp(sInv_t) 
	+ xi_F*exp(vac_F_t) + xi_I*exp(vac_I_t);

	[name='29. Main monetary policy Taylor rule *']
	1+i_t = (1+i_ss)*((1+pi_t)/(1+pi_obj))^(phi_pi_i) 
	* exp(sig_i*e_i_t);

	[name='30. Real GDP']
	exp(gdp_t) = (exp(c_t)+exp(Inv_t))+exp(G_t);

	[name='31. Labor aggregator']
	exp(y_L_t) = (exp(y_F_t)^((eps_L-1)/eps_L) 
	+ exp(y_I_t)^((eps_L-1)/eps_L))^(eps_L/(eps_L-1));

	[name='32. Inflation evolution']
	(1+pi_t)^(1-eps_p) = phi_p + (1-phi_p)*(1+pi_opt_t)^(1-eps_p);

	[name='33. Household labor endowment']
	exp(L_F_t) + exp(L_I_t) + exp(un_t) + exp(O_t) = 1;

	[name='34. Formal output']
	exp(y_F_t) = a_F_t*exp(L_F_t);

	[name='35. Informal output']
	exp(y_I_t) = a_I*exp(L_I_t);

	[name='36. Consumption transaction costs']
	exp(sc_t) = Ac_t*exp(vc_t) + B/exp(vc_t) - 2*sqrt(A*B);

	[name='37. A coefficient']
	Ac_t = S_TC_t*A;

	[name='38. Consumption transaction costs']
	exp(sInv_t) = AInv_t*exp(vInv_t) + B/exp(vInv_t) - 2*sqrt(A*B);

	[name='39. A coefficient']
	AInv_t = S_TC_t*A;

	[name='40. LGF']
	exp(LGF_t) = (exp(m_c_t))^(theta);

	[name='41. LGF derivative w.r.t. exp(m_c_t)']
	exp(LGF_m_c_t) = theta*(exp(m_c_t))^(theta-1);
	
	[name='42. Government expenditure']
	exp(G_t) = g_t*exp(gdp_t);

	[name='43. Rate at which a vac. matches with a worker in the formal sector']
	exp(q_F_t) = (exp(un_t)/exp(vac_F_t))^(omega);

	[name='44. Rate at which a vac. matches with a worker in the informal sector']
	exp(q_I_t) = (exp(un_t)/exp(vac_I_t))^(omega);

	[name='45. Rate at which a worker finds a vacancy in the formal sector']
	exp(pr_F_t) = (exp(un_t)/exp(vac_F_t))^(omega-1);

	[name='46. Rate at which a worker finds a vacancy in the informal sector']
	exp(pr_I_t) = (exp(un_t)/exp(vac_I_t))^(omega-1);

	[name='47. AR(1) General TFP shock']
	log(a_t) = (rho_a)*log(a_t(-1)) + (1-rho_a)*log(1) + sig_a*e_a_t;

	[name='48. AR(1) Formal sector productivity shock']
	log(a_F_t) = rho_a_F*log(a_F_t(-1)) + e_a_F_t*sig_a_F;

	[name='49. AR(1) General TC demand shock']
	log(S_TC_t) = rho_S_TC*log(S_TC_t(-1)) + e_S_TC_t;

	[name='50. Investment process']
	exp(k_t) = (1-delta)*exp(k_t(-1)) 
	+ exp(Inv_t)*(1-(phi_Inv/2)*(exp(Inv_t)/exp(Inv_t(-1))-1)^2);

	[name='51. Households stochastic discount factor']
	exp(SDF_t) = betta*exp(lambda_c_t(+1))/exp(lambda_c_t);

	[name='52. Jaimovich-Rebelo dynamic shifter']
	exp(X_t) = exp(c_t)^(gama)*(exp(X_t(-1)))^(1-gama);

	[name='53. Marginal utility of contemporaneous consumption']
	exp(u_c_t) = z_u_t/(exp(c_t) - pssi*exp(X_t)*((exp(L_t))^(1+chi))/(1+chi));

	[name='54. Formal labor law of motion']
	exp(L_F_t) = (1-eta_F)*exp(L_F_t(-1)) + exp(q_F_t)*exp(vac_F_t);

	[name='55. Informal labor law of motion']
	exp(L_I_t) = (1-eta_I)*exp(L_I_t(-1)) + exp(q_I_t)*exp(vac_I_t);

	[name='56. Consumption based velocity']
	exp(vc_t) = (1+tau_c_t)*exp(c_t)/exp(LGF_t);

	[name='57. Investment based velocity']
	exp(vInv_t) = exp(Inv_t)/exp(LGF_t);

	[name='58. Labor force']
	exp(L_t) = exp(L_I_t) + exp(L_F_t);

	[name='59. AR(1) Preferences shock process']
	log(z_u_t) = rho_u*log(z_u_t(-1)) + e_z_u_t*sig_u;

	[name='60. AR(1) Formal sector payroll taxes']
	tau_w_t = rho_tau_w*tau_w_t(-1) + (1-rho_tau_w)*tau_w_ss;

	[name='61. AR(1) VAT']
	tau_c_t = rho_tau_c*tau_c_t(-1) + (1-rho_tau_c)*tau_c_ss;

	[name='62. AR(1) Capital return rate']
	tau_r_t = rho_tau_r*tau_r_t(-1) + (1-rho_tau_r)*tau_r_ss;

	[name='63. AR(1) Discretional government expenditure']
	g_t = rho_g*g_t(-1) + (1-rho_g)*g_ss + sig_g*e_g_t;

	[name='64. Unemployment rate definition']
	unemp_rate_t = exp(un_t)/(exp(L_t)+exp(un_t));

	[name='65. Informality rate']
	infor_rate_t = exp(L_I_t)/exp(L_t);
	
	[name='66. Labor force particiation rate']
	lab_force_part_rate_t = exp(L_t) + exp(un_t);

	[name='67. Government budget constraint']
	tau_ls_t + tau_c_t*exp(c_t) + tau_r_t*r_t*exp(k_t(-1)) 
	+ tau_w_t*exp(w_F_t)*exp(L_F_t) + exp(m_c_t) + d_t = 
	(1+pi_t)^(-1)*((1+i_t(-1))*d_t(-1) + exp(m_c_t(-1))) + exp(G_t);

	[name='68. Fiscal rule']
	tau_ls_t + tau_c_t*exp(c_t) + tau_r_t*r_t*exp(k_t(-1)) 
	+ tau_w_t*exp(w_F_t)*exp(L_F_t) 
	-(tau_ls_ss+tau_c_ss*exp(c_ss)+tau_r_ss*r_ss*exp(k_ss)
		+tau_w_ss*exp(w_F_ss)*exp(L_F_ss)) = 
	gamma_g*((1+pi_t)^(-1)*((1+i_t(-1))*d_t(-1)+exp(m_c_t(-1)))
		- ((1+pi_obj)^(-1)*((1+i_ss)*d_ss+exp(m_c_ss))));

	[name='69. Aggregate matchings']
	exp(match_t) = (exp(un_t))^(omega)*(exp(vac_F_t))^(1-omega) 
	+ (exp(un_t))^(omega)*(exp(vac_I_t))^(1-omega);

	[name='70. Total factor productivity']
	exp(TFP_t) = exp(Y_t)/((exp(k_t(-1)))^(alpha_k)*(exp(L_t))^(1-alpha_k));

	[name='71. Average transaction costs distortion']
	exp(tc_avg_t) = ((1+2*(Ac_t*exp(vc_t)-sqrt(A*B)))*exp(c_t)
		+(1+2*(AInv_t*exp(vInv_t)-sqrt(A*B)))*exp(Inv_t))
	/(exp(c_t)+exp(Inv_t));

	[name='72. Average real wage of the economy']
	exp(w_avg_t) = exp(w_F_t)*(1-infor_rate_t)+exp(w_I_t)*infor_rate_t;
end;

steady;
shocks;
	var e_i_t; stderr 1; 
	var e_S_TC_t; stderr 0; 
	var e_a_t; stderr 1; 
	var e_a_F_t; stderr 1; 
	var e_z_u_t; stderr 1; 
	var e_g_t; stderr 1; 
end;


options_.noprint = 1;
options_.nograph = 1;

stoch_simul(hp_filter=1600, order=1, irf=20);

nocbdc_full_ss = oo_.dr.ys;
save nocbdc_full_ss.mat nocbdc_full_ss;
