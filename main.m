% This script runs the subroutines that generate the results and plots shown in
% the document.
% The .m files were run in MATLAB R2021b, and .mod files were run using Dynare
% 4.6.4

clear;close;clc;
tic

x_alpha_k = 1/3;
x_eps_p = 5;
x_delta = (1.1)^(1/4)-1;
x_betta = 0.99;
x_chi = 1;
x_pi_obj = (1.03)^(1/4)-1; % Banco de la República target inflation
x_i_ss = (1+x_pi_obj)/x_betta-1; % Compl. w. Banco de la República tgt. inflat.
x_eta_I = 1;
x_sig_S_TC = 0;
x_rho_S_TC = 0.9;
x_rho_a_F = 0.9;
x_rho_u = 0.9;
x_rho_a = 0.9; 
x_sig_a_F = 0.001; 
x_rho_g = 0.9; 
x_sig_u = 0.0; 
x_rho_tau_w = 0.9;
x_rho_tau_c = 0.9;
x_rho_tau_r = 0.9;
x_omega = 0.5; % Blanchard and Galí (2010)
x_muu = 0.162; % González et al. (2012)
x_eps_L = 8; % Restrepo-Echavarría (2014)
x_theta = 0.95; % Barrdear-Kumhof (2020) 
x_sig_i = 0.001; % González et al. (2012)
x_B = 0.07524; % Schmitt-Grohé and Uribe 2004
x_A = 0.111; % Schmitt-Grohé and Uribe 2004
x_phi_p = 0.202; % González et al. (2012)
x_a_I = 0.67; % Alberola and Urrutia (2020)
x_g_ss = 0.143787095699914; % Average govt. expend. to GDP ratio
x_tau_w_ss = 0.02366461; % Calibrated using Mendoza's (1994) methodology
x_tau_c_ss = 0.032094133; % Calibrated using Mendoza's (1994) methodology
x_tau_r_ss = 0.030558439; % Calibrated using Mendoza's (1994) methodology
% Policy coefficients
x_phi_pi_i = 2;
x_phi_y_i = 0;
x_phi_pi_2 = 1;
x_phi_y_2 = 0;
x_gamma_g = 1; % CBDC parameters
x_i_spread = 0.015;
x_spr = x_i_ss/(x_i_ss-x_i_spread);
x_T = 1.5;
x_m_dc_to_gdp = 0.2;

save fixed_params.mat x_alpha_k x_eps_p x_delta x_betta x_chi x_pi_obj...
x_i_ss x_eta_I x_sig_S_TC x_rho_S_TC x_rho_a_F x_rho_u x_rho_a x_sig_a_F...
x_rho_g x_sig_u x_rho_tau_w x_omega x_muu x_eps_L x_theta x_sig_i x_B...
x_A x_phi_p x_a_I x_g_ss x_tau_w_ss x_tau_c_ss x_tau_r_ss x_rho_tau_c...
x_rho_tau_r x_phi_pi_i x_phi_y_i x_phi_pi_2 x_phi_y_2 x_gamma_g x_spr x_T...
x_i_spread x_m_dc_to_gdp

%%%% IRFs analysis
business_cycle_analysis
%% Compute cbdc quant. Rule steady state and store it for computation 
%% efficiency
% dynare('cbdc_quant_rule.mod')
% cbdc_ss_ig = oo_.dr.ys
% save cbdc_qr_ss_initial_guess.mat cbdc_ss_ig

%%%% Sacrifice ratios analysis
sacrifice_ratios
% dynare cbdc_price_rule

toc
