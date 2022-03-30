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


clear;% close;
% dynare('nocbdc.mod','nowarn')
% nocbdc_oo = oo_;
% dynare('cbdc_price_rule.mod','nowarn')
% cbdc_price_rule_oo = oo_;
% dynare('cbdc_quant_rule.mod','nowarn')
% cbdc_quant_rule_oo = oo_;

%% Set plot parameters
n_sr = 4; % Cumulative periods to calculate the sacrifice ratio
plot_line_width = 1.5;
title_font_size = 13;
n_points_qr = 30;
lege_font_size = 10;
label_font_size = 12;

% Compute the sacrifice ratio for the non-CBDC economy
load baseline_irfs_results.mat
gdp_cumul = cumsum(nocbdc_oo.irfs.gdp_t_e_i_t);
unemp_rate_cumul = cumsum(nocbdc_oo.irfs.unemp_rate_t_e_i_t);
pi_cumul = cumsum(nocbdc_oo.irfs.pi_t_e_i_t);
nocbdc_sacrifice_ratio = gdp_cumul(n_sr)/pi_cumul(n_sr)

load results_sacrifice_ratios.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the sacrifice ratio for different T and i_spread CBDC price rule
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dynare('cbdc_price_rule.mod','nowarn')
n_points = 10;
prop_factor = 100;
% spreads = linspace(0.008,0.0176,n_points);
% cbdc_pr_sacrifice_ratio_i_sp = NaN(1,n_points);
% % Loop through i_spread values CBDC-pr
% for k=1:n_points
% 	set_param_value("i_spread",spreads(k));
% 	[info, oo_, options_] = stoch_simul(M_, options_, oo_, var_list_);
% 	i_spread
% 	if info==0
% 		gdp_cumul_cbdc_pr = cumsum(oo_.irfs.gdp_t_e_i_t);
% 		pi_cumul_cbdc_pr = cumsum(oo_.irfs.pi_t_e_i_t);
% 		cbdc_pr_sacrifice_ratio_i_sp(k) = gdp_cumul_cbdc_pr(n_sr)/...
% 		pi_cumul_cbdc_pr(n_sr);
% 	else
% 		get_error_message(info)
% 	end
% end

% Plot GDP sacrifice ratio vs. i_spread values
figure;
subplot(1,2,1);
plot(prop_factor*spreads,cbdc_pr_sacrifice_ratio_i_sp,"blue",...
	"LineWidth",plot_line_width);
title("\textbf{Price rule sacrifice ratio}","Interpreter","latex",...
	"fontsize",title_font_size);
ylabel("GDP-Inflation sacrifice ratio","FontSize",label_font_size);
xlabel("% spread","FontSize",label_font_size);
yline(nocbdc_sacrifice_ratio,"--r","LineWidth",plot_line_width);
xlim(prop_factor*[spreads(1) spreads(length(spreads))]);
legend("CBDC - Price rule","No CBDC",...
	'Location',"southeast",'FontSize',lege_font_size);
legend boxon box on;grid on;


%% Loop through T values CBDC-pr
% Ts = linspace(1,2,n_points);
% cbdc_pr_sacrifice_ratio_T = NaN(1,n_points);
% for k=1:n_points
% 	set_param_value("T",Ts(k));
% 	[info, oo_, options_] = stoch_simul(M_, options_, oo_, var_list_);
% 	i_spread
% 	if info==0
% 		gdp_cumul_cbdc_pr = cumsum(oo_.irfs.gdp_t_e_i_t);
% 		pi_cumul_cbdc_pr = cumsum(oo_.irfs.pi_t_e_i_t);
% 		cbdc_pr_sacrifice_ratio_T(k) = gdp_cumul_cbdc_pr(n_sr)/...
% 		pi_cumul_cbdc_pr(n_sr);
% 	else
% 		get_error_message(info)
% 	end
% end

subplot(1,2,2);
plot(Ts,cbdc_pr_sacrifice_ratio_T,"blue","LineWidth",plot_line_width);
title("\textbf{Price rule sacrifice ratio}","Interpreter","latex",...
	"fontsize",title_font_size);
ylabel("GDP-Inflation sacrifice ratio","FontSize",label_font_size);
xlabel("Liquidity level of CBDC","FontSize",label_font_size);
yline(nocbdc_sacrifice_ratio,"--r","LineWidth",plot_line_width);
xlim([Ts(1) Ts(length(Ts))]);
legend("CBDC - Price rule","No CBDC",...%'Orientation','horizontal',...
	'Location',"southwest",'FontSize',lege_font_size);
legend boxon;box on;grid on;

set(gcf,"PaperPosition",[0 0 10 3.8]);
% saveas(gcf,strcat("price_rule_sr_i_sp_and_T",".eps"),"epsc");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the sacrifice ratio for different spreads CBDC quant. rule
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dynare('cbdc_quant_rule.mod','nowarn')
% cbdc_to_gdp_ratios = linspace(0.01,0.3,n_points_qr);
% cbdc_qr_sacrifice_ratio_T = NaN(1,n_points_qr);
% % Loop through m_dc_to_gdp values CBDC-qr
% for k=1:n_points_qr
% 	set_param_value("m_dc_to_gdp",cbdc_to_gdp_ratios(k));
% 	[info, oo_, options_] = stoch_simul(M_, options_, oo_, var_list_);
% 	cbdc_to_gdp_ratios(k)
% 	if info==0
% 		gdp_cumul_cbdc_pr = cumsum(oo_.irfs.gdp_t_e_g_t);
% 		pi_cumul_cbdc_pr = cumsum(oo_.irfs.pi_t_e_g_t);
% 		cbdc_qr_sacrifice_ratio_T(k) = gdp_cumul_cbdc_pr(n_sr)/...
% 		pi_cumul_cbdc_pr(n_sr);
% 	else
% 		get_error_message(info)
% 	end
% end
figure;
% subplot(1,1,3);
plot(cbdc_to_gdp_ratios,cbdc_qr_sacrifice_ratio_T,"blue","LineWidth",plot_line_width);
title("\textbf{Quantitative rule sacrifice ratio}","Interpreter","latex",...
	"fontsize",title_font_size);
ylabel("GDP-Inflation sacrifice ratio","FontSize",label_font_size);
xlabel("s.s. CBDC to GDP ratio","FontSize",label_font_size);
yline(nocbdc_sacrifice_ratio,"--r","LineWidth",plot_line_width);
xlim([cbdc_to_gdp_ratios(1) cbdc_to_gdp_ratios(length(cbdc_to_gdp_ratios))]);
legend("CBDC - Quant. rule","No CBDC",'Orientation','horizontal',...
	'Location',"southeast",'FontSize',lege_font_size);
legend boxon;grid on;
set(gcf,"PaperPosition",(1/1.3)*[0 0 6 5]);
% saveas(gcf,strcat("quant_rule_sr_m_dc_to_gdp",".eps"),"epsc");

save results_sacrifice_ratios.mat cbdc_pr_sacrifice_ratio_i_sp...
cbdc_pr_sacrifice_ratio_T cbdc_qr_sacrifice_ratio_T nocbdc_sacrifice_ratio...
spreads Ts cbdc_to_gdp_ratios

