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


% Computes and plots IRFs at calibrated size of shocks 
% clear;close;clc;
clear;


%%%% IRFs analysis
% Set up variables controlling aspects of the plots
irf_line_width = 1.5;
zero_line_width = 0.8;
title_font_size = 11;
prop_factor = 100;

nocbdc_irf_color = "red";
cbdc_irf_pr_color = "blue";
cbdc_irf_qr_color = "green";
nocbdc_irf_ls = "-";
cbdc_irf_pr_ls = "--";
cbdc_irf_qr_ls = "-.";

x_label = "Quarter";
y_label = "% s.s. dev.";
x_label_font_size = 9;
y_label_font_size = 9;

%% Run the three model variants and get irfs
% dynare('nocbdc.mod','nowarn')
% nocbdc_oo = oo_;
% dynare('cbdc_price_rule.mod','nowarn')
% cbdc_price_rule_oo = oo_;
% dynare('cbdc_quant_rule.mod','nowarn')
% cbdc_quant_rule_oo = oo_;
%% Save the results for avoiding re-running the models when testing
% save baseline_irfs_results.mat nocbdc_oo cbdc_price_rule_oo cbdc_quant_rule_oo
load baseline_irfs_results.mat

variable_names_main={"gdp_t","pi_t","c_t","Inv_t","i_t","unemp_rate_t",...
"infor_rate_t","match_t","LGF_t","tc_avg_t","r_t","lab_force_part_rate_t"};
variable_names_main_titles={"\textbf{GDP}","\textbf{Inflation}",...
"\textbf{Consumption}","\textbf{Investment}",...
"\textbf{Main policy interest}","\textbf{Unemployment rate}",...
"\textbf{Informality rate}","\textbf{Job creation}",...
"\textbf{LGF}",...
"\textbf{Avg. transac. costs distor.}",...
"\textbf{Real interest}","\textbf{LFPR}"};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% IRFs of main variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close all;
% shock_names = {"_e_a_t","_e_i_t","_e_a_F_t","_e_g_t"};

% for s=1:length(shock_names)
% 	figure;
% 	for k=1:length(variable_names_main)
% 		subplot(3,4,k);
% 		yline(0,":k","LineWidth",zero_line_width);
% 		hold on;
% 		% No-CBDC IRFs
% 		plot(1:20,prop_factor*nocbdc_oo.irfs.(...
% 			strcat(char(variable_names_main{k}),...
% 			char(shock_names{s}))),...
% 			nocbdc_irf_ls,"LineWidth",irf_line_width,"Color",nocbdc_irf_color);
% 		xlim([1 20]);
% 		hold on;
% 		% Price rule IRFs
% 		plot(1:20,prop_factor*cbdc_price_rule_oo.irfs.(strcat(...
% 			char(variable_names_main{k}),char(shock_names{s}))),...
% 			cbdc_irf_pr_ls,"LineWidth",irf_line_width,"Color",...
% 			cbdc_irf_pr_color);
% 		xlim([1 20]);
% 		hold on;
% 		% Quant. rule IRFs
% 		plot(1:20,prop_factor*cbdc_quant_rule_oo.irfs.(strcat(...
% 			char(variable_names_main{k}),char(shock_names{s}))),...
% 			cbdc_irf_qr_ls,"LineWidth",irf_line_width,"Color",...
% 			cbdc_irf_qr_color);
% 		xlim([1 20]);
% 		% Figure titles
% 		title(variable_names_main_titles(k)," ","Interpreter","latex",...
% 			"fontsize",title_font_size);%,'Position');, [0],...
% 			% "Units","normalized");
% 		xlabel(x_label,"FontSize",x_label_font_size);
% 		ylabel(y_label,"FontSize",y_label_font_size);
% 		grid on;box on;hold off;
% 	end
% 	legen{1,1}="";legen{1,2}="No CBDC";legen{1,3}="CBDC - Price rule";
% 	legen{1,4}="CBDC - Quant. rule";
% 	lege=legend(legen,'Orientation','horizontal',...
% 		'Position',[0.5 0.04 0 0],'FontSize',10);
% 	set(lege,'Box','on');
% 	set(gcf,"PaperPosition",[0 0 10.5 7.75]);
% 	% saveas(gcf,strcat("irfs",shock_names{s},".eps"),"epsc");
% end
% % close all;


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%% IRFs of other variables
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% variable_names_2nd={"m_c_t","m_dc_t","w_I_t","w_F_t","w_avg_t","V_I_h_t",...
% "V_F_h_t","V_I_f_t","V_F_f_t","i_dc_t","mc_t","vac_t"};
% variable_names_2nd_titles={"\textbf{Cash}","\textbf{CBDC}",...
% "\textbf{Informal wage}","\textbf{Formal wage}","\textbf{Average wage}",...
% "\textbf{V. acc. Inf. HH.}","\textbf{V. acc. For. HH.}",...
% "\textbf{V. acc. Inf. firm}","\textbf{V. acc. For. firm}",...
% "\textbf{CBDC interest}","\textbf{Real marginal cost}",...
% "\textbf{Aggregate vacancies}"};
% shock_names = {"_e_a_t","_e_i_t","_e_a_F_t","_e_g_t"};
% for s=1:length(shock_names)
% 	figure;
% 	for k=1:length(variable_names_2nd)
% 		subplot(3,4,k);
% 		yline(0,":k","LineWidth",zero_line_width);
% 		hold on;
% 		% No-CBDC IRFs
% 		if variable_names_2nd{k}~="m_dc_t" && variable_names_2nd{k}~="i_dc_t"
% 			plot(1:20,prop_factor*nocbdc_oo.irfs.(...
% 				strcat(char(variable_names_2nd{k}),...
% 				char(shock_names{s}))),...
% 				nocbdc_irf_ls,"LineWidth",irf_line_width,"Color",nocbdc_irf_color);
% 			xlim([1 20]);
% 			hold on;
% 		end
% 		% Price rule IRFs
% 		plot(1:20,prop_factor*cbdc_price_rule_oo.irfs.(strcat(...
% 			char(variable_names_2nd{k}),char(shock_names{s}))),...
% 			cbdc_irf_pr_ls,"LineWidth",irf_line_width,"Color",...
% 			cbdc_irf_pr_color);
% 		xlim([1 20]);
% 		hold on;
% 		% Quant. rule IRFs
% 		plot(1:20,prop_factor*cbdc_quant_rule_oo.irfs.(strcat(...
% 			char(variable_names_2nd{k}),char(shock_names{s}))),...
% 			cbdc_irf_qr_ls,"LineWidth",irf_line_width,"Color",...
% 			cbdc_irf_qr_color);
% 		xlim([1 20]);
% 		% Figure titles
% 		title(variable_names_2nd_titles(k)," ","Interpreter","latex",...
% 			"fontsize",title_font_size);%,'Position');, [0],...
% 			% "Units","normalized");
% 		xlabel(x_label,"FontSize",x_label_font_size);
% 		ylabel(y_label,"FontSize",y_label_font_size);
% 		grid on;box on;hold off;
% 	end
% 	legen{1,1}="";legen{1,2}="No CBDC";legen{1,3}="CBDC - Price rule";
% 	legen{1,4}="CBDC - Quant. rule";
% 	lege=legend(legen,'Orientation','horizontal',...
% 		'Position',[0.5 0.05 0 0],'FontSize',10);
% 	set(lege,'Box','on');
% 	set(gcf,"PaperPosition",[0 0 10.5 7.75]);
% 	saveas(gcf,strcat("irfs",shock_names{s},"_remain",".eps"),"epsc")
% end
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% IRFs to e_a_t with high T
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
disp("Cumulative IRFs to e_a_t changing T and i_spread");
disp("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
dynare('cbdc_price_rule.mod','nowarn')
set_param_value("T",5)
[info, cbdc_price_rule_oo_high_T, options_] = ...
stoch_simul(M_, options_, oo_, var_list_);
set_param_value("T",1.5)
set_param_value("i_spread",0.0015)
[info, cbdc_price_rule_oo_low_i_sp, options_] = ...
stoch_simul(M_, options_, oo_, var_list_);
set_param_value("i_spread",0.015)
% figure;
% for k=1:length(variable_names)
% 	subplot(4,3,k);
% 	yline(0,":k","LineWidth",zero_line_width);
% 	hold on;
% 	plot(1:20,prop_factor*nocbdc_oo.irfs.(...
% 		strcat(char(variable_names{k}),...
% 		char("_e_a_t"))),...
% 		nocbdc_irf_ls,"LineWidth",irf_line_width,"Color",nocbdc_irf_color);
% 	xlim([1 20]);
% 	hold on;
% 	plot(1:20,prop_factor*cbdc_price_rule_oo_high_T.irfs.(strcat(...
% 		char(variable_names{k}),char("_e_a_t"))),...
% 		cbdc_irf_pr_ls,"LineWidth",irf_line_width,"Color",...
% 		cbdc_irf_pr_color);
% 	xlim([1 20]);
% 	title(variable_names_titles(k)," ","Interpreter","latex",...
% 		"fontsize",title_font_size);%,'Position');, [0],...
% 		% "Units","normalized");
% 	xlabel(x_label,"FontSize",x_label_font_size);
% 	ylabel(y_label,"FontSize",y_label_font_size);
% 	grid on;
% 	hold off;
% end
% legen{1,1}="";legen{1,2}="No CBDC";legen{1,3}="CBDC - Price rule";
% lege=legend(legen,'Orientation','horizontal',...
% 	'Position',[0.5 0.05 0 0],'FontSize',10);
% set(lege,'Box','on');
% set(gcf,"PaperPosition",[0 0 8 8]);
% saveas(gcf,strcat("irfs","_e_a_t","_high_T",".eps"),"epsc")

% Create table with relative cumulative effects of shock to a_t
% No CBDC
cumul_gdp_nocbdc = cumsum(nocbdc_oo.irfs.gdp_t_e_a_t);
cumul_inflation_nocbdc = cumsum(nocbdc_oo.irfs.pi_t_e_a_t);
cumul_unemp_rate_nocbdc = cumsum(nocbdc_oo.irfs.unemp_rate_t_e_a_t);
cumul_a_t_nocbdc = cumsum(nocbdc_oo.irfs.a_t_e_a_t);
disp("No-CBDC gdp_t/a_t cumulative response")
cumul_gdp_nocbdc(4)/cumul_a_t_nocbdc(4)
disp("No-CBDC pi_t/a_t cumulative response")
cumul_inflation_nocbdc(4)/cumul_a_t_nocbdc(4)
disp("No-CBDC unemp_rate_t/a_t cumulative response")
cumul_unemp_rate_nocbdc(4)/cumul_a_t_nocbdc(4)

% Price rule with baseline T=1.5 and i_spread=1.5%
cumul_gdp_cbdc_base_T = cumsum(cbdc_price_rule_oo.irfs.gdp_t_e_a_t);
cumul_inflation_cbdc_base_T = cumsum(cbdc_price_rule_oo.irfs.pi_t_e_a_t);
cumul_unemp_rate_cbdc_base_T = ...
cumsum(cbdc_price_rule_oo.irfs.unemp_rate_t_e_a_t);
cumul_a_t_cbdc_base_T = cumsum(cbdc_price_rule_oo.irfs.a_t_e_a_t);
disp("CBDC-pr gdp_t/a_t cumulative response T=1.5")
cumul_gdp_cbdc_base_T(4)/cumul_a_t_cbdc_base_T(4)
disp("CBDC-pr pi_t/a_t cumulative response T=1.5")
cumul_inflation_cbdc_base_T(4)/cumul_a_t_cbdc_base_T(4)
disp("CBDC-pr unemp_rate_t/a_t cumulative response T=1.5")
cumul_unemp_rate_cbdc_base_T(4)/cumul_a_t_cbdc_base_T(4)

% Price rule with T=5, i_spread=1.5%
cumul_gdp_cbdc_high_T = cumsum(cbdc_price_rule_oo_high_T.irfs.gdp_t_e_a_t);
cumul_inflation_cbdc_high_T = cumsum(cbdc_price_rule_oo_high_T.irfs.pi_t_e_a_t);
cumul_unemp_rate_cbdc_high_T = ...
cumsum(cbdc_price_rule_oo_high_T.irfs.unemp_rate_t_e_a_t);
cumul_a_t_cbdc_high_T = cumsum(cbdc_price_rule_oo_high_T.irfs.a_t_e_a_t);
disp("CBDC-pr gdp_t/a_t cumulative response T=5, i_spread=1.5%")
cumul_gdp_cbdc_high_T(4)/cumul_a_t_cbdc_high_T(4)
disp("CBDC-pr pi_t/a_t cumulative response T=5, i_spread=1.5%")
cumul_inflation_cbdc_high_T(4)/cumul_a_t_cbdc_high_T(4)
disp("CBDC-pr gdp_t/a_t cumulative response T=5, i_spread=1.5%")
cumul_unemp_rate_cbdc_high_T(4)/cumul_a_t_cbdc_high_T(4)

% Price rule with i_spread=0.15%, T=1.5
cumul_gdp_cbdc_low_i_sp = cumsum(cbdc_price_rule_oo_low_i_sp.irfs.gdp_t_e_a_t);
cumul_inflation_cbdc_low_i_sp = ...
cumsum(cbdc_price_rule_oo_low_i_sp.irfs.pi_t_e_a_t);
cumul_unemp_rate_cbdc_low_i_sp = ...
cumsum(cbdc_price_rule_oo_low_i_sp.irfs.unemp_rate_t_e_a_t);
cumul_a_t_cbdc_low_i_sp = cumsum(cbdc_price_rule_oo_low_i_sp.irfs.a_t_e_a_t);
disp("CBDC-pr gdp_t/a_t cumulative response T=1.5, i_spread=0.15%")
cumul_gdp_cbdc_low_i_sp(4)/cumul_a_t_cbdc_low_i_sp(4)
disp("CBDC-pr pi_t/a_t cumulative response T=1.5, i_spread=0.15%")
cumul_inflation_cbdc_low_i_sp(4)/cumul_a_t_cbdc_low_i_sp(4)
disp("CBDC-pr gdp_t/a_t cumulative response T=1.5, i_spread=0.15%")
cumul_unemp_rate_cbdc_low_i_sp(4)/cumul_a_t_cbdc_low_i_sp(4)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Cumulative IRFs to e_g_t changing T and i_spread
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
disp("Cumulative IRFs to e_g_t changing T and i_spread");
disp("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
variable_names_sens={"gdp_t","pi_t","c_t","Inv_t","i_t","unemp_rate_t",...
"infor_rate_t","match_t","LGF_t","tc_avg_t","r_t","lab_force_part_rate_t",...
"i_dc_t","m_c_t","m_dc_t","w_avg_t","vac_t","O_t","vInv_t","vc_t",...
"lambda_F_t","lambda_I_t","LGF_m_c_t","LGF_m_dc_t","L_t"};
variable_names_sens_titles={"\textbf{GDP}","\textbf{Inflation}",...
"\textbf{Consumption}","\textbf{Investment}",...
"\textbf{Main policy interest}","\textbf{Unemployment rate}",...
"\textbf{Informality rate}","\textbf{Job creation}",...
"\textbf{LGF}",...
"\textbf{Avg. transac. costs distor.}",...
"\textbf{Real interest}","\textbf{LFPR}",...
"\textbf{CBDC interest}","\textbf{Cash}","\textbf{CBDC}",...
"\textbf{Avg. wage}","\textbf{Agg. vacancies}",...
"\textbf{Out of the lf}","\textbf{Inv. vel.}","\textbf{Con. vel.}",...
"\textbf{MU Inf. lab}","\textbf{MU For. lab}","\textbf{ML $m^{c}$}",...
"\textbf{ML $m^{dc}$}","Labor supply"};
shock_names_sens = {"_e_a_t","_e_g_t","_e_i_t"};
for s=1:length(shock_names_sens)
	figure;
	for k=1:length(variable_names_sens)
		subplot(5,5,k);
		yline(0,":k","LineWidth",zero_line_width);
		hold on;
		if variable_names_sens{k}~="i_dc_t"&&variable_names_sens{k}~="m_dc_t"&&...
			variable_names_sens{k}~="LGF_m_dc_t"
			plot(1:20,prop_factor*nocbdc_oo.irfs.(...
				strcat(char(variable_names_sens{k}),...
				char(shock_names_sens{s}))),...
				nocbdc_irf_ls,"LineWidth",irf_line_width,"Color",nocbdc_irf_color);
			xlim([1 20]);hold on;
		end
		plot(1:20,prop_factor*cbdc_price_rule_oo.irfs.(...
			strcat(char(variable_names_sens{k}),...
			char(shock_names_sens{s}))),...
			cbdc_irf_pr_ls,"LineWidth",irf_line_width,"Color",cbdc_irf_pr_color);
		xlim([1 20]);hold on;
		plot(1:20,prop_factor*cbdc_price_rule_oo_high_T.irfs.(strcat(...
			char(variable_names_sens{k}),char(shock_names_sens{s}))),...
			"-.","LineWidth",irf_line_width,"Color","#77AC30");
		xlim([1 20]);hold on;
		plot(1:20,prop_factor*cbdc_price_rule_oo_low_i_sp.irfs.(strcat(...
			char(variable_names_sens{k}),char(shock_names_sens{s}))),...
			"-o","LineWidth",irf_line_width,"Color","#D95319");
		xlim([1 20]);
		title(variable_names_sens_titles(k)," ","Interpreter","latex",...
			"fontsize",title_font_size);
		xlabel(x_label,"FontSize",x_label_font_size);
		ylabel(y_label,"FontSize",y_label_font_size);
		grid on;box on;hold off;
	end
	% suptitle("Response senti")
	legen{1,1}="";legen{1,2}="No CBDC";legen{1,3}="CBDC - Price rule";
	lege=legend({"","No CBDC","CBDC-Price rule base","CBDC-Price rule T=5",...
		"CBDC-Price rule spread=0.15%"},...
		'Orientation','horizontal',...
		'Position',[0.5 0.05 0 0],'FontSize',10);
	set(lege,'Box','on');
	set(gcf,"PaperPosition",[0 0 15 13]);
	% saveas(gcf,strcat("irfs_sens",shock_names_sens{s},".eps"),"epsc");
end
close all;

% Create table with relative cumulative effects of shock to g_t
% No CBDC
cumul_gdp_nocbdc = cumsum(nocbdc_oo.irfs.gdp_t_e_g_t);
cumul_inflation_nocbdc = cumsum(nocbdc_oo.irfs.pi_t_e_g_t);
cumul_unemp_rate_nocbdc = cumsum(nocbdc_oo.irfs.unemp_rate_t_e_g_t);
cumul_consump_nocbdc = cumsum(nocbdc_oo.irfs.c_t_e_g_t);
cumul_Inv_nocbdc = cumsum(nocbdc_oo.irfs.Inv_t_e_g_t);
cumul_g_t_nocbdc = cumsum(nocbdc_oo.irfs.g_t_e_g_t);
disp("No-CBDC gdp_t/g_t cumulative response")
cumul_gdp_nocbdc(4)/cumul_g_t_nocbdc(4)
disp("No-CBDC pi_t/g_t cumulative response")
cumul_inflation_nocbdc(4)/cumul_g_t_nocbdc(4)
disp("No-CBDC unemp_rate_t/g_t cumulative response")
cumul_unemp_rate_nocbdc(4)/cumul_g_t_nocbdc(4)
disp("No-CBDC c_t/g_t cumulative response")
cumul_consump_nocbdc(4)/cumul_g_t_nocbdc(4)
disp("No-CBDC Inv_t/g_t cumulative response")
cumul_Inv_nocbdc(4)/cumul_g_t_nocbdc(4)

% Price rule with baseline T=1.5 and i_spread=1.5%
cumul_gdp_cbdc_base_T = cumsum(cbdc_price_rule_oo.irfs.gdp_t_e_g_t);
cumul_inflation_cbdc_base_T = cumsum(cbdc_price_rule_oo.irfs.pi_t_e_g_t);
cumul_unemp_rate_cbdc_base_T = ...
cumsum(cbdc_price_rule_oo.irfs.unemp_rate_t_e_g_t);
cumul_consump_base_T = cumsum(cbdc_price_rule_oo.irfs.c_t_e_g_t);
cumul_Inv_cbdc_base_T = cumsum(cbdc_price_rule_oo.irfs.Inv_t_e_g_t);
cumul_g_t_cbdc_base_T = cumsum(cbdc_price_rule_oo.irfs.g_t_e_g_t);
disp("CBDC-pr gdp_t/g_t cumulative response T=1.5")
cumul_gdp_cbdc_base_T(4)/cumul_g_t_cbdc_base_T(4)
disp("CBDC-pr pi_t/g_t cumulative response T=1.5")
cumul_inflation_cbdc_base_T(4)/cumul_g_t_cbdc_base_T(4)
disp("CBDC-pr unemp_rate_t/g_t cumulative response T=1.5")
cumul_unemp_rate_cbdc_base_T(4)/cumul_g_t_cbdc_base_T(4)
disp("CBDC-pr c_t/g_t cumulative response T=1.5")
cumul_consump_base_T(4)/cumul_g_t_cbdc_base_T(4)
disp("CBDC-pr Inv_t/g_t cumulative response T=1.5")
cumul_Inv_cbdc_base_T(4)/cumul_g_t_cbdc_base_T(4)

% Price rule with T=5, i_spread=1.5%
cumul_gdp_cbdc_high_T = cumsum(cbdc_price_rule_oo_high_T.irfs.gdp_t_e_g_t);
cumul_inflation_cbdc_high_T = cumsum(cbdc_price_rule_oo_high_T.irfs.pi_t_e_g_t);
cumul_unemp_rate_cbdc_high_T = ...
cumsum(cbdc_price_rule_oo_high_T.irfs.unemp_rate_t_e_g_t);
cumul_consump_high_T = cumsum(cbdc_price_rule_oo_high_T.irfs.c_t_e_g_t);
cumul_Inv_cbdc_high_T = cumsum(cbdc_price_rule_oo_high_T.irfs.Inv_t_e_g_t);
cumul_g_t_cbdc_high_T = cumsum(cbdc_price_rule_oo_high_T.irfs.g_t_e_g_t);
disp("CBDC-pr gdp_t/g_t cumulative response T=5, i_spread=1.5%")
cumul_gdp_cbdc_high_T(4)/cumul_g_t_cbdc_high_T(4)
disp("CBDC-pr pi_t/g_t cumulative response T=5, i_spread=1.5%")
cumul_inflation_cbdc_high_T(4)/cumul_g_t_cbdc_high_T(4)
disp("CBDC-pr unemp_rate_t/g_t cumulative response T=5, i_spread=1.5%")
cumul_unemp_rate_cbdc_high_T(4)/cumul_g_t_cbdc_high_T(4)
disp("CBDC-pr c_t/g_t cumulative response T=5, i_spread=1.5%")
cumul_consump_high_T(4)/cumul_g_t_cbdc_high_T(4)
disp("CBDC-pr Inv_t/g_t cumulative response T=5, i_spread=1.5%")
cumul_Inv_cbdc_high_T(4)/cumul_g_t_cbdc_high_T(4)

% Price rule with i_spread=0.15%, T=1.5
cumul_gdp_cbdc_low_i_sp = cumsum(cbdc_price_rule_oo_low_i_sp.irfs.gdp_t_e_g_t);
cumul_inflation_cbdc_low_i_sp = ...
cumsum(cbdc_price_rule_oo_low_i_sp.irfs.pi_t_e_g_t);
cumul_unemp_rate_cbdc_low_i_sp = ...
cumsum(cbdc_price_rule_oo_low_i_sp.irfs.unemp_rate_t_e_g_t);
cumul_consump_low_i_sp = cumsum(cbdc_price_rule_oo_low_i_sp.irfs.c_t_e_g_t);
cumul_Inv_low_i_sp = cumsum(cbdc_price_rule_oo_low_i_sp.irfs.Inv_t_e_g_t);
cumul_g_t_cbdc_low_i_sp = cumsum(cbdc_price_rule_oo_low_i_sp.irfs.g_t_e_g_t);
disp("CBDC-pr gdp_t/g_t cumulative response T=1.5, i_spread=0.15%")
cumul_gdp_cbdc_low_i_sp(4)/cumul_g_t_cbdc_low_i_sp(4)
disp("CBDC-pr pi_t/g_t cumulative response T=1.5, i_spread=0.15%")
cumul_inflation_cbdc_low_i_sp(4)/cumul_g_t_cbdc_low_i_sp(4)
disp("CBDC-pr unemp_rate_t/g_t cumulative response T=1.5, i_spread=0.15%")
cumul_unemp_rate_cbdc_low_i_sp(4)/cumul_g_t_cbdc_low_i_sp(4)
disp("CBDC-pr c_t/g_t cumulative response T=1.5, i_spread=0.15%")
cumul_consump_low_i_sp(4)/cumul_g_t_cbdc_low_i_sp(4)
disp("CBDC-pr Inv_t/g_t cumulative response T=1.5, i_spread=0.15%")
cumul_Inv_low_i_sp(4)/cumul_g_t_cbdc_low_i_sp(4)
