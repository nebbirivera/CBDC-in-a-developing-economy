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

function [ys,params,check] = cbdc_price_rule_steadystate(ys,exo,M_,options_)

    % read out parameters to access them with their name
    NumberOfParameters = M_.param_nbr;
    for ii = 1:NumberOfParameters
        paramname = M_.param_names{ii};
        eval([ paramname ' = M_.params(' int2str(ii) ');']);
    end
    % initialize indicator
    check = 0;
    
    % Known closed form steady states
    pi_t = pi_obj; 
    i_t = i_ss;
    i_dc_ss = i_t/(i_ss/(i_ss-i_spread));
    i_dc_t = i_dc_ss;
    pi_opt_t = (((1+pi_t)^(1-eps_p)-phi_p)/(1-phi_p))^(1/(1-eps_p))-1;
    aux_mc = (((1+pi_t)*(1-phi_p*(1+pi_t)^(eps_p-1)*betta))/(1-phi_p*...
        (1+pi_t)^(eps_p)*betta))^(-1);
    mc_t = ((eps_p-1)/eps_p)*(1+pi_opt_t)*aux_mc;
    ups_p_t = ((1-phi_p)*(1+pi_opt_t)^(-eps_p)*(1+pi_t)^(eps_p))/...
        (1-phi_p*(1+pi_t)^(eps_p));
    SDF_t = betta;
    a_I_t = a_I;
    Ac_t = A;
    AInv_t = A;
    z_L_t = 1;
    z_u_t = 1;
    S_TC_t = 1;
    S_TCc_t = 1;
    S_TCInv_t = 1;
    tau_w_t = tau_w_ss;
    tau_c_t = tau_c_ss;
    tau_r_t = tau_r_ss;
    g_t = g_ss;

    % Labor part steady state targets
    % This new ss targets generate the same calibrated parameters 
    % as in the nocbdc variant
    load cbdc_pr_labor_rates_ss
    lab_force_part_rate = cbdc_pr_lfpr_ss;
    inf_rate = cbdc_pr_infor_rate_ss;
    un_rate = cbdc_pr_unemp_rate_ss;
    lab_force_part_rate_t = lab_force_part_rate;
    un_t = un_rate*lab_force_part_rate;
    L_t = lab_force_part_rate - un_t;
    L_t + un_t - lab_force_part_rate_t;
    L_I_t = inf_rate*L_t;
    L_F_t = L_t - L_I_t;
    a_F_t = 1;
    a_t = 1; 
    y_I_t = a_I*L_I_t;
    y_F_t = a_F_t*L_F_t;
    y_L_t = (y_F_t^((eps_L-1)/eps_L) + y_I_t^((eps_L-1)/eps_L))^(eps_L/(eps_L-1));
    unemp_rate_t = un_rate;
    infor_rate_t = inf_rate;
    O_t = 1 - (L_t + un_t);

    % Set vector of fixed values for numerical solver
    fixed(1) = A;
    fixed(2) = B;
    fixed(3) = tau_c_t;
    fixed(4) = delta;
    fixed(5) = theta;
    fixed(6) = i_t;
    fixed(7) = i_dc_t;
    fixed(8) = ups_p_t;
    fixed(9) = betta;
    fixed(10) = alpha_k;
    fixed(11) = mc_t;
    fixed(12) = a_t;
    fixed(13) = y_L_t;
    fixed(14) = tau_r_t;
    fixed(15) = un_t;
    fixed(16) = g_t;
    fixed(17) = xi_F;
    fixed(18) = eta_F;
    fixed(19) = xi_I;
    fixed(20) = eta_I;
    fixed(21) = omega;
    fixed(22) = T;
    fixed(23) = L_F_t;
    fixed(24) = L_I_t;
    
    % Start parameters for refining calibration and steady state computation
    flag_accuracy = false;
    last_step_size = 1;
    xi_I_prev = xi_I;
    xi_F_prev = xi_F;
    zzeta_prev = zzeta;
    iter_count = 0;
    min_last_step_size = 1e-4;
    max_num_iter = 100;

    while flag_accuracy ~= true

        if abs(last_step_size) < min_last_step_size
            % If the difference between the current value and the last
            % iteration is greater less than the tolerance value then
            % the routine stops the refining of the calibration
            xi_F = xi_F_prev;
            xi_I = xi_I_prev;
            zzeta = zzeta_prev;
            flag_accuracy = true;
        elseif iter_count > max_num_iter
            % If the number of iterations are grater than a previously defined
            % maximum value the iteration is halted
            display('Max iter count reached');
            flag_accuracy = true;
        elseif abs(last_step_size) >= min_last_step_size
            iter_count = iter_count+1;
            fixed(17) = xi_F_prev;
            fixed(19) = xi_I_prev;
            zzeta = zzeta_prev;

            % Compute initial guess values for numeric solver
            [c_0 k_0 m_dc_0] = cbdc_price_rule_initial_guess(fixed, ones(3,1));
            x0 = [double(c_0) double(k_0) double(m_dc_0)];

            % Set numeric solver options
            options = optimset('Display', 'off', 'TolX', 1e-8, 'TolFun',...
                1e-8, 'MaxFunEvals',1e4,'MaxIter',1e4);
            [x1 Fval] = fsolve(@(x) cbdc_price_rule_num_eqs(x, fixed), x0,...
                options);

        end
        % Assign results from numeric solver
        c_t = x1(1);
        k_t = x1(2);
        m_dc_t = x1(3);
        Inv_t = delta*k_t;
        Y_t = (a_t/ups_p_t)*(k_t)^(alpha_k)*(y_L_t)^(1-alpha_k);
        m_c_t = m_dc_t*((i_t-i_dc_t)/(T^(theta)*i_t))^(1/(1-theta));
        LGF_t = (m_c_t)^(theta) + (T*m_dc_t)^(theta);
        LGF_m_c_t = theta*(m_c_t)^(theta-1);
        LGF_m_dc_t = theta*(m_dc_t)^(theta-1)*T^(theta);
        vc_t = (1+tau_c_t)*c_t/LGF_t;
        vInv_t = Inv_t/LGF_t;
        sc_t = Ac_t*vc_t + B/vc_t - 2*sqrt(A*B);
        sInv_t = AInv_t*vInv_t + B/vInv_t - 2*sqrt(A*B);
        vac_F_t = (L_F_t*eta_F/(un_t^omega))^(1/(1-omega));
        vac_I_t = (L_I_t*eta_I/(un_t^omega))^(1/(1-omega));
        q_F_t = (un_t/vac_F_t)^(omega);
        q_I_t = (un_t/vac_I_t)^(omega);
        pr_F_t = (un_t/vac_F_t)^(omega-1);
        pr_I_t = (un_t/vac_I_t)^(omega-1);
        X_t = c_t;
        u_c_t = z_u_t/(c_t - pssi*z_L_t*X_t*((L_t)^(1+chi))/(1+chi));
        lambda_X_t = -(pssi*(L_t)^(1+chi)*u_c_t)/((1+chi)*(1-betta*(1-gama)));
        lambda_c_t = (u_c_t + gama*lambda_X_t)/((1+tau_c_t)*(1+2*(Ac_t*vc_t-...
            sqrt(A*B))));
        p_L_t = (1-alpha_k)*mc_t*Y_t/y_L_t;
        p_F_t = p_L_t*(y_L_t/y_F_t)^(1/eps_L);
        p_I_t = p_L_t*(y_L_t/y_I_t)^(1/eps_L);
        V_I_f_t = (muu/(1-(1-eta_I)*betta))*(p_I_t*a_I_t - pssi*c_t*...
            L_t^(chi)*u_c_t/lambda_c_t);
        w_I_t = p_I_t*a_I_t - V_I_f_t*(1-(1-eta_I)*betta);
        V_I_h_t = (w_I_t - pssi*c_t*L_t^(chi)*u_c_t/lambda_c_t)/...
            (1-(1-eta_I)*betta);
        lambda_I_t = V_I_h_t*lambda_c_t;
        V_F_f_t = muu*(p_F_t*a_F_t - (1+tau_w_t)*pssi*c_t*L_t^(chi)*u_c_t/...
            lambda_c_t)/((1+tau_w_t*(1-muu))*(1-(1-eta_F)*betta));
        w_F_t = (p_F_t*a_F_t - (1-(1-eta_F)*betta)*V_F_f_t)/(1+tau_w_t);
        V_F_h_t = (w_F_t - pssi*c_t*L_t^(chi)*u_c_t/lambda_c_t)/...
            (1-(1-eta_F)*betta);

        xi_I = V_I_f_t*q_I_t;
        xi_F = V_F_f_t*q_F_t;
        lambda_Inv_t = lambda_c_t*(1+2*(AInv_t*vInv_t-sqrt(A*B)));
        lambda_F_t = V_F_h_t*lambda_c_t;
        lambda_I_t = V_I_h_t*lambda_c_t;
        zzeta = (pr_F_t*lambda_F_t + pr_I_t*lambda_I_t)/un_t;
        last_step_size = abs(xi_I - xi_I_prev)+abs(xi_F - xi_F_prev)+...
            abs(zzeta - zzeta_prev);
        xi_F_prev = xi_F;
        xi_I_prev = xi_I;
        zzeta_prev = zzeta;

    end

    aux_p_1_t = (mc_t*Y_t)/(1-phi_p*SDF_t*(1+pi_t)^(eps_p));
    aux_p_2_t = Y_t/(1-phi_p*SDF_t*(1+pi_t)^(eps_p-1));
    r_t = (1+2*(A*(Inv_t/(LGF_t))-sqrt(A*B)))*((1/betta)-1+delta)/(1-tau_r_t);
    gdp_t = (c_t+Inv_t)/(1-g_t);
    G_t = g_t*gdp_t;
    d_ss = 0;
    d_t = d_ss;
    tau_ls_ss = (1+pi_t)^(-1)*((i_t-pi_t)*d_t+(-pi_t)*m_c_t+(i_dc_t-pi_t)*...
        m_dc_t) + G_t - (tau_c_t*c_t + tau_r_t*r_t*k_t + tau_w_t*w_F_t*L_F_t);
    tau_ls_t = tau_ls_ss;
    vac_t = vac_F_t + vac_I_t;
    match_t = (un_t)^(omega)*(vac_F_t)^(1-omega)+(un_t)^(omega)*...
        (vac_I_t)^(1-omega);
    TFP_t = Y_t/((k_t)^(alpha_k)*(L_t)^(1-alpha_k));
    tc_avg_t = ((1+2*(Ac_t*vc_t-sqrt(A*B)))*c_t+...
        (1+2*(AInv_t*vInv_t-sqrt(A*B)))*Inv_t)/(c_t+Inv_t);
    w_avg_t = w_F_t*(1-infor_rate_t)+w_I_t*infor_rate_t;

    c_t = log(c_t);
    Inv_t = log(Inv_t);
    X_t = log(X_t);
    k_t = log(k_t);
    L_F_t = log(L_F_t);
    L_I_t = log(L_I_t);
    un_t = log(un_t);
    O_t = log(O_t);
    m_c_t = log(m_c_t);
    m_dc_t = log(m_dc_t);
    w_F_t = log(w_F_t);
    w_I_t = log(w_I_t);
    lambda_c_t = log(lambda_c_t);
    lambda_Inv_t = log(lambda_Inv_t);
    lambda_F_t = log(lambda_F_t);
    lambda_I_t = log(lambda_I_t);
    u_c_t = log(u_c_t);
    sc_t = log(sc_t);
    sInv_t = log(sInv_t);
    vc_t = log(vc_t);
    vInv_t = log(vInv_t);
    LGF_t = log(LGF_t);
    LGF_m_c_t = log(LGF_m_c_t);
    LGF_m_dc_t = log(LGF_m_dc_t);
    V_F_h_t = log(V_F_h_t);
    V_F_f_t = log(V_F_f_t);
    V_I_h_t = log(V_I_h_t);
    V_I_f_t = log(V_I_f_t);
    vac_F_t = log(vac_F_t);
    vac_I_t = log(vac_I_t);
    q_F_t = log(q_F_t);
    q_I_t = log(q_I_t);
    pr_F_t = log(pr_F_t);
    pr_I_t = log(pr_I_t);
    Y_t = log(Y_t);
    p_L_t = log(p_L_t);
    p_F_t = log(p_F_t);
    p_I_t = log(p_I_t);
    y_L_t = log(y_L_t);
    mc_t = log(mc_t);
    SDF_t = log(SDF_t);
    y_F_t = log(y_F_t);
    y_I_t = log(y_I_t);
    ups_p_t = log(ups_p_t);
    aux_p_1_t = log(aux_p_1_t);
    aux_p_2_t = log(aux_p_2_t);
    L_t = log(L_t);
    gdp_t = log(gdp_t);
    gdp_ss = gdp_t;
    G_t = log(G_t);
    c_ss = c_t;
    k_ss = k_t;
    r_ss = r_t;
    w_F_ss = w_F_t;
    L_F_ss = L_F_t;
    m_dc_ss = m_dc_t;
    m_c_ss = m_c_t; 
    vac_t = log(vac_t);
    match_t = log(match_t);
    TFP_t = log(TFP_t);
    tc_avg_t = log(tc_avg_t);
    w_avg_t = log(w_avg_t);

    %% end own model equations

    params=NaN(NumberOfParameters,1);
    for iter = 1:length(M_.params) 
        eval([ 'params(' num2str(iter) ') = ' M_.param_names{iter} ';' ])
    end

    NumberOfEndogenousVariables = M_.orig_endo_nbr; 
    for ii = 1:NumberOfEndogenousVariables
        varname = M_.endo_names{ii};
        eval(['ys(' int2str(ii) ') = ' varname ';']);
    end 
