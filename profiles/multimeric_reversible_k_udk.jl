k_ons = [1.0]
E_bs = [-2.0]
E_Ls = [-1.0]
E_Ps = [2.0]
E_P_H3_H4s = [2.0]
ratios = [1,0.1,0.01,0.001,0.0001]
k_undocking_coupleds = [10.0^i for i in -12:1:0]
label = "reversible_k_udk"
pars = [
    Par(k_on = k_on, E_c = E_c, E_L = E_L, E_P_H2d = E_P_H2d, E_P_H3_H4 = E_P_H2d, k_undocking_simple = k_undocking * ratio, k_undocking_coupled = k_undocking)
    for k_on in k_ons, E_c in E_bs, E_L in E_Ls, E_P_H2d in E_Ps,  k_undocking in k_undocking_coupleds, ratio in ratios
        if true # reshape high-dimensional arrays to 1D array
    ]