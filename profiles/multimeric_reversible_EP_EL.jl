k_ons = [1.0]
E_bs = [-2.0]
E_Ls = [-1.0, -2, -3, -4, -5]
E_Ps = [i for i in 0:1:15]
E_P_H3_H4s = [2.0]
k_undockings = [0.01]
label = "reversible_EP_EL"
pars = [
    Par(k_on = k_on, E_c = E_c, E_L = E_L, E_P_H2d = E_P_H2d, E_P_H3_H4 = E_P_H2d, k_undocking_simple = k_undocking, k_undocking_coupled = k_undocking)
    for k_on in k_ons, E_c in E_bs, E_L in E_Ls, E_P_H2d in E_Ps,  k_undocking in k_undockings
        if true # reshape high-dimensional arrays to 1D array
    ]