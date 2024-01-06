k_ons = [1.0]
E_bs = [-2.0]
E_Ls = [ε for ε in -10:1:-1]
E_P_H2ds = [2.0]
E_P_H3_H4s = [2.0]
k_undockings = [10.0^i for i in -8:3:2]
label = "reversible_EL"
pars = [
    Par(k_on = k_on, E_c = E_c, E_L = E_L, E_P_H2d = E_P_H2d, E_P_H3_H4 = E_P_H3_H4, k_undocking_simple = k_undocking, k_undocking_coupled = k_undocking)
    for k_on in k_ons, E_c in E_bs, E_L in E_Ls, E_P_H2d in E_P_H2ds, E_P_H3_H4 in E_P_H3_H4s, k_undocking in k_undockings
        if true # reshape high-dimensional arrays to 1D array
    ]