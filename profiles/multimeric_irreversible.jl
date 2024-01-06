k_ons = [1.0]
E_bs = [-2.0]
E_Ls = [-20, -15, -10, -5, 0]
E_P_H2ds = [Inf]
E_P_H3_H4s = [Inf]
k_undocking_simples = [0.001]
k_undocking_coupleds = [10.0^i for i in -8:1:-1]
k_undockings = [10.0^i for i in -13:1:-1]
label = "irreversible"
pars = [
    Par(k_on = k_on, E_c = E_c, E_L = E_L, E_P_H2d = E_P_H2d, E_P_H3_H4 = E_P_H3_H4, k_undocking_simple = k_undocking, k_undocking_coupled = k_undocking)
    for k_on in k_ons, E_c in E_bs, E_L in E_Ls, E_P_H2d in E_P_H2ds, E_P_H3_H4 in E_P_H3_H4s, k_undocking in k_undockings
        if true # reshape high-dimensional arrays to 1D array
    ]