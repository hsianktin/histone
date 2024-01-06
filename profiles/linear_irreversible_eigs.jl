k_ons = [1.0]
E_bs = [log(10) * i for i in -4:0.5:-1]
E_Ls = [0.0]
E_Ps = [0.0]
E_P_H3_H4s = [2.0]
k_undockings = [1.0]
label = "irreversible_eigs"
pars = [
    Par(k_on = k_on, E_c = E_c, E_L = E_L, E_P_H2d = E_P_H2d, E_P_H3_H4 = E_P_H2d, k_undocking_simple = k_undocking, k_undocking_coupled = k_undocking)
    for k_on in k_ons, E_c in E_bs, E_L in E_Ls, E_P_H2d in E_Ps,  k_undocking in k_undockings
        if true # reshape high-dimensional arrays to 1D array
    ]