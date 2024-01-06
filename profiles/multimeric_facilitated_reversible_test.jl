k_ons = [1.0]
E_bs = [-2.0]
E_Ls = [-1.0]
E_P_H2ds = [Inf]
E_Ps = [i for i in 3:3] #∪ [Inf]
k_undockings = [0.01]
q_ons = [10.0^i for i in -4:0.25:2]
q_offs = [0.01]

pars = [
    Par(k_on = k_on, E_c = E_c, E_L = E_L, E_P_H2d = E_P, E_P_H3_H4 = E_P, k_undocking_simple = k_undocking, k_undocking_coupled = k_undocking, q_on = q_on, q_off = q_off)
    for k_on in k_ons, E_c in E_bs, E_P in E_Ps, k_undocking in k_undockings, q_on in q_ons, q_off in q_offs, E_L in E_Ls
        if true
]

label = "reversible_test"