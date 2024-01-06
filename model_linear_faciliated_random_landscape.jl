using SparseArrays
using LinearAlgebra
using ProgressMeter

using DataFrames
using CSV
using Statistics

BLAS.set_num_threads(20)

E_b_landscape(N) = -(rand(N).*0.5 .+ 1.5)
    

function data_acquire(N, k_on, E_bs, q_ons, q_offs)
    NL = 1
    NR = 1
    NM = N - NL - NR
    # only one macro states is allowed, therefore NL, NR values are not used
    MacroStates = [[2,2,1]] 
    ∑ = sum
    # construct the elementary state space [i,j] pair which i + j ≤ N-1
    function Φ(N)
        X = Array{Int64,1}[]
        for i in 0:N-1
            for j in 0:N-1-i
                push!(X,[i,j])
            end
        end
        return X
    end
    # map the macrostates ∈ MacroStates to its corresponding characterization of microstates
    function macro_state2sizes(i,j,k) 
        nc = 0
        nl = 0
        nr = 0
        if i == 0
            # do nothing
        elseif i == 1
            nl += NL
        elseif i == 2
            nc += NL
        else
            println("i domain error")
        end
        if j == 0
            # do nothing
        elseif j == 1
            nr += NR
        elseif j == 2
            nc += NR
        else
            println("j domain error")
        end
        if k == 0
            # do nothing
        elseif k == 1
            nc += NM
        else
            println("k domain error")
        end
        list = [x for x in [nl,nc,nr] if x > 0]
        # macro_state2sizes(2,2,1) = [14] when N = 14
        # macro_state2sizes(1, 2, 1) = [1, 13] when N = 14, NL = 1
        return list
    end
    
    # map the macrostates to its corresponding characterization of microstates
    function macro_state2exact_sizes(i,j,k) 
        nc = 0
        nl = 0
        nr = 0
        if i == 0
            # do nothing
        elseif i == 1
            nl += NL
        elseif i == 2
            nl += NL
        else
            println("i domain error")
        end
        if j == 0
            # do nothing
        elseif j == 1
            nr += NR
        elseif j == 2
            nr += NR
        else
            println("j domain error")
        end
        if k == 0
            # do nothing
        elseif k == 1
            nc += NM
        else
            println("k domain error")
        end
        # macro_state2exact_sizes(2,2,1) = [1,12,1] when N = 14, NL = 1, NR = 1
        # it does not reflect the connection topology between different modules
        return [nl,nc,nr]
    end
    
    function state_mapping(state)
        X,Y,Z = state[1:3]
        sizes = macro_state2sizes(X,Y,Z)
        exact_sizes = macro_state2exact_sizes(X,Y,Z)
        if length(sizes) == 0
            exact_state = zeros(N)
        elseif length(sizes) == 1
            position_flag = [x > 0 for x in exact_sizes]
            undetermined_state = zeros(sizes[1])
            i,j = state[4]
            undetermined_state[i+1:sizes[1]-j] .= 1
            exact_state = Array{Int,1}()
            put_in_flag = true
            for i in 1:3
                flag = position_flag[i]
                if flag == 1 && put_in_flag
                    exact_state = [exact_state;undetermined_state]
                    put_in_flag = false
                end
                if flag == 0
                    exact_state = [exact_state;zeros([NL,NM,NR][i])]
                end
            end
        elseif length(sizes) == 2
            undetermined_state_1 = zeros(sizes[1])
            i1,j1 = state[4]
            undetermined_state_1[i1+1:sizes[1]-j1] .= 1
            undetermined_state_2 = zeros(sizes[2])
            i2,j2 = state[5]
            undetermined_state_2[i2+1:sizes[2]-j2] .= 1
            undetermined_states = Array{Any,1}()
            push!(undetermined_states,undetermined_state_1)
            push!(undetermined_states,undetermined_state_2)
            exact_state = Array{Int,1}()
            cursor = 1
            position_flag = [x > 0 for x in exact_sizes]
            for i in 1:3
                flag = position_flag[i]
                if flag == 1
                    if cursor < 3
                        exact_state = [exact_state; undetermined_states[cursor]]
                        cursor += 1
                    end
                elseif flag == 0
                    exact_state = [exact_state;zeros([NL,NM,NR][i])]
                end
            end
        elseif length(sizes) == 3
            L_state = zeros(NL)
            M_state = zeros(NM)
            R_state = zeros(NR)
            i1,j1 = state[4]
            i2,j2 = state[5]
            i3,j3 = state[6]
            L_state[i1+1:NL-j1] .= 1
            M_state[i2+1:NM-j2] .= 1
            R_state[i3+1:NR-j3] .= 1
            exact_state = [L_state; M_state; R_state]
        end
        return exact_state
    end
    
    # get plausible helicase postion for a given microstate
    function get_plausible_helicase_postion(microstate)
        list = Array{Int,1}[]
        i=j= 0
        while i+1 < N && microstate[i+1] == 0
            i += 1
        end
        # restrict to the case where at least one site is occupied
        while j+1 < N && microstate[N-j] == 0
            j += 1
        end
        for a in 0:1:i, b in 0:1:j
            push!(list,[a,b])
        end
        return list
    end
    
    # the same function as state_mapping, but state ∈ S0 also contains the helicase position
    function state_mapping_with_remodeller(state)
        X,Y,Z = state[1:3]
        sizes = macro_state2sizes(X,Y,Z)
        exact_sizes = macro_state2exact_sizes(X,Y,Z)
        if length(sizes) == 0
            exact_state = zeros(N)
            state = [0]
        elseif length(sizes) == 1
            position_flag = [x > 0 for x in exact_sizes]
            undetermined_state = zeros(sizes[1])
            i,j = state[4]
            undetermined_state[i+1:sizes[1]-j] .= 1
            exact_state = Array{Int,1}()
            put_in_flag = true
            for i in 1:3
                flag = position_flag[i]
                if flag == 1 && put_in_flag
                    exact_state = [exact_state;undetermined_state]
                    put_in_flag = false
                end
                if flag == 0
                    exact_state = [exact_state;zeros([NL,NM,NR][i])]
                end
            end
        elseif length(sizes) == 2
            undetermined_state_1 = zeros(sizes[1])
            i1,j1 = state[4]
            undetermined_state_1[i1+1:sizes[1]-j1] .= 1
            undetermined_state_2 = zeros(sizes[2])
            i2,j2 = state[5]
            undetermined_state_2[i2+1:sizes[2]-j2] .= 1
            undetermined_states = Array{Any,1}()
            push!(undetermined_states,undetermined_state_1)
            push!(undetermined_states,undetermined_state_2)
            exact_state = Array{Int,1}()
            cursor = 1
            position_flag = [x > 0 for x in exact_sizes]
            for i in 1:3
                flag = position_flag[i]
                if flag == 1
                    if cursor < 3
                        exact_state = [exact_state; undetermined_states[cursor]]
                        cursor += 1
                    end
                elseif flag == 0
                    exact_state = [exact_state;zeros([NL,NM,NR][i])]
                end
            end
        elseif length(sizes) == 3
            L_state = zeros(NL)
            M_state = zeros(NM)
            R_state = zeros(NR)
            i1,j1 = state[4]
            i2,j2 = state[5]
            i3,j3 = state[6]
            L_state[i1+1:NL-j1] .= 1
            M_state[i2+1:NM-j2] .= 1
            R_state[i3+1:NR-j3] .= 1
            exact_state = [L_state; M_state; R_state]
        end
        return [exact_state,state[end]]
    end
    
    function is_adjacent(x,y)
        X = S0[x]
        Y = S0[y]
        state_X, helicase_X = state_mapping_with_remodeller(X)
        state_Y, helicase_Y  = state_mapping_with_remodeller(Y)
        if X == [0,0,0] || Y == [0,0,0]
            if ∑(abs.(state_X-state_Y)) > 1
                return false
            else
                return true
            end
        elseif ∑(abs.(state_X-state_Y))+∑(abs.(helicase_X-helicase_Y)) > 1
            return false
        else
            return true
        end
    end
    
    # internal transition from y to x
    function internal_on(x::Int,y::Int)
        X = S0[x]
        Y = S0[y]
        if length(X) == length(Y)
            Δ = X - Y
            if Δ[1:3] == [0,0,0]
            state_X, helicase_X = state_mapping_with_remodeller(X)
            state_Y, helicase_Y  = state_mapping_with_remodeller(Y)
                if ∑(abs.(state_X-state_Y))+∑(abs.(helicase_X-helicase_Y)) == 1
                    if ∑(state_X-state_Y) > 0
                        return 1
                    else
                        return 0
                    end
                else
                    return 0
                end
            else
                return 0
            end
        else 
            return 0
        end
    end
    
    function internal_off(x::Int, y::Int)
        X = S0[x]
        Y = S0[y]
        if length(X) == length(Y)
            Δ = X - Y
            if Δ[1:3] == [0,0,0]
                state_X, helicase_X = state_mapping_with_remodeller(X)
                state_Y, helicase_Y  = state_mapping_with_remodeller(Y)
                if ∑(abs.(state_X-state_Y))+∑(abs.(helicase_X-helicase_Y)) == 1
                    if ∑(state_X-state_Y) < 0
                        # return the index i at which state_X[i] < state_Y[i]
                        return findfirst(state_X .< state_Y)
                    else
                        return 0
                    end
                else
                    return 0
                end
            else
                return 0
            end
        else 
            return 0
        end
    end
    
    function helicase_on(x::Int,y::Int)
        X = S0[x]
        Y = S0[y]
        if length(X) == length(Y)
            Δ = X - Y
            if Δ[1:3] == [0,0,0]
            state_X, helicase_X = state_mapping_with_remodeller(X)
            state_Y, helicase_Y  = state_mapping_with_remodeller(Y)
                if ∑(abs.(state_X-state_Y))+∑(abs.(helicase_X-helicase_Y)) == 1
                    if ∑(helicase_X-helicase_Y) > 0
                        return 1
                    else
                        return 0
                    end
                else
                    return 0
                end
            else
                return 0
            end
        else 
            return 0
        end
    end
    
    function helicase_off(x::Int,y::Int)
        X = S0[x]
        Y = S0[y]
        if length(X) == length(Y)
            Δ = X - Y
            if Δ[1:3] == [0,0,0]
            state_X, helicase_X = state_mapping_with_remodeller(X)
            state_Y, helicase_Y  = state_mapping_with_remodeller(Y)
                if ∑(abs.(state_X-state_Y))+∑(abs.(helicase_X-helicase_Y)) == 1
                    if ∑(helicase_X-helicase_Y) < 0
                        return 1
                    else
                        return 0
                    end
                else
                    return 0
                end
            else
                return 0
            end
        else 
            return 0
        end
    end
    
    function dissociation(x::Int,y::Int)
        X = S0[x]
        Y = S0[y]
        state_X, helicase_X = state_mapping_with_remodeller(X)
        state_Y, helicase_Y= state_mapping_with_remodeller(Y) # restriction on microstates is imposed by the fact that different macrostates have different restrictions on microstates.
        if ∑(abs.(state_X-state_Y)) == 1
            inds = [i for i in 1:3 if X[i] != Y[i]]
            if length(inds) == 1 && 3 ∉ inds # lost of one of (H2A-H2B)
                ind = inds[1]
                if X[ind] == 0 && Y[ind] == 1
                    return findfirst(state_X .< state_Y)
                else
                    return 0
                end
            elseif length(inds) == 1 && 3 ∈ inds # lost of (H3-H4)₂, note that it CAN occur when X = [0,0,0], Y = [0,1,0] and CANNOT occur when either one of H2A-H2B attaches to (H3-H4)₂.
                ind = inds[1]
                if maximum(X[1:3]) < 2 && maximum(Y[1:3]) < 2 && X[3] == 0 && Y[3] == 1 # no connection is established between (H3-H4)₂ and (H2A-H2B)
                    return findfirst(state_X .< state_Y)
                else
                    return 0
                end
            elseif length(inds) == 2 && sort(Y[inds]) == [1,2] && sort(X[inds]) == [0,0] && 3 ∈ inds 
                return findfirst(state_X .< state_Y)
            elseif length(inds) == 3 && Y[inds] == [2,2,1] && X[inds] == [0,0,0]
                return findfirst(state_X .< state_Y)
            else
                return 0
            end
        else
            return 0
        end
    end

    ## initialize state space
    S0 = Array{Any,1}[]
    # From MacroStates to Possible Pure microstates, macrostates pair 
    # From these pair to construct extended form microstate (0/1s) and possible helicase position
    for MS in MacroStates
        X,Y,Z = MS
        numbers = macro_state2sizes(X,Y,Z)
        if length(numbers) == 0
            state = [X,Y,Z]
            microstate = state_mapping(state)
            for υ in get_plausible_helicase_postion(microstate)
                push!(S0,[X,Y,Z,υ])
            end
        elseif length(numbers) == 1
            a = numbers[1]
            for ϕ in Φ(a)
                state = [X,Y,Z,ϕ]
                microstate = state_mapping(state)
                for υ in get_plausible_helicase_postion(microstate)
                    push!(S0,[X,Y,Z,ϕ,υ])
                end
            end
        elseif length(numbers) == 2
            a,b = numbers
            for ϕ in Φ(a), ψ in Φ(b)
                state = [X,Y,Z,ϕ,ψ]
                microstate = state_mapping(state)
                for υ in get_plausible_helicase_postion(microstate)
                    push!(S0,[X,Y,Z,ϕ,ψ,υ])
                end
            end
        elseif length(numbers) == 3
            a,b,c = numbers
            for ϕ in Φ(a), ψ in Φ(b), χ in Φ(c)
                state = [X,Y,Z,ϕ,ψ,χ]
                microstate = state_mapping(state)
                for υ in get_plausible_helicase_postion(microstate)
                    push!(S0,[X,Y,Z,ϕ,ψ,χ,υ])
                end
            end
        end
    end

    push!(MacroStates,[0,0,0])# add the empty state
    push!(S0,[0,0,0]) # add the empty state

    N_total = length(S0)


    # if isfile("simple_facilitated_matrix_data.jld")
    #     d = load("simple_facilitated_matrix_data.jld")
    #     M_helicase_on = d["helicase_on"]
    #     M_helicase_off = d["helicase_off"]
    #     M_internal_on = d["on"]
    #     M_internal_off = d["off"]
    #     M_dissociation = d["dissociation"]
    # else
        M_helicase_on = zeros(Int32,N_total, N_total)
        Threads.@threads for x in 1:N_total
            for y in 1:N_total
                if is_adjacent(x,y)
                    if helicase_on(x,y) > 0
                        M_helicase_on[x,y] = helicase_on(x, y)
                    end
                end
            end
        end
        M_helicase_off = zeros(Int32,N_total, N_total)
        Threads.@threads for x in 1:N_total
            for y in 1:N_total
                # println("x=$x, y=$y")
                if is_adjacent(x,y)
                    if helicase_off(x,y) > 0
                        M_helicase_off[x,y] = helicase_off(x, y)
                    end
                end
            end
        end
        M_internal_on = zeros(Int32,N_total, N_total)
        Threads.@threads for x in 1:N_total
            for y in 1:N_total
                if is_adjacent(x,y)
                    if internal_on(x,y) > 0
                        M_internal_on[x,y] = internal_on(x, y)
                    end
                end
            end
        end
        M_internal_off = zeros(Int32,N_total, N_total)
        Threads.@threads for x in 1:N_total
            for y in 1:N_total
                if is_adjacent(x,y)
                    if internal_off(x, y) > 0
                        M_internal_off[x,y] = internal_off(x, y)
                    end
                end
            end
        end
        M_dissociation = zeros(Int32,N_total, N_total)
        Threads.@threads for x in 1:N_total
            for y in 1:N_total
                if is_adjacent(x,y)
                    if dissociation(x,y) > 0
                        M_dissociation[x,y] = dissociation(x, y)
                    end
                end
            end
        end
    #     save("simple_facilitated_matrix_data.jld","helicase_on",M_helicase_on,"helicase_off",M_helicase_off,"on",M_internal_on,"off",M_internal_off,"dissociation",M_dissociation)
    # end


    function TransitionKernel(k_on,E_c,q_on,q_off)
        # derived parameters
        # we assume that k_on is identical for different E_b_landscape
        # only k_off is different
        k_off = exp.(E_c) .* k_on
        k = [k_on,q_on,q_off]
        M = [M_internal_on,M_helicase_on,M_helicase_off]
        Q = sum([k[i]*M[i] for i in eachindex(k)])
        # the nonzero value of M_internal_on and M_dissociation
        # represents the index of the corresponding histone-DNA bond
        for ind in eachindex(Q)
            if M_internal_off[ind] > 0
                Q[ind] += k_off[M_internal_off[ind]]
            elseif M_dissociation[ind] > 0
                Q[ind] += k_off[M_dissociation[ind]]
            end
        end
        # rectification
        for x in 1:N_total
            Q[x,x] = - sum([Q[y,x] for y in 1:N_total if y != x])
        end
        # the last state is the fully dissociated state, which we don't need
        return Q[1:(N_total-1),1:(N_total-1)]
    end

    # λ₀ calculate the principal eigenvalue of the transition kernel
    function λ₀(k_on,E_c,q_on,q_off)
        Q = TransitionKernel(k_on,E_c,q_on,q_off)
        eigs = eigvals(Q)
        return abs(eigs[end])
    end

    # λ̂₀ provides an elementary estiamte of the principal eigenvalue of the transition kernel
    function λ̂₀(k_on,E_c,q_on,q_off)
        
        k_off = exp(mean(E_c)) * k_on
        # steady state estimates
        if q_off == 0
            ss = k_off
        else
            ss = exp(N*mean(E_c))*∑([
                (N-l)*(l+1)*(q_on/q_off)^l for l ∈ 0:1:N-1
                ]) / ∑([
                    (l+1)*(exp(mean(E_c))*(q_on/q_off))^l for l in 0:1:N-1]
                )
        end
        # perturbation estimates
        perturb = minimum([
            k_off,
            (N * k_on * exp(N*mean(E_c)) + q_on * ∑(
                [(k+1) * exp(k*mean(E_c)) for k ∈ 1:1:N-1]
            )) / ∑(
                [(k+1) * exp(k*mean(E_c)) for k ∈ 0:1:N]
            )
        ])
        return minimum([exp(mean(E_c))*k_on,ss,perturb])
        # return minimum([exp(E_c)*q_on,exp(E_c)*k_on,maximum([N*(q_on/q_off)^(N-1)*exp(N*E_c),N*exp(N*E_c)])])
    end

    function f(k_on,E_c,q_on,q_off)
        k_off = exp(mean(E_c))*k_on
        b = k_off + q_off - exp(mean(E_c))*q_on
        frac = 2*exp(mean(E_c))*q_on/(b+√(b^2+4*exp(mean(E_c))*q_on*q_off))
        interm = exp(N*mean(E_c))*(∑([(N-l)*(l+1)*(q_on/q_off)^l for l in 0:1:N-2])+N * (q_on/q_off)^(N-2)* frac/exp(mean(E_c)))/(∑([(l+1)*(exp(mean(E_c))*(q_on/q_off))^l for l in 0:1:N-2])+ N * exp((N-2)*mean(E_c))* (q_on/q_off)^(N-2)*frac )
        return interm
    end
    f′(k_on,E_c,q_on,q_off) = exp(N*mean(E_c))*∑([(N-l)*(l+1)*(q_on/q_off)^l for l in 0:1:N-1])/∑([(l+1)*(exp(mean(E_c))*(q_on/q_off))^l for l in 0:1:N-1])

    # q_off = 0.01
    # plot(legend=false,xaxis=:log, yaxis=:log)
    # plot!(q_onss,[f(k_on,E_c,q_on,q_off) for q_on in q_onss])
    # plot!(q_onss,[f′(k_on,E_c,q_on,q_off) for q_on in q_onss])

    function higher_limit_estimate(k_on,E_c,q_on,q_off)
        k_off = exp(mean(E_c))*k_on
        δ = q_off/(q_on*exp(mean(E_c)))
        return minimum([q_on*(exp(mean(E_c)))*(1- ((δ^N)/(1-δ))^(1/(N-2))),k_off])
    end

    function dλ(k_on,E_c,q_on,q_off)
        k_off = exp(mean(E_c))*k_on
        q̃_on = q_on*k_off/k_on
        frac = q̃_on/q_off
        frac = 1/(1/frac - log(q_on)/N)
        J = ∑([(N-k)*(k+1) * k_on * (k_off/k_on)^(N-k) * (frac)^k for k in 0:1:N-1])
        Ω = ∑([(k+1)*(frac)^k for k in 0:1:N-1])
        return J/Ω
    end

    df = DataFrame(
        k_on=Float64[], 
        E_c=Float64[], 
        q_on=Float64[], 
        q_off=Float64[], 
        N=Int[], 
        λ₀=Float64[],
        id = Int[],
        )
    id = 0
    @showprogress "calculating eigenvalues" for E_c in E_bs
        id += 1
            for q_on in q_ons, q_off in q_offs
                @show q_on, q_off
               push!(df, [k_on, mean(E_c), q_on, q_off, N, λ₀(k_on,E_c,q_on,q_off), id])
            end
        end
    
    CSV.write("linear_facilitated_results_$(N)_random_landscape.csv", df)

    df = DataFrame(
        k_on=Float64[], 
        E_c=Float64[], 
        q_on=Float64[], 
        q_off=Float64[], 
        N=Int[], 
        λ̂₀=Float64[],
        id = Int[],
        )
    id = 0
    # save results to DataFrame
    @showprogress "estimating eigenvalues" for E_c in E_bs 
        id += 1
            for q_on in q_ons, q_off in q_offs, E_c in E_bs
                push!(df, [k_on, mean(E_c), q_on, q_off, N, λ̂₀(k_on,E_c,q_on,q_off), id])
            end
    end

    CSV.write("linear_facilitated_results_$(N)_random_landscape_approximate.csv", df)
end

k_on = 1
E_bs = [E_b_landscape(14) for i in 1:5]
q_offs = [10.0^(-3)]
q_ons = [10.0^(i) for i in -14:0.5:1]
id = 0
for N ∈ 14:14
    println("N = $N")
    data_acquire(N, k_on, E_bs, q_ons, q_offs)
end