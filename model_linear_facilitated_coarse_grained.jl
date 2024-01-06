using SparseArrays
using LinearAlgebra
using ProgressBars
using Parameters
using DataFrames
using CSV
∑ = sum
@with_kw struct Par
    N::Int64 = 14
    k_on::Float64 = 1.
    E_c::Float64 = -2.
    q_on::Float64 = 0.1
    q_off::Float64 = 0.01
end
BLAS.set_num_threads(20)

# let x, y ∈ states. the rates between them is given by:
function transition_rate(
    pair::Pair{Tuple{Int64, Int64}, Tuple{Int64, Int64}};
    par=par)
    x, y = pair
    if abs.(x.-y) |> ∑ > 1
        return 0
    else
        Δ = y .- x
        if Δ == (1,0) || Δ == (0,1)
            return par.q_on * exp(par.E_c)
        elseif Δ == (-1,0) || Δ == (0,-1)
            return par.q_off
        else
            return 0
        end
    end
end

function flux_out(state::Tuple; par=par)
    N = par.N
    k_on = par.k_on
    E_c = par.E_c
    n = N - ∑(state)
    return n * k_on * exp(n*E_c)
end


function principal_eigenvalue(par::Par)
    N = par.N
    # numerate the coarse grained states 
    # the following is a vector of all possible states
    states = [(i,j) for i ∈ 0:1:N-1, j ∈ 0:1:N-1 if i+j ≤ N-1]
    # build the transition matrix
    transition_matrix = zeros(length(states), length(states))
    for (i, x) in enumerate(states)
        for (j, y) in enumerate(states)
            transition_matrix[i,j] = transition_rate(x=>y; par=par)
        end
    end
    for (i, x) in enumerate(states)
        transition_matrix[i,i] = -(∑(transition_matrix[i,:]) 
                                    + flux_out(x; par=par))
    end
    # calculate the principal eigenvalue
    λ₀ = eigvals(transition_matrix)[end]
    return λ₀ |> abs
end


N = 14
k_on = 1
E_bs = [-2]
q_offs = [0, 10.0^(-5), 10.0^(-3),  10.0^(-1), 10.0^1]
q_ons = [10.0^(i) for i in -14:0.2:3]


pars = [Par(N=N, k_on=k_on, E_c=E_c, q_on=q_on, q_off=q_off) 
        for E_c in E_bs, q_on in q_ons, q_off in q_offs if true]

# "if true" argument turns pars into 1d array

λ₀s = zeros(length(pars))

for it in ProgressBar(1:length(pars))
    par = pars[it]
    λ₀s[it] = principal_eigenvalue(par)
end

# collect the results into a DataFrame
result = DataFrame(
    k_on= [par.k_on for par in pars],
    E_c= [par.E_c for par in pars],
    q_on= [par.q_on for par in pars],
    q_off= [par.q_off for par in pars],
    N= [par.N for par in pars],
    λ₀= λ₀s,
    )

CSV.write("linear_facilitated_results_coarse_grained.csv", result)
