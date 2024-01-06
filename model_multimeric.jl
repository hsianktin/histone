
∑(A) = sum(A)
using LinearAlgebra
using Plots
using JLD
using ArgParse
### Parameter settings and methods
using Parameters
using DataFrames, CSV
@with_kw struct Par
    k_on::Float64 = 1
    E_c::Float64 = - 2
    E_L::Float64 = - 0.5
    E_P_H2d::Float64 = Inf
    E_P_H3_H4::Float64 = Inf
    k_undocking_simple::Float64 = 0.01
    k_undocking_coupled::Float64 = 0.01
end

BLAS.set_num_threads(24)
## specifying macrostates 

s = ArgParseSettings()

@add_arg_table! s begin
    "profile"
    help = "input profile file"
    default = "profiles/multimeric_reversible_EL.jl"
    required = false

    "--calculate"
    help = "whether to numerically calculate the eigenvalue"
    default = false
    arg_type = Bool
    required = false
end

profile = parse_args(ARGS, s)["profile"]
calculate_flag = parse_args(ARGS, s)["calculate"]
include(profile)

N = 14
NL = 4
NR = 4
NM = N - NL - NR

MacroStates = vcat(vec([(i,j,1) for i in 0:2, j in 0:2]),vec([(i,j,0) for i in 0:1, j in 0:1 if i + j >0]))
push!(MacroStates,(0,0,0))

include("utils.jl")

## initialize state space
S0 = Array{Any,1}[]

for MS in MacroStates
    X,Y,Z = MS
    numbers = macro_state2sizes(X,Y,Z)
    if length(numbers) == 0
        push!(S0,[X,Y,Z])
    elseif length(numbers) == 1
        a = numbers[1]
        for ϕ in Φ(a)
            push!(S0,[X,Y,Z,ϕ])
        end
    elseif length(numbers) == 2
        a,b = numbers
        for ϕ in Φ(a), ψ in Φ(b)
            push!(S0,[X,Y,Z,ϕ,ψ])
        end
    elseif length(numbers) == 3
        a,b,c = numbers
        for ϕ in Φ(a), ψ in Φ(b), χ in Φ(c)
            push!(S0,[X,Y,Z,ϕ,ψ,χ])
        end
    end
end

N_total = length(S0)
# find the index in S0 corresponding to  [2,2,1,(0,0)]
s₀ = findall(x->state_mapping(x)==ones(N),S0)
s₁ = findall(x->∑(state_mapping(x))==N-1,S0)
for i in s₁
    @show S0[i]
end

# # test validity of the mapping algorithm, which is the key of the model.
# state = S0[rand(1:length(S0))]
# bar(state_mapping(state),title=state,size = (800,200))

## Defining energetics (existence of detailed equilibrium)
# E_c = log(k_off/k_on)
# E_L = log(k_undocking/k_docking)
# E_L′ = log(k_undocking/k_docking′)
# E_L′′ = log(k_undocking/k_docking′′)
# E_P = log(k_on/k_association)
# E_P′ = log(k_on/k_association′)
# E_P′′ = log(k_on/k_association′′)

# EQ1 = (E_L′-E_L-E_P)
# EQ2 = (E_L′′-2E_L-2E_P)
# EQ3 = (E_P′ - E_L - 2E_P)
# EQ4 = (E_P′′ - 2E_L - 3E_P)
# EQS = [EQ1,EQ2,EQ3,EQ4]
## Defining dynamics


## contruct rate matrices
if calculate_flag
if isfile("matrix_data_spontaneous_multimeric.jld")
    d = load("matrix_data_spontaneous_multimeric.jld")
    M_internal_on = d["on"]
    M_internal_off = d["off"]
    M_docking_simple = d["docking_simple"]
    M_docking_H2d = d["docking_H2d"]
    M_docking_H3_H4 = d["docking_H3_H4"]
    M_docking_HEX = d["docking_HEX"]
    M_undocking_simple = d["undocking_simple"]
    M_undocking_coupled = d["undocking_coupled"]
    M_association_H2d = d["association_H2d"]
    M_association_H3_H4 = d["association_H3_H4"]
    M_association_HEX = d["association_HEX"]
    M_association_OCT = d["association_OCT"]
    M_dissociation = d["dissociation"]
else
    M_internal_on = zeros(N_total, N_total)
    Threads.@threads for x in 1:N_total
        for y in 1:N_total
            if is_adjacent(x,y)
                M_internal_on[x,y] = internal_on(x, y)
            end
        end
    end
    M_internal_off = zeros(N_total, N_total)
    Threads.@threads for x in 1:N_total
        for y in 1:N_total
            if is_adjacent(x,y)
                M_internal_off[x,y] = internal_off(x, y)
            end
        end
    end
    M_docking_simple = zeros(N_total, N_total)
    Threads.@threads for x in 1:N_total
        for y in 1:N_total
            if is_adjacent(x,y)
                M_docking_simple[x,y] = docking(x, y,"simple")
            end
        end
    end
    M_docking_H2d = zeros(N_total, N_total)
    Threads.@threads for x in 1:N_total
        for y in 1:N_total
            if is_adjacent(x,y)
                M_docking_H2d[x,y] = docking(x, y,"H2d")
            end
        end
    end
    M_docking_H3_H4 = zeros(N_total, N_total)
    Threads.@threads for x in 1:N_total
        for y in 1:N_total
            if is_adjacent(x,y)
                M_docking_H3_H4[x,y] = docking(x, y,"H3_H4")
            end
        end
    end
    M_docking_HEX = zeros(N_total, N_total)
    Threads.@threads for x in 1:N_total
        for y in 1:N_total
            if is_adjacent(x,y)
                M_docking_HEX[x,y] = docking(x, y,"HEX")
            end
        end
    end
    M_undocking_simple = zeros(N_total, N_total)
    Threads.@threads for x in 1:N_total
        for y in 1:N_total
            if is_adjacent(x,y)
                M_undocking_simple[x,y] = undocking(x, y,"simple")
            end
        end
    end
    M_undocking_coupled = zeros(N_total, N_total)
    Threads.@threads for x in 1:N_total
        for y in 1:N_total
            if is_adjacent(x,y)
                M_undocking_coupled[x,y] = undocking(x, y, "coupled")
            end
        end
    end
    M_dissociation = zeros(N_total, N_total)
    Threads.@threads for x in 1:N_total
        for y in 1:N_total
            if is_adjacent(x,y)
                M_dissociation[x,y] = dissociation(x, y)
            end
        end
    end
    M_association_H2d = zeros(N_total, N_total)
    Threads.@threads for x in 1:N_total
        for y in 1:N_total
            if is_adjacent(x,y)
                M_association_H2d[x,y] = association(x, y, "H2d")
            end
        end
    end
    M_association_H3_H4 = zeros(N_total, N_total)
    Threads.@threads for x in 1:N_total
        for y in 1:N_total
            if is_adjacent(x,y)
                M_association_H3_H4[x,y] = association(x, y,"H3_H4")
            end
        end
    end
    M_association_HEX = zeros(N_total, N_total)
    Threads.@threads for x in 1:N_total
        for y in 1:N_total
            if is_adjacent(x,y)
                M_association_HEX[x,y] = association(x, y,"HEX")
            end
        end
    end
    M_association_OCT = zeros(N_total, N_total)
    Threads.@threads for x in 1:N_total
        for y in 1:N_total
            if is_adjacent(x,y)
                M_association_OCT[x,y] = association(x, y,"OCT")
            end
        end
    end
    save("matrix_data_spontaneous_multimeric.jld","on",M_internal_on,"off",M_internal_off,"docking_simple",M_docking_simple,"docking_H2d",M_docking_H2d,"docking_H3_H4",M_docking_H3_H4,"docking_HEX",M_docking_HEX,"undocking_simple",M_undocking_simple,"undocking_coupled",M_undocking_coupled,"association_H2d",M_association_H2d,"association_H3_H4",M_association_H3_H4,"association_HEX",M_association_HEX,"association_OCT",M_association_OCT,"dissociation",M_dissociation)
end
end
# controlling parameters

function TransitionKernel(k_on,E_c,E_L,E_P_H2d,E_P_H3_H4,k_undocking_simple,k_undocking_coupled)
    # derived parameters
    k_off = exp(E_c) * k_on
    k_docking_simple = exp(-E_L) * k_undocking_simple
    k_docking_H2d = exp(-E_L-E_P_H2d) * k_undocking_coupled
    k_docking_H3_H4 = exp(-E_L-E_P_H3_H4) * k_undocking_coupled
    k_docking_HEX = exp(-2E_L - E_P_H2d - E_P_H3_H4) * k_undocking_coupled
    k_association_H2d = exp(-E_P_H2d) * k_on
    k_association_H3_H4 = exp(-E_P_H3_H4) * k_on
    k_association_HEX = exp(-E_L - E_P_H2d - E_P_H3_H4) * k_on
    k_association_OCT = exp(-2E_L - 2E_P_H2d - E_P_H3_H4) * k_on
    k = [k_on,k_off,k_docking_simple,k_docking_H2d,k_docking_H3_H4,k_docking_HEX,k_undocking_simple,k_undocking_coupled,k_association_H2d,k_association_H3_H4,k_association_HEX,k_association_OCT,k_off]
    M = [M_internal_on,M_internal_off,M_docking_simple,M_docking_H2d,M_docking_H3_H4,M_docking_HEX,M_undocking_simple,M_undocking_coupled,M_association_H2d,M_association_H3_H4,M_association_HEX, M_association_OCT,M_dissociation]
    Q = sum([k[i]*M[i] for i in eachindex(k)])
    # rectification
    for x in 1:N_total
        Q[x,x] = - sum([Q[y,x] for y in 1:N_total if y != x])
    end
    return Q
end
function E(state,E_c,E_L,E_P_H2d,E_P_H3_H4)
    i,j,k = state[1:3]
    function interpret(i)
        if i == 0
            return [0,0]
        elseif i == 1
            return [1,0]
        else
            return [1,1]
        end
    end
    P_l, L_l = interpret(i)
    P_m = k
    P_r, L_r = interpret(j)
    microstate = state_mapping(state)
    Energy = E_P_H2d*(P_l + P_r) + E_P_H3_H4*P_m + E_L*(L_l + L_r) + E_c*∑(microstate)
    return Energy
end

# estimate of λ₁ based on flow
function flow_parameter(k_on,E_c,E_L,E_P_H2d,E_P_H3_H4,k_undocking_simple,k_undocking_coupled)
    # println("------------------------------------------------")
    k_off = exp(E_c) * k_on
    k_docking_simple = exp(-E_L) * k_undocking_simple
    k_docking_H2d = exp(-E_L-E_P_H2d) * k_undocking_coupled
    k_docking_H3_H4 = exp(-E_L-E_P_H3_H4) * k_undocking_coupled
    k_docking_HEX = exp(-2E_L - E_P_H2d - E_P_H3_H4) * k_undocking_coupled
    k_association_H2d = exp(-E_P_H2d) * k_on
    k_association_H3_H4 = exp(-E_P_H3_H4) * k_on
    k_association_HEX = exp(-E_L - E_P_H2d - E_P_H3_H4) * k_on
    k_association_OCT = exp(-2E_L - 2E_P_H2d - E_P_H3_H4) * k_on
    # println(k_docking_H2d)
    (k_docking_simple +  exp(NL*E_c + E_P_H2d) * k_docking_H2d)
    2k_association_H2d + k_association_OCT + 2k_association_HEX
    k_docking_H2d+k_association_H2d
    # println(minimum([k_docking_H2d+k_association_H2d,exp(-(E_P_H2d+E_L))*k_on]))
    # println(minimum([k_association_H2d + k_association_OCT + k_association_HEX,exp(-(E_P_H2d+E_c))*k_on]))

    # println(minimum([k_docking_H2d,exp(-(E_P_H2d+E_L))*k_on]),minimum([k_association_H2d,exp(-(E_P_H2d+E_c))*k_on]))
    if E_P_H2d + NL * E_c <= 0
        return maximum(
            [
                minimum(
                    [
                        k_docking_simple +  exp(NL*E_c + E_P_H2d) * k_docking_H2d, 
                        (k_docking_H2d+k_association_H2d),
                        # exp(-(E_P_H2d+E_L))*k_on, 
                        2k_association_H2d + k_association_OCT + 2k_association_HEX,
                        # exp(-(E_P_H2d+E_c))*k_on,
                        (k_off + k_undocking_simple)/(exp(E_c) + exp(E_L)),
                        k_on
                    ]),
                N*exp((N-1)*E_c)*k_off,
            ]
        )
    elseif 2E_P_H2d + E_P_H3_H4 + N*E_c + 2E_L > 0
        return minimum([NM * exp(NM*E_c) *k_on])
    else
        return maximum(
            [
                minimum(
                    [
                        k_docking_simple, 
                        k_docking_H2d+k_association_H2d,
                        # exp(-(E_P_H2d+E_L))*k_on,
                        k_association_H2d + k_association_OCT + k_association_HEX,
                        # exp(-(E_P_H2d+E_c))*k_on,
                        k_on,
                        ]
                    ), 
                N*exp((N-1)*E_c)*k_off
            ]
            )
    end
end

using Arpack
function λs(k_on,E_c,E_L,E_P_H2d,E_P_H3_H4,k_undocking_simple,k_undocking_coupled)
    # @info "Computing the largest eigenvalue and eigenvector..."
    λs = eigvals(TransitionKernel(k_on,E_c,E_L,E_P_H2d,E_P_H3_H4,k_undocking_simple,k_undocking_coupled))
    # λs, vs = eigs(
    #     TransitionKernel(
    #         k_on,
    #         E_c,
    #         E_L,
    #         E_P_H2d,
    #         E_P_H3_H4,
    #         k_undocking_simple,
    #         k_undocking_coupled
    #         );
    #         nev = 2, # number of eigenvalues to compute
    #         which = :LR, # largest real part
    #         maxiter = 10000000,
    #         tol = 1e-12,
    #         )
    # sort eigenvalues in descending order by real part
    # λs = sort(λs, by = x -> real(x), rev = true)
    # @info λs[1:2]
    return -1 .* real.(
        λs[end:-1:end-1]
    )
end



function λs(p::Par)
    return λs(
            p.k_on,
            p.E_c,
            p.E_L,
            p.E_P_H2d,
            p.E_P_H3_H4,
            p.k_undocking_simple,
            p.k_undocking_coupled
        )
end

result_df = DataFrame(
    k_on = Float64[],
    E_c = Float64[],
    E_L = Float64[],
    E_P_H2d = Float64[],
    E_P_H3_H4 = Float64[],
    k_undocking_simple = Float64[],
    k_undocking_coupled = Float64[],
    λ₀ = Float64[],
    λ₁ = Float64[]
)



using ProgressBars
if calculate_flag
    prog = ProgressBar(1:length(pars))
    for id in prog
        p = pars[id]
        λ = λs(p)
        push!(result_df, (p.k_on, p.E_c, p.E_L, p.E_P_H2d, p.E_P_H3_H4, p.k_undocking_simple, p.k_undocking_coupled, λ[1] ,λ[2]))
        set_multiline_postfix(prog, "k_on = $(p.k_on)\n E_c = $(p.E_c)\n E_L = $(p.E_L)\n k_undocking_simple = $(p.k_undocking_simple)\n k_undocking_coupled = $(p.k_undocking_coupled)\n λ₀ = $(λ[1])\n λ₁ = $(λ[2]) ")
    end

    CSV.write("eigvals_model_multimeric_$(label).csv", result_df)
end
function λ̂₁(p::Par)
    # estimate based on mean first passage time
    k_on = p.k_on
    E_c = p.E_c
    k_undocking_simple = p.k_undocking_simple
    k_undocking_coupled = p.k_undocking_coupled
    k_docking = k_undocking_simple * exp(-p.E_L)
    p_docked = k_docking / (k_undocking_simple + k_docking)
    p_reattach = k_on / (k_on + k_undocking_coupled)
    τ_L = 1/(k_on * exp(NL * E_c))
    τ_M = 1/(NM* k_on * exp(NM * E_c))
    return maximum([minimum([(τ_L / (1 - p_docked * p_reattach) + τ_M)^(-1), k_undocking_coupled]), N * exp((N)*E_c)*k_on])
end

estimate_df = DataFrame(
    k_on = Float64[],
    E_c = Float64[],
    E_L = Float64[],
    E_P_H2d = Float64[],
    E_P_H3_H4 = Float64[],
    k_undocking_simple = Float64[],
    k_undocking_coupled = Float64[],
    λ̂₁ = Float64[]
)

Z(p) = ∑([exp(-E(S0[x],p.E_c,p.E_L,p.E_P_H2d,p.E_P_H3_H4)) for x in eachindex(S0)])
π₀(p) = [exp(-E(S0[x],p.E_c,p.E_L,p.E_P_H2d,p.E_P_H3_H4))/Z(p) for x in eachindex(S0)]
for p in (pars)
    begin
        # Z = ∑([exp(-E(S0[x],E_c,E_L,E_P_H2d,E_P_H3_H4)) for x in eachindex(S0)])
        # p₀=[exp(-E(S0[x],E_c,E_L,E_P_H2d,E_P_H3_H4))/Z for x in eachindex(S0)]

        λ̂ = λ̂₁(p)
        push!(estimate_df, (p.k_on, p.E_c, p.E_L, p.E_P_H2d, p.E_P_H3_H4, p.k_undocking_simple, p.k_undocking_coupled, λ̂))
    end
end

CSV.write("eigvals_model_multimeric_$(label)_estimate.csv", estimate_df)

# wrapper function for flow parameter
function flow_parameter(p::Par)
    return flow_parameter(p.k_on, p.E_c, p.E_L, p.E_P_H2d, p.E_P_H3_H4, p.k_undocking_simple, p.k_undocking_coupled)
end

estimate_df = DataFrame(
    k_on = Float64[],
    E_c = Float64[],
    E_L = Float64[],
    E_P_H2d = Float64[],
    E_P_H3_H4 = Float64[],
    k_undocking_simple = Float64[],
    k_undocking_coupled = Float64[],
    λ̂₁ = Float64[]
)

for p in (pars)
    begin
        # Zₚ = Z(p)
        # πₚ = π₀(p)
        # @show p
        # @show ∑([πₚ[i] for i  in s₁]) / (1 - ∑([πₚ[i] for i  in s₀]))
        # @show πₚ[end]
        λ̂ = flow_parameter(p)
        push!(estimate_df, (p.k_on, p.E_c, p.E_L, p.E_P_H2d, p.E_P_H3_H4, p.k_undocking_simple, p.k_undocking_coupled, λ̂))
    end
end

CSV.write("eigvals_model_multimeric_$(label)_flow_parameter.csv", estimate_df)