
∑(A) = sum(A)
using SparseArrays
using LinearAlgebra
using Arpack
using Plots
using Hwloc
using JLD
using Parameters
using DataFrames, CSV, ProgressBars
using BSON: @save, @load
using ArgParse

BLAS.set_num_threads(num_physical_cores())
@show BLAS.get_num_threads()
@with_kw struct Par
    k_on::Float64 = 1
    E_c::Float64 = - 2
    E_L::Float64 = - 0.5
    E_P_H2d::Float64 = Inf
    E_P_H3_H4::Float64 = Inf
    k_undocking_simple::Float64 = 0.01
    k_undocking_coupled::Float64 = 0.01
    q_on::Float64 = 0.01
    q_off::Float64 = 0.0001
end
s = ArgParseSettings()

@add_arg_table! s begin
    "profile"
    help = "input profile file"
    default = "profiles/multimeric_facilitated_irreversible.jl"
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
## specifying microstates
include("utils.jl")

N = 14
NL = 4
NR = 4
NM = N - NL - NR

MacroStates = vcat(vec([(i,j,1) for i in 0:2, j in 0:2]),vec([(i,j,0) for i in 0:1, j in 0:1]))


## initialize state space
S0 = Array{Any,1}[]




###################### state space construction ######################
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
########################################################################


N_total = length(S0)

if calculate_flag
if isfile("matrix_data_facilitated_multimeric.jld")
    @info "Loading transition matrices from file..."
    d = load("matrix_data_facilitated_multimeric.jld")
    M_helicase_on = d["helicase_on"]
    M_helicase_off = d["helicase_off"]
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
    M_dissociation = d["dissociation_with_remodeller"]
else
    @info "Building the transition matrices..."
    M_helicase_on = zeros(Int32,N_total, N_total)
    @time Threads.@threads for x in 1:N_total
        for y in 1:N_total
            if is_adjacent_with_remodeller(x,y)
                M_helicase_on[x,y] = helicase_on(x, y)
            end
        end
    end
    M_helicase_on = sparse(M_helicase_on)
    M_helicase_off = zeros(Int32,N_total, N_total)
    @time Threads.@threads for x in 1:N_total
        for y in 1:N_total
            if is_adjacent_with_remodeller(x,y)
                M_helicase_off[x,y] = helicase_off(x, y)
            end
        end
    end
    M_helicase_off = sparse(M_helicase_off)
    M_internal_on = zeros(Int32,N_total, N_total)
    @time Threads.@threads for x in 1:N_total
        for y in 1:N_total
            if is_adjacent_with_remodeller(x,y)
                M_internal_on[x,y] = internal_on_with_remodeller(x, y)
            end
        end
    end
    M_internal_on = sparse(M_internal_on)
    M_internal_off = zeros(Int32,N_total, N_total)
    @time Threads.@threads for x in 1:N_total
        for y in 1:N_total
            if is_adjacent_with_remodeller(x,y)
                M_internal_off[x,y] = internal_off_with_remodeller(x, y)
            end
        end
    end
    M_internal_off = sparse(M_internal_off)
    M_docking_simple = zeros(Int32,N_total, N_total)
    @time Threads.@threads for x in 1:N_total
        for y in 1:N_total
            if is_adjacent_with_remodeller(x,y)
                M_docking_simple[x,y] = docking_with_remodeller(x, y,"simple")
            end
        end
    end
    M_docking_simple = sparse(M_docking_simple)
    M_docking_H2d = zeros(Int32,N_total, N_total)
    @time Threads.@threads for x in 1:N_total
        for y in 1:N_total
            if is_adjacent_with_remodeller(x,y)
                M_docking_H2d[x,y] = docking_with_remodeller(x, y,"H2d")
            end
        end
    end
    M_docking_H2d = sparse(M_docking_H2d)
    M_docking_H3_H4 = zeros(Int32,N_total, N_total)
    @time Threads.@threads for x in 1:N_total
        for y in 1:N_total
            if is_adjacent_with_remodeller(x,y)
                M_docking_H3_H4[x,y] = docking_with_remodeller(x, y,"H3_H4")
            end
        end
    end
    M_docking_H3_H4 = sparse(M_docking_H3_H4)
    M_docking_HEX = zeros(Int32,N_total, N_total)
    @time Threads.@threads for x in 1:N_total
        for y in 1:N_total
            if is_adjacent_with_remodeller(x,y)
                M_docking_HEX[x,y] = docking_with_remodeller(x, y,"HEX")
            end
        end
    end
    M_docking_HEX = sparse(M_docking_HEX)
    M_undocking_simple = zeros(Int32,N_total, N_total)
    @time Threads.@threads for x in 1:N_total
        for y in 1:N_total
            if is_adjacent_with_remodeller(x,y)
                M_undocking_simple[x,y] = undocking_with_remodeller(x, y,"simple")
            end
        end
    end
    M_undocking_simple = sparse(M_undocking_simple)
    M_undocking_coupled = zeros(Int32,N_total, N_total)
    @time Threads.@threads for x in 1:N_total
        for y in 1:N_total
            if is_adjacent_with_remodeller(x,y)
                M_undocking_coupled[x,y] = undocking_with_remodeller(x, y, "coupled")
            end
        end
    end
    M_undocking_coupled = sparse(M_undocking_coupled)
    M_dissociation = zeros(Int32,N_total, N_total)
    @time Threads.@threads for x in 1:N_total
        for y in 1:N_total
            if is_adjacent_with_remodeller(x,y) 
                M_dissociation[x,y] = dissociation_with_remodeller(x, y)
            end
        end
    end
    M_dissociation = sparse(M_dissociation)
    M_association_H2d = zeros(Int32,N_total, N_total)
    @time Threads.@threads for x in 1:N_total
        for y in 1:N_total
            if is_adjacent_with_remodeller(x,y)
                M_association_H2d[x,y] = association_with_remodeller(x, y, "H2d")
            end
        end
    end
    M_association_H2d = sparse(M_association_H2d)
    M_association_H3_H4 = zeros(Int32,N_total, N_total)
    @time Threads.@threads for x in 1:N_total
        for y in 1:N_total
            if is_adjacent_with_remodeller(x,y)
                M_association_H3_H4[x,y] = association_with_remodeller(x, y,"H3_H4")
            end
        end
    end
    M_association_H3_H4 = sparse(M_association_H3_H4)
    M_association_HEX = zeros(Int32,N_total, N_total)
    @time Threads.@threads for x in 1:N_total
        for y in 1:N_total
            if is_adjacent_with_remodeller(x,y)
                M_association_HEX[x,y] = association_with_remodeller(x, y,"HEX")
            end
        end
    end
    M_association_HEX = sparse(M_association_HEX)
    M_association_OCT = zeros(Int32,N_total, N_total)
    @time Threads.@threads for x in 1:N_total
        for y in 1:N_total
            if is_adjacent_with_remodeller(x,y)
                M_association_OCT[x,y] = association_with_remodeller(x, y,"OCT")
            end
        end
    end
    M_association_OCT = sparse(M_association_OCT)
    save("matrix_data_facilitated_multimeric.jld","helicase_on",M_helicase_on,"helicase_off",M_helicase_off,"on",M_internal_on,"off",M_internal_off,"docking_simple",M_docking_simple,"docking_H2d",M_docking_H2d,"docking_H3_H4",M_docking_H3_H4,"docking_HEX",M_docking_HEX,"undocking_simple",M_undocking_simple,"undocking_coupled",M_undocking_coupled,"association_H2d",M_association_H2d,"association_H3_H4",M_association_H3_H4,"association_HEX",M_association_HEX,"association_OCT",M_association_OCT,"dissociation_with_remodeller",M_dissociation)
end
end

function TransitionKernel(
    k_on,
    E_c,
    E_L,
    E_P_H2d,
    E_P_H3_H4,
    k_undocking_simple,
    k_undocking_coupled,
    q_on,
    q_off,
)
    # derived parameters
    k_off = k_on * exp(E_c)
    k_docking_simple = exp(-E_L) * k_undocking_simple
    k_docking_H2d = exp(-E_L-E_P_H2d) * k_undocking_coupled
    k_docking_H3_H4 = exp(-E_L-E_P_H3_H4) * k_undocking_coupled
    k_docking_HEX = exp(-2E_L - E_P_H2d - E_P_H3_H4) * k_undocking_coupled
    k_association_H2d = exp(-E_P_H2d) * k_on
    k_association_H3_H4 = exp(-E_P_H3_H4) * k_on
    k_association_HEX = exp(-E_L - E_P_H2d - E_P_H3_H4) * k_on
    k_association_OCT = exp(-2E_L - 2E_P_H2d - E_P_H3_H4) * k_on
    # collect the rates 
    rates = [
        k_on,
        k_off,
        k_docking_simple,
        k_docking_H2d,
        k_docking_H3_H4,
        k_docking_HEX,
        k_undocking_simple,
        k_undocking_coupled,
        k_association_H2d,
        k_association_H3_H4,
        k_association_HEX,
        k_association_OCT,
        q_on,
        q_off,
        k_off
    ]
    # collect the associated matrices
    matrices = [
        M_internal_on,
        M_internal_off,
        M_docking_simple,
        M_docking_H2d,
        M_docking_H3_H4,
        M_docking_HEX,
        M_undocking_simple,
        M_undocking_coupled,
        M_association_H2d,
        M_association_H3_H4,
        M_association_HEX,
        M_association_OCT,
        M_helicase_on,
        M_helicase_off,
        M_dissociation
    ]
    Q = ∑(rateᵢ * matrixᵢ for (rateᵢ, matrixᵢ) in zip(rates, matrices))
    # rectification
    for x in 1:size(Q)[1]
        Q[x,x] = - sum([Q[y,x] for y in 1:size(Q)[1] if y != x])
    end
    return Q
end

# Transition matrix for the Ω₁ states, i.e. states where there is at least one DNA-histone bond.
function TransitionKernel₁(
    k_on,
    E_c,
    E_L,
    E_P_H2d,
    E_P_H3_H4,
    k_undocking_simple,
    k_undocking_coupled,
    q_on,
    q_off,
)
   Q = TransitionKernel(
    k_on,
    E_c,
    E_L,
    E_P_H2d,
    E_P_H3_H4,
    k_undocking_simple,
    k_undocking_coupled,
    q_on,
    q_off,
    )
    Q = Q[
        [i for i ∈ eachindex(S0) if state_mapping(S0[i]) ≠ zeros(N)],
        [i for i ∈ eachindex(S0) if state_mapping(S0[i]) ≠ zeros(N)]    
    ]
    return Q
end


############# Calculation ############
id = rand()
function λ₁(θ::Par)
    nev = 1
    ncv=max(20,2*nev+1)
    which = :LR
    tol = 0.0
    maxiter = 10_000_000_000
    M = TransitionKernel₁(
        θ.k_on,
        θ.E_c,
        θ.E_L,
        θ.E_P_H2d,
        θ.E_P_H3_H4,
        θ.k_undocking_simple,
        θ.k_undocking_coupled,
        θ.q_on,
        θ.q_off,
    )
    λs,vs,nconv, niter,nmult,resid = eigs(M; 
            nev=nev, ncv=ncv, which=which, tol=tol, maxiter=maxiter)
    λs = real.(λs)
    # sort λs in descending order
    v0 = real.(vs[:,1] ./ ∑(vs[:,1]))
    @save "v0_$(id).bson" v0
    λs = sort(λs, rev=true)
    return abs.(λs[end])
end

function λ₁(θ::Par, v0::Vector{Float64})
    # v0 is the initial guess of the eigenvector
    # providing a proper v0 can speed up the computation
    nev = 1
    ncv=max(20,2*nev+1)
    which = :LR
    tol = 0.0
    maxiter = 10_000_000_000
    M = TransitionKernel₁(
        θ.k_on,
        θ.E_c,
        θ.E_L,
        θ.E_P_H2d,
        θ.E_P_H3_H4,
        θ.k_undocking_simple,
        θ.k_undocking_coupled,
        θ.q_on,
        θ.q_off,
    )
    λs,vs,nconv, niter,nmult,resid = eigs(M; 
            nev=nev, ncv=ncv, which=which, tol=tol, maxiter=maxiter, v0 = v0)
    λs = real.(λs)
    v0 = real.(vs[:,1] ./ ∑(vs[:,1]))
    @save "v0_$(id).bson" v0
    # sort λs in descending order
    λs = sort(λs, rev=true)
    return abs.(λs[end])
end

if calculate_flag
# check if the saving file exists
# if exists, load it
# and then check if the parameters are already in the file
# if not, compute the eigenvalue and save it

if isfile("model_multimeric_facilitated_disassembly_$(label).csv")
    df = CSV.read("model_multimeric_facilitated_disassembly_$(label).csv", DataFrame)
else
    df = DataFrame(
        k_on = Float64[],
        E_c = Float64[],
        E_L = Float64[],
        E_P_H2d = Float64[],
        E_P_H3_H4 = Float64[],
        k_undocking_simple = Float64[],
        k_undocking_coupled = Float64[],
        q_on = Float64[],
        q_off = Float64[],
        λ₁ = Float64[],
    )
end

if nrow(df) > 0
    calculated_pars = [
        Par(
            k_on = df.k_on[i],
            E_c = df.E_c[i],
            E_L = df.E_L[i],
            E_P_H2d = df.E_P_H2d[i],
            E_P_H3_H4 = df.E_P_H3_H4[i],
            k_undocking_simple = df.k_undocking_simple[i],
            k_undocking_coupled = df.k_undocking_coupled[i],
            q_on = df.q_on[i],
            q_off = df.q_off[i],
        )
        for i in 1:size(df, 1)
    ]
else
    calculated_pars = []
end

remaining_pars = setdiff(pars, calculated_pars)

@info "$(length(remaining_pars)) parameter(s) left to be calculated."

for par in ProgressBar(remaining_pars)
    if isfile("v0_$(id).bson")
        @load "v0_$(id).bson" v0
        push!(df, (
            k_on = par.k_on,
            E_c = par.E_c,
            E_L = par.E_L,
            E_P_H2d = par.E_P_H2d,
            E_P_H3_H4 = par.E_P_H3_H4,
            k_undocking_simple = par.k_undocking_simple,
            k_undocking_coupled = par.k_undocking_coupled,
            q_on = par.q_on,
            q_off = par.q_off,
            λ₁ = λ₁(par, v0),
        ))    
    else
        push!(df, (
            k_on = par.k_on,
            E_c = par.E_c,
            E_L = par.E_L,
            E_P_H2d = par.E_P_H2d,
            E_P_H3_H4 = par.E_P_H3_H4,
            k_undocking_simple = par.k_undocking_simple,
            k_undocking_coupled = par.k_undocking_coupled,
            q_on = par.q_on,
            q_off = par.q_off,
            λ₁ = λ₁(par),
        ))
    end
end

CSV.write("model_multimeric_facilitated_disassembly_$(label).csv", df)

end



############# Estimate ###############
using ProgressBars

function λ̂₀ₗ(k_on,E_c,q_on,q_off;N=N)
    k_off = exp(E_c) * k_on
    # steady state estimates
    if q_off == 0
        ss = k_off
    else
        ss = exp(N*E_c)*∑([
            (N-l)*(l+1)*(q_on/q_off)^l for l ∈ 0:1:N-1
            ]) / ∑([
                (l+1)*(exp(E_c)*(q_on/q_off))^l for l in 0:1:N-1]
            )
    end
    # perturbation estimates
    perturb = minimum([
        k_off,
        (N * k_on * exp(N*E_c) + q_on * ∑(
            [(k+1) * exp(k*E_c) for k ∈ 1:1:N-1]
        )) / ∑(
            [(k+1) * exp(k*E_c) for k ∈ 0:1:N]
        )
    ])
    return minimum([exp(E_c)*k_on,ss,perturb])
    # return minimum([exp(E_c)*q_on,exp(E_c)*k_on,maximum([N*(q_on/q_off)^(N-1)*exp(N*E_c),N*exp(N*E_c)])])
end

function λ̂₀ₗₛ(k_on,E_c,q_on,q_off;N=N) # one-sided perturbation
    k_off = exp(E_c) * k_on
    # steady state estimates
    if q_off == 0
        ss = k_off
    else
        ss = k_on * exp(N*E_c)*∑([
            (q_on/q_off)^l for l ∈ 0:1:N-1
            ]) / ∑([
               (exp(E_c)*(q_on/q_off))^l for l in 0:1:N-1]
            )
    end
    # perturbation estimates
    perturb = minimum([
        k_off,
        (k_on * exp(N*E_c) + q_on * ∑(
            [exp(k*E_c) for k ∈ 1:1:N-1]
        )) / ∑(
            [exp(k*E_c) for k ∈ 0:1:N]
        )
    ])
    return minimum([exp(E_c)*k_on,ss,perturb])
    # return minimum([exp(E_c)*q_on,exp(E_c)*k_on,maximum([N*(q_on/q_off)^(N-1)*exp(N*E_c),N*exp(N*E_c)])])
end

function λ̂₁(θ::Par)
    k_on = θ.k_on
    E_c = θ.E_c
    q_on = θ.q_on
    q_off = θ.q_off
    E_f_ = minimum([log(q_off/q_on),0])
    k_undocking_simple = θ.k_undocking_simple
    k_undocking_coupled = θ.k_undocking_coupled
    k_docking = k_undocking_simple * exp(-θ.E_L)
    p_docked = k_docking / (k_undocking_simple + k_docking)
    k̃_on = k_on#/(1+ exp((NL -1)*(E_c - E_f_)))
    p_reattach = k̃_on / (k̃_on + k_undocking_coupled + q_on)
    τ_L = 1/λ̂₀ₗₛ(k_on,E_c,q_on,q_off; N=NL)
    τ_M = 1/λ̂₀ₗ(k_on,E_c,q_on,q_off; N=NM)
    # if E_c - E_f_ < 0
        return maximum([minimum([(τ_L / (1 - p_docked * p_reattach) + τ_M)^(-1), k_undocking_coupled]),
        λ̂₀_linear(k_on,E_c,q_on,q_off) ])
    # else 
    #     return (τ_L / (1 - p_docked * p_reattach) + τ_M)^(-1)
    # end
end
function λ̂₀_linear(k_on,E_c,q_on,q_off)
    k_off = exp(E_c) * k_on
    # steady state estimates
    if q_off == 0
        ss = k_off
    else
        ss = exp(N*E_c)*∑([
            (N-l)*(l+1)*(q_on/q_off)^l for l ∈ 0:1:N-1
            ]) / ∑([
                (l+1)*(exp(E_c)*(q_on/q_off))^l for l in 0:1:N-1]
            )
    end
    # perturbation estimates
    perturb = minimum([
        k_off,
        (N * k_on * exp(N*E_c) + q_on * ∑(
            [(k+1) * exp(k*E_c) for k ∈ 1:1:N-1]
        )) / ∑(
            [(k+1) * exp(k*E_c) for k ∈ 0:1:N]
        )
    ])
    return minimum([exp(E_c)*k_on,ss,perturb])
    # return minimum([exp(E_c)*q_on,exp(E_c)*k_on,maximum([N*(q_on/q_off)^(N-1)*exp(N*E_c),N*exp(N*E_c)])])
end
function λ_perturb(θ::Par; N=N)
    k_on = θ.k_on
    E_c = θ.E_c
    q_on = θ.q_on
    q_off = θ.q_off
    k_off = exp(E_c) * k_on
    perturb = minimum([
            k_off,
            (N * k_on * exp(N*E_c) + q_on * ∑(
                [(k+1) * exp(k*E_c) for k ∈ 1:1:N-1]
            )) / ∑(
                [(k+1) * exp(k*E_c) for k ∈ 0:1:N]
            )
        ])
    return perturb
end
function J(p::Par)
    Ξ₁ = vcat(vec([(i,j,1) for i in 0:2, j in 0:2]),vec([(i,j,0) for i in 0:1, j in 0:1 if i + j >0]))
    function E(macrostate, microstate; p = Par())
        i,j,k = macrostate
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
        Energy = p.E_P_H2d*(P_l + P_r) + p.E_P_H3_H4*P_m + p.E_L*(L_l + L_r) + (p.E_c- minimum([0, log(p.q_off/p.q_on)]))*∑(microstate)
        return Energy
    end
    j = N * p.k_on  * (
             ∑(
                [
                    Nᵣ(ξ) / N *  exp(- E(ξ, zeros(N); p=p))
                    for ξ in Ξ₁
                ]
              )
            /∑(
                [
                    exp(- E(ξ, zeros(N); p=p) 
                            - Nᵣ(ξ)*(p.E_c- minimum([0, log(p.q_off/p.q_on)])))
                    for ξ in Ξ₁
                ]
              )
        )
    return minimum([j, maximum([
                    p.k_undocking_coupled,
                    λ_perturb(p;N=N)])])
end

λ̂_d(θ) = minimum([λ_perturb(θ; N=N),J(θ)])
# λ̂_d(θ) = J(θ)

df = DataFrame(
    k_on = Float64[],
    E_c = Float64[],
    E_L = Float64[],
    E_P_H2d = Float64[],
    E_P_H3_H4 = Float64[],
    k_undocking_simple = Float64[],
    k_undocking_coupled = Float64[],
    q_on = Float64[],
    q_off = Float64[],
    λ̂₁ = Float64[],
    λ̂_d = Float64[],
)

for par in ProgressBar(pars)
    push!(df, (
        k_on = par.k_on,
        E_c = par.E_c,
        E_L = par.E_L,
        E_P_H2d = par.E_P_H2d,
        E_P_H3_H4 = par.E_P_H3_H4,
        k_undocking_simple = par.k_undocking_simple,
        k_undocking_coupled = par.k_undocking_coupled,
        q_on = par.q_on,
        q_off = par.q_off,
        λ̂₁ = λ̂₁(par), # not used, this is only for the irreversible case
        λ̂_d = λ̂_d(par),
    ))
end

CSV.write("model_multimeric_facilitated_disassembly_$(label)_estimate.csv", df)
