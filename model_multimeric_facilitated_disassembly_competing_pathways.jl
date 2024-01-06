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
    default = true
    arg_type = Bool
    required = false
end

profile = parse_args(ARGS, s)["profile"]
calculate_flag = true
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
    # first identify adjacent states
    @time atlas = [(x,y) for x in 1:N_total, y in 1:N_total 
                    if is_adjacent_with_remodeller(x,y)]

    M_helicase_on = zeros(Int32,N_total, N_total)
    @time Threads.@threads for (x,y) ∈ atlas
                M_helicase_on[x,y] = helicase_on(x, y)
    end
    M_helicase_on = sparse(M_helicase_on)
    M_helicase_off = zeros(Int32,N_total, N_total)
    @time Threads.@threads for (x,y) ∈ atlas
                M_helicase_off[x,y] = helicase_off(x, y)
    end
    M_helicase_off = sparse(M_helicase_off)
    M_internal_on = zeros(Int32,N_total, N_total)
    @time Threads.@threads for (x,y) ∈ atlas
                M_internal_on[x,y] = internal_on_with_remodeller(x, y)
    end
    M_internal_on = sparse(M_internal_on)
    M_internal_off = zeros(Int32,N_total, N_total)
    @time Threads.@threads for (x,y) ∈ atlas
                M_internal_off[x,y] = internal_off_with_remodeller(x, y)
    end
    M_internal_off = sparse(M_internal_off)
    M_docking_simple = zeros(Int32,N_total, N_total)
    @time Threads.@threads for (x,y) ∈ atlas
                M_docking_simple[x,y] = docking_with_remodeller(x, y,"simple")
    end
    M_docking_simple = sparse(M_docking_simple)
    M_docking_H2d = zeros(Int32,N_total, N_total)
    @time Threads.@threads for (x,y) ∈ atlas
                M_docking_H2d[x,y] = docking_with_remodeller(x, y,"H2d")
    end
    M_docking_H2d = sparse(M_docking_H2d)
    M_docking_H3_H4 = zeros(Int32,N_total, N_total)
    @time Threads.@threads for (x,y) ∈ atlas
                M_docking_H3_H4[x,y] = docking_with_remodeller(x, y,"H3_H4")
    end
    M_docking_H3_H4 = sparse(M_docking_H3_H4)
    M_docking_HEX = zeros(Int32,N_total, N_total)
    @time Threads.@threads for (x,y) ∈ atlas
                M_docking_HEX[x,y] = docking_with_remodeller(x, y,"HEX")
    end
    M_docking_HEX = sparse(M_docking_HEX)
    M_undocking_simple = zeros(Int32,N_total, N_total)
    @time Threads.@threads for (x,y) ∈ atlas
                M_undocking_simple[x,y] = undocking_with_remodeller(x, y,"simple")
    end
    M_undocking_simple = sparse(M_undocking_simple)
    M_undocking_coupled = zeros(Int32,N_total, N_total)
    @time Threads.@threads for (x,y) ∈ atlas
                M_undocking_coupled[x,y] = undocking_with_remodeller(x, y, "coupled")
    end
    M_undocking_coupled = sparse(M_undocking_coupled)
    M_dissociation = zeros(Int32,N_total, N_total)
    @time Threads.@threads for (x,y)  ∈ atlas
                M_dissociation[x,y] = dissociation_with_remodeller(x, y)
    end
    M_dissociation = sparse(M_dissociation)
    M_association_H2d = zeros(Int32,N_total, N_total)
    @time Threads.@threads for (x,y) ∈ atlas
                M_association_H2d[x,y] = association_with_remodeller(x, y, "H2d")
    end
    M_association_H2d = sparse(M_association_H2d)
    M_association_H3_H4 = zeros(Int32,N_total, N_total)
    @time Threads.@threads for (x,y) ∈ atlas
                M_association_H3_H4[x,y] = association_with_remodeller(x, y,"H3_H4")
    end
    M_association_H3_H4 = sparse(M_association_H3_H4)
    M_association_HEX = zeros(Int32,N_total, N_total)
    @time Threads.@threads for (x,y) ∈ atlas
                M_association_HEX[x,y] = association_with_remodeller(x, y,"HEX")
    end
    M_association_HEX = sparse(M_association_HEX)
    M_association_OCT = zeros(Int32,N_total, N_total)
    @time Threads.@threads for (x,y) ∈ atlas
                M_association_OCT[x,y] = association_with_remodeller(x, y,"OCT")
    end
    M_association_OCT = sparse(M_association_OCT)
    save("matrix_data_facilitated_multimeric.jld","helicase_on",M_helicase_on,"helicase_off",M_helicase_off,"on",M_internal_on,"off",M_internal_off,"docking_simple",M_docking_simple,"docking_H2d",M_docking_H2d,"docking_H3_H4",M_docking_H3_H4,"docking_HEX",M_docking_HEX,"undocking_simple",M_undocking_simple,"undocking_coupled",M_undocking_coupled,"association_H2d",M_association_H2d,"association_H3_H4",M_association_H3_H4,"association_HEX",M_association_HEX,"association_OCT",M_association_OCT,"dissociation_with_remodeller",M_dissociation)
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

# function Alt_TransitionKernel(
#     θ
# )

############# Calculation ############
function Pₐ(θ::Par)
    𝐐 = TransitionKernel(
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
    ∂Ω = [i for i in eachindex(S0) if (state_mapping_with_remodeller(S0[i])[1] |> ∑ == 1)]
    Sₐ = [S0[i][1:3] for i in eachindex(S0) if (state_mapping_with_remodeller(S0[i])[1] |> ∑ == 1)] |> unique
    
    Ωᵢₙ = [i for i ∈ eachindex(S0) if state_mapping(S0[i]) ≠ zeros(N)]
    x₀ = findfirst(x -> S0[x] == [2,2,1,[0,0],[0,0]], Ωᵢₙ)
    ∂Ωᵢₙ = [i for i ∈ eachindex(Ωᵢₙ) if Ωᵢₙ[i] ∈ ∂Ω] # coordiante ∂Ω in terms of Ωᵢₙ
    𝐐ᵢₙ_ᵢₙ = 𝐐[Ωᵢₙ, Ωᵢₙ]
    𝐐ₐ_ᵢₙ = zeros(length(Sₐ), length(Ωᵢₙ))
    for i in ∂Ωᵢₙ
        j = findfirst(x -> x == S0[Ωᵢₙ[i]][1:3], Sₐ)
        𝐐ₐ_ᵢₙ[j, i] = - ∑(𝐐ᵢₙ_ᵢₙ[:, i])
    end
    # sparse matrix
    𝐐ₐ_ᵢₙ = sparse(𝐐ₐ_ᵢₙ)
    Sₐᵣ = [i for i in eachindex(Sₐ) if ∑(Sₐ[i] .> 0) == 3]
    𝐐ₐᵣ_ᵢₙ = 𝐐ₐ_ᵢₙ[Sₐᵣ, :]
    𝟏 = ones(length(Sₐᵣ))
    Pₐ_ᵢₙ = - 𝐐ᵢₙ_ᵢₙ' \ (𝐐ₐᵣ_ᵢₙ' * 𝟏)
    # 𝐐ᵢₙ_ᵢₙ' * ones(length(Ωᵢₙ)) + 𝐐ₐ_ᵢₙ' * ones(length(Sₐ)) .≈ 0
    return Pₐ_ᵢₙ[x₀]
end    


# check if the saving file exists
# if exists, load it
# and then check if the parameters are already in the file
# if not, compute the eigenvalue and save it

if isfile("corrected_model_multimeric_facilitated_disassembly_competing_pathway_$(label).csv")
    df = CSV.read("corrected_model_multimeric_facilitated_disassembly_competing_pathway_$(label).csv", DataFrame)
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
        Pₐₗₓ₀ = Float64[],
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
            Pₐₗₓ₀ = Pₐ(par),
        ))    
end

CSV.write("corrected_model_multimeric_facilitated_disassembly_competing_pathway_$(label).csv", df)

