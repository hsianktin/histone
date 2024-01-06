using CSV, DataFrames
∑ = sum
df = CSV.read(
    "model_multimeric_facilitated_disassembly_reversible.csv",
    DataFrame)

# choose a specific E_P 
E_P = 0.0
df = df[df[!,:E_P_H2d] .== E_P, :]
df = sort(df, [:q_on])

using Plots
pgfplotsx()
plot(
    df.q_on, 
    df.λ₁, 
    label="λ₁", 
    xlabel="q_on", 
    ylabel="λ₁", 
    legend=:topleft, 
    xaxis=:log, 
    yaxis=:log,
    size=(400, 400),
    st=:scatter,
)


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

# estimate by linear, facilitated model:
E_c = df.E_c[1]
k_on = df.k_on[1]
q_off = df.q_off[1]
N = 14
plot!(
    df.q_on, 
    [λ̂₀ₗ(k_on,E_c,q_on,q_off) for q_on in df.q_on], 
    label="λ̂_linear_facilitated", 
    xlabel="q_on", 
    ylabel="λ₁", 
    legend=:topleft, 
    xaxis=:log, 
    yaxis=:log,
    size=(400, 400),
    )

comp_df = CSV.read("linear_facilitated_results_14.csv", DataFrame)
comp_df = comp_df[comp_df[!,:q_off] .== 0.1, :]
comp_df.q_on ./= 10
comp_df = comp_df[comp_df[!,:q_on] .<= 1e2, :]
comp_df = comp_df[comp_df[!,:q_on] .≥ 1e-4, :]
plot!(
    comp_df.q_on, 
    comp_df.λ₀, 
    label="λ̂_linear_facilitated", 
    xlabel="q_on", 
    ylabel="λ₀", 
    legend=:topleft, 
    xaxis=:log, 
    yaxis=:log,
    size=(400, 400),
)
