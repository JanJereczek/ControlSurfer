using DrWatson
@quickactivate "Surfer.jl"

using DifferentialEquations:ODEProblem, solve
using RobustAndOptimalControl
using NLsolve:mcpsolve
using Enzyme:jacobian, Forward
using CairoMakie
using ControlSystems
using LinearAlgebra

include(srcdir("model.jl"))
include(srcdir("generate_linear_models.jl"))

##############################################################
# Generate models linearised around [280, 600] ppm
##############################################################

function carboncylce!(du,u,p,t)
    M_A, M_U, M_D, M_L = u
    
    Emissions, Injections = p
    cumulative_emission = 0

    du[1] = Emissions(t) - k_AU*(M_A - (mA/(W_U*K0))*B(M_U)*M_U) - k_AL*(β_L*M_A_PI*(1 - M_A_PI/M_A) - (M_L - M_L_PI))
    du[2] = k_AU*(M_A - (mA/(W_U*K0))*B(M_U)*M_U) - k_UD*(M_U - M_D/(δ*δDIC))
    du[3] = k_UD*(M_U - M_D/(δ*δDIC))
    du[4] = k_AL*(β_L*M_A_PI*(1 - M_A_PI/M_A) - (M_L - M_L_PI))

end

sfmodels = [eqmodel, eqcarboncylce]
sfmodels! = [eqmodel!, eqcarboncylce!]
choose_model = 1
sfmodel = sfmodels[choose_model]
sfmodel! = sfmodels![choose_model]

# initial guess as pre-industrial equilibrium
if choose_model == 1
    x_eq_preindustrial = [
        M_A_PI,
        M_U_PI,
        M_D_PI,
        M_L_PI,
        0.0,
        0.0,
        1.0,
        1.0,
        0.0,
    ]
elseif choose_model == 2
    x_eq_preindustrial = [
        M_A_PI,
        M_U_PI,
        M_D_PI,
        M_L_PI,
    ]
    initial_co2 = sum(x_eq_preindustrial) + M_D_PI
end

nx = length(x_eq_preindustrial)
ppmCO2 = 280.0:10.0:600.0
Xeq = zeros(nx, length(ppmCO2))
Xeq[:, 1] = x_eq_preindustrial
transfer_functions = Dict{Real, TransferFunction}()
transfer_functions[ppmCO2[1]] = get_transfer_function(Xeq[:, 1], sfmodel)

for i in 2:length(ppmCO2)
    GtCO2 = ppmtoGtC(ppmCO2[i])
    Xeq[:, i] = get_equilibrium(GtCO2, vcat(GtCO2, Xeq[2:end, i-1]), sfmodel!)
    transfer_functions[ppmCO2[i]] = get_transfer_function(Xeq[:, i], sfmodel)

    fzero_check = isapprox(sfmodel(Xeq[:, i]), zeros(nx), atol = 1e-10)
    println("Linearisation for $(ppmCO2[i]) completed, zero derivative $(fzero_check)")
end

ss_tests = [linear_ss_surfer(Xeq[:, i], sfmodel) for i in eachindex(ppmCO2)]
ranks = [rank(ss_.A) for ss_ in ss_tests]

##############################################################
# Plot uncertainty from 280ppm to 600ppm
##############################################################

w_plot = 10 .^ range(-5, stop = 2, length = 500)
ppm_nominal = 450.0
G = minreal(transfer_functions[ppm_nominal])

fig = Figure(resolution = (1600, 900))
ax1 = Axis(
    fig[1, 1],
    ylabel = L"Amplitude response (dB) $\,$",
    xscale = log10,
    yscale = log10,
)
ax2 = Axis(
    fig[2, 1],
    xlabel = L"Frequency (rad/s) $\,$",
    ylabel = L"Additive uncertainty $l_A$ (dB)",
    xscale = log10,
    yscale = log10,
)

for i in eachindex(ppmCO2)
    transfer_function = transfer_functions[ppmCO2[i]]
    mag, phase, w = bode(transfer_function, w_plot)
    delta_tf = (transfer_function - G)/G
    delta_mag, delta_phase, w = bode(delta_tf, w_plot)
    line_color = (ppmCO2[i] == ppm_nominal ? :blue : :gray)
    lines!(ax1, w, mag[1, 1, :], color = line_color)
    lines!(ax2, w, delta_mag[1, 1, :] .+ 1e-8, color = :gray)
end

M_S = 1.0
A_S = 0.001
W_S_star = 1
WS = tf([1/M_S, W_S_star], [1, W_S_star * A_S])

WU = tf([1], [1])

M_T = 0.6
A_T = 0.001
W_T_star = 0.1
WT = tf([1, W_T_star / M_T], [A_T, W_T_star])     # multplicative uncertainty, weight of T

mag_ws, phase_ws, w_ws = bode(1/WS, w_plot)
mag_wu, phase_wu, w_wu = bode(1/WU, w_plot)
mag_wt, phase_wt, w_wt = bode(1/WT, w_plot)
lines!(ax2, w_plot, mag_ws[1, 1, :] .+ 1e-8, color = :blue, label = L"$w_S^{-1}$")
lines!(ax2, w_plot, mag_wu[1, 1, :] .+ 1e-8, color = :black, label = L"$w_U^{-1}$")
lines!(ax2, w_plot, mag_wt[1, 1, :] .+ 1e-8, color = :red, label = L"$w_T^{-1}$")
ylims!(ax2, (1e-4, 1.1e0))
axislegend(ax2)

save(plotsdir("model_family.png"), fig)
save(plotsdir("model_family.pdf"), fig)

# Form augmented P dynamics in state-space
P = hinfpartition(G, WS, WU, WT)
# rank(P.A - P.B2 * pinv(P.D12) * P.C1) < size(P.A, 1)

# Check that the assumptions are satisfied
flag = hinfassumptions(P)

# Synthesize the H-infinity optimal controller
C, γ = hinfsynthesize(P)
Pcl, S, CS, T = hinfsignals(P, G, C)
# Kr, hs, infor = baltrunc_coprime(C)
# n = findlast(RobustAndOptimalControl.error_bound(hs) .> 2/3) # 2/3 e sets the robustness margin
# Ksr, hs, infor = baltrunc_coprime(C; n)
# ncfmargin(G_origin, Ksr)[1] >= 2/3
controller_reduction_plot(G_origin, C)

x_perturbation = [
    M_A_PI + 10,
    M_U_PI + 10,
    M_D_PI + 100,
    M_L_PI + 1,
    0.0 + 0.1,
    0.0 + 0.1,
    1.0 - 0.1,
    1.0 - 0.1,
    0.0 + 0.1,
]

# GG = named_ss( transfer_functions[ppm_nominal], x=:xG, u=:uG, y=:yG)
# KK = named_ss( C, x=:xK, u=:uK, y=:yK)
# dif = sumblock("uC = yR - yP") # Sum node before C

# connections = [
#     :yG => :yG # Output to input
#     :uP => :uP
# ]

stp_fig = Figure()
stp_ax = Axis(stp_fig[1,1])
step_sol = step(Pcl, 1e1)
[lines!(stp_ax, step_sol.t, step_sol.y[i, :]) for i in axes(step_sol.y, 1)]
stp_fig

# t = collect(0e0:.1:1e3)
# u(x,t) = lsim(tf(C), )
# Gclosed = feedback(G_origin * C)
# y, t, x, uout = lsim(G_origin, u, t; x0 = x_perturbation - x_eq_preindustrial)
# lines(t, y[1,:])

# function controled_surfer()
#     dxu = K.A * xu + K.B * e
#     u = K.C * xu
#     dx = G.A * x + G.B * u
#     y = G.C * x
#     return y
# end
