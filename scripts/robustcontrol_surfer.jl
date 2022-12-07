using DrWatson
@quickactivate "Surfer.jl"

using DifferentialEquations
using RobustAndOptimalControl
using NLsolve
using Enzyme
using CairoMakie
using ControlSystems
using LinearAlgebra
using Test

include(srcdir("model.jl"))
include(srcdir("generate_linear_models.jl"))

##############################################################
# Generate models linearised around [280, 600] ppm
##############################################################

sfmodel = eqmodel
sfmodel! = eqmodel!
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
@test sum(ranks .== length(x_eq_preindustrial) ) .== length(ranks)

##############################################################
# Plot uncertainty from 280ppm to 600ppm
##############################################################

w_plot = 10 .^ range(-5, stop = 2, length = 500)
ppm_nominal = 280.0
G_nominal = minreal(transfer_functions[ppm_nominal])

logticks_val = -4:4
ticksval = 10.0 .^ logticks_val
ticksstr = [L"$10^{%$(logtick_val)}$" for logtick_val in logticks_val]
power_ticks = (ticksval, ticksstr)

fig = Figure(resolution = (1600, 900))
ax1 = Axis(
    fig[1, 1],
    ylabel = L"Amplitude response (dB) $\,$",
    xscale = log10,
    yscale = log10,
    xticks = power_ticks,
    yticks = power_ticks,
    xminorticks = IntervalsBetween(9),
    xminorgridvisible = true,
)
ax2 = Axis(
    fig[2, 1],
    xlabel = L"Frequency (rad/s) $\,$",
    ylabel = L"Additive uncertainty $l_A$ (dB)",
    xscale = log10,
    yscale = log10,
    xticks = power_ticks,
    yticks = power_ticks,
    xminorticks = IntervalsBetween(9),
    xminorgridvisible = true,
)

for i in eachindex(ppmCO2)
    transfer_function = transfer_functions[ppmCO2[i]]
    mag, phase, w = bode(transfer_function, w_plot)
    delta_tf = (transfer_function - G_nominal)/G_nominal
    delta_mag, delta_phase, w = bode(delta_tf, w_plot)
    line_color = (ppmCO2[i] == ppm_nominal ? :blue : :gray)
    lines!(ax1, w, mag[1, 1, :], color = line_color)
    lines!(ax2, w, delta_mag[1, 1, :] .+ 1e-8, color = :gray)
end

M_S = 1.0
A_S = 0.0001
W_S_star = .01
WS = tf([1/M_S, W_S_star], [1, W_S_star * A_S])

WU = tf([1], [1])

if ppm_nominal == 450
    M_T = 0.6
    A_T = 0.001
    W_T_star = 0.1
elseif ppm_nominal == 280
    M_T = 2.0
    A_T = 0.001
    W_T_star = 0.5
end
WT = tf([1, W_T_star / M_T], [A_T, W_T_star])     # multplicative uncertainty, weight of T

mag_ws, phase_ws, w_ws = bode(1/WS, w_plot)
mag_wu, phase_wu, w_wu = bode(1/WU, w_plot)
mag_wt, phase_wt, w_wt = bode(1/WT, w_plot)
lines!(ax2, w_plot, mag_ws[1, 1, :] .+ 1e-8, color = :royalblue, label = L"$w_S^{-1}$")
lines!(ax2, w_plot, mag_wu[1, 1, :] .+ 1e-8, color = :mediumpurple4, label = L"$w_U^{-1}$")
lines!(ax2, w_plot, mag_wt[1, 1, :] .+ 1e-8, color = :darkorange, label = L"$w_T^{-1}$")
ylims!(ax2, (1e-3, 2e0))
axislegend(ax2)

save(plotsdir("model_family.png"), fig)
save(plotsdir("model_family.pdf"), fig)

##############################################################
# Robust control design and reduction
##############################################################

P = hinfpartition(G_nominal, WS, WU, WT)    # Form augmented P dynamics in state-space
flag = hinfassumptions(P)           # Check that the assumptions are satisfied
K, γ = hinfsynthesize(P)            # Synthesize the H-infinity optimal controller
Pcl, S, CS, T = hinfsignals(P, G_nominal, K)

x_perturbation = [10, 0, 0, 0, 0, 0, 0, 0, 0]
yw = ppmtoGtC(ppm_nominal)
xu0 = zeros(length(K.B))
p = Dict{String, Any}()
p["yw"] = yw
p["K"] = K

perturbation_case = 2

if perturbation_case == 1
    z0 = vcat(x_eq_preindustrial + x_perturbation, xu0)
elseif perturbation_case == 2
    z0 = vcat(Xeq[:, 8], xu0)
end

include(srcdir("model.jl"))
prob = ODEProblem(controlled_model, z0, (0, 1e3), p)
@time ctrl_sol = solve(prob, Rosenbrock23(), reltol=1e-3, abstol=1e-12)
Xctrl = hcat(ctrl_sol.u...)

auto_sol = solve(prob, Rosenbrock23(), reltol=1e-3, abstol=1e-12)
Xauto = hcat(auto_sol.u...)

ylabels = [
    "Atmospheric CO2 (GtC)",
    "Upper ocean CO2 (GtC)",
    "Deep ocean CO2 (GtC)",
    "Land CO2 (GtC)",
    "Upper ocean temperature anomaly (K)",
    "Deep ocean temperature anomaly (K)",
    "GrIS volume (a.u.)",
    "AIS volume (a.u.)",
    "Sea-level change (m)",
]
time_resp_fig = Figure( resolution = (1600, 900) )
nrows, ncols = 3, 3
for i in 1:nrows
    for j in 1:ncols
        k = (i-1)*ncols + j
        time_resp_ax = Axis(
            time_resp_fig[i,j],
            xlabel = i == 3 ? "Time (yr)" : " ",
            ylabel = ylabels[k],
        )
        lines!(time_resp_ax, auto_sol.t, Xauto[k, :], label = "autonomous")
        lines!(time_resp_ax, ctrl_sol.t, Xctrl[k, :], label = "controlled")
    end
end
time_resp_fig
figname = "tresponse_ppmnominal$(ppm_nominal)_perturbation$(perturbation_case)"
save(plotsdir(string(figname, ".png")), time_resp_fig)
save(plotsdir(string(figname, ".pdf")), time_resp_fig)

# function emission(x::Vector{T}, xu::Vector{T}, K::StateSpace{Continuous, T}) where {T<:Real}
#     e = x[1] - yw
#     dxu = K.A * xu + K.B * e
#     xu += dxu * delta_t
#     u = K.C * xu
#     return xu, u[1]
# end
# xu, U = emission(x_eq_preindustrial, xu0, K)

# Kr, hs, infor = baltrunc_coprime(C)
# n = findlast(RobustAndOptimalControl.error_bound(hs) .> 2/3) # 2/3 e sets the robustness margin
# Ksr, hs, infor = baltrunc_coprime(C; n)
# ncfmargin(G_origin, Ksr)[1] >= 2/3
# controller_reduction_plot(G_origin, C)
