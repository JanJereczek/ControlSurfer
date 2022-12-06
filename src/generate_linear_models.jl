using Enzyme:jacobian, Forward
using ControlSystems:ss, tf, TransferFunction
using NLsolve:mcpsolve
include("model.jl")

function get_transfer_function(
    x_eq::Vector{T},
    f::Function,
) where {T<:Real}
    return tf(linear_ss_surfer(x_eq, f))
end

function linear_ss_surfer(
    x_eq::Vector{T},
    f::Function,
) where {T<:Real}
    filler_length = length(x_eq) - 1
    A = jacobian(Forward, f, x_eq)
    # Generate small perturbation of the dep ocean dynamics
    A[3, 2] += 1e-5 * rand()
    A[3, 3] += 1e-5 * rand()
    B = vcat([1.0], zeros(filler_length))
    C = vcat([1.0], zeros(filler_length))'
    D = [0.0]
    check_nan = sum(isnan.(A))
    if check_nan > 0
        println("Detected $check_nan NaN entries")
    end
    return ss(A,B,C,D)
end

sources = [t -> 0.0, t -> 0.0]      # [co2rate, so2rate]
function eqmodel!(dx, x)
    model!(dx, x, sources, 0.0)
end

function eqmodel(x)
    dx = similar(x)
    eqmodel!(dx, x)
    return dx
end

function eqcarboncylce!(dx, x)
    carboncylce!(dx, x, sources, 0.0)
end

function eqcarboncylce(x)
    dx = similar(x)
    eqcarboncylce!(dx, x)
    return dx
end

lower_bound = [ ppmtoGtC(150), 1e3, 30e3, 2e3, -1e1, -1e1, 0e0, 0e0, -Inf ]
upper_bound = [ ppmtoGtC(1000), 2e3, 40e3, 3e3, 1e1, 1e1, 2e0, 2e0, Inf ]

# lower_bound = [ ppmtoGtC(150), 1e3, 2e3, -1e1, -1e1, 0e0, 0e0, -Inf ]
# upper_bound = [ ppmtoGtC(1000), 2e3, 3e3, 1e1, 1e1, 2e0, 2e0, Inf ]

function get_equilibrium(
    GtCO2::T,
    guess::Vector{T},
    f!::Function,
) where {T<:Real}
    n = length(guess)
    lb = lower_bound[2:n]
    ub = upper_bound[2:n]
    x_eq = mcpsolve(
        f!,
        vcat(GtCO2-1e-8, lb),   # tight bounds on CO2 --> prescribed
        vcat(GtCO2+1e-8, ub),
        guess,
        reformulation = :minmax,
        # autodiff = :finite,
        ftol = 1e-11,
        iterations = 100_000,
    )
    println(x_eq.iterations)
    return x_eq.zero
end