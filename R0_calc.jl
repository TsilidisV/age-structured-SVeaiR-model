using GLMakie

using Integrals, QuadGK
using Dierckx, CSV, DataFrames

using FunctionWrappers
import FunctionWrappers: FunctionWrapper

# For a function that sends (x1::T1, x2::T2, ...) -> ::TN, you use the FunctionWrapper{TN, Tuple{T1, T2, ...}}.
Base.@kwdef struct SVeaiR_Param
    N0::Float64
    μ ::Float64
    prA::Float64
    prI::Float64
    c ::FunctionWrapper{Float64, Tuple{Float64}}
    βA::FunctionWrapper{Float64, Tuple{Float64}} = x -> prA * c(x) / N0
    βI::FunctionWrapper{Float64, Tuple{Float64}} = x -> prI * c(x) / N0
    p ::Float64
    ε ::Float64
    ζ ::Float64
    k ::FunctionWrapper{Float64, Tuple{Float64}}
    q ::FunctionWrapper{Float64, Tuple{Float64}}
    ξ ::FunctionWrapper{Float64, Tuple{Float64}}
    χ ::FunctionWrapper{Float64, Tuple{Float64}}
    γA::FunctionWrapper{Float64, Tuple{Float64}}
    γI::FunctionWrapper{Float64, Tuple{Float64}}
end

Base.@kwdef struct SVeaiR_ICs
    S::Float64 = 20.0*10^6
    V::Float64 = 20.0*10^6
    E::Float64
    e::Float64 = E/(ω*90)
    a::Float64 = E/(ω*90)
    i::Float64 = E/(ω*90)
end


##### Data loader #####
ω=360
function χF(θ)
    if θ/ω < 30
        return 1/5
    elseif 30 <= θ/ω <= 40
        return 1/5.8
    elseif 40 <= θ/ω <= 50
        return 1/5.8
    elseif 50 <= θ/ω <= 60
        return 1/6.5
    elseif 60 <= θ/ω <= 70
        return 1/4.1
    elseif 70 <= θ/ω
        return 1/7
    end
end

kF(θ) = χF(θ) / (1 - χF(θ))

dfAsmptPrcnt = DataFrame(CSV.File("data/dataQ.csv"))
qY = Spline1D(dfAsmptPrcnt.Age, dfAsmptPrcnt.Percentage ./ 100)
qF(θ) = qY(θ/ω)

dfContacts = DataFrame(CSV.File("data/dataContacts.csv"))
contactsY = Spline1D(dfContacts.Age, dfContacts.Percentage)
contacts(θ) = contactsY(θ/ω)


##### R_0 calculator #####
intgr(f, a, b) = solve(IntegralProblem((x, p) -> f(x), a, b), QuadGKJL(); reltol = 1e-10, abstol = 1e-10).u

# function that takes as an argument a set of parameters and spits out the R_0
function R0(P::SVeaiR_Param)
    ΕΝΔ = 90.0 * 360.0

    f1(s) = exp(-intgr(τ -> P.k(τ) + P.μ, 0.0, s)[1] )
    f2(s) = exp(-intgr(τ  -> P.γA(τ) * P.ξ(τ) + P.χ(τ) * (1 - P.ξ(τ)) + P.μ, 0.0, s)[1] )
    f3(s) = exp(-intgr(τ  -> P.γI(τ) + P.μ, 0.0, s)[1] )

    cint = intgr(s -> P.k(s) * P.q(s) * f1(s), 0.0, ΕΝΔ)[1]

    RA = cint * intgr(s -> P.βA(s)  * f2(s), 0.0, ΕΝΔ)[1]
    RI =   (quadgk(s -> P.k(s) * f1(s), 0.0, ΕΝΔ)[1]  - cint + 
              cint *
              quadgk(s -> P.χ(s) * (1 - P.ξ(s)) * f2(s), 0.0, ΕΝΔ)[1]
          ) * quadgk(s -> P.βI(s) * f3(s), 0.0, ΕΝΔ)[1]

    return ( (P.μ * P.N0) / (P.p + P.μ) ) * (1 + P.p * (1-P.ε) / (P.ε * P.ζ + P.μ)  ) * (RA + RI)
end

