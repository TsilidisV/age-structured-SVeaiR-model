@views function SVeaiR_Solver(P::SVeaiR_Param, IC::SVeaiR_ICs)
    T = 1500  # total time
    dt = 0.05  # time step
    tp = 0:dt:T  # time points
    tL = length(tp) # number of iterations

    θł = 360*90 # maximum age
    θp = 0:dt:θł  # age points
    θL = length(θp) # number of age-steps needed to reach maximum age

    # The parameteres with age are only multiplied with with variables with age.
    # So they need to be the same size as those variables, i.e. it must be length(θpleft) = length(eleft) = length(aleft) = length(ileft)
    θpleft = θp[1:end-1]

    # Parameter initialization
    N0 = P.N0
    μ = P.μ
    βA = [ P.βA(θ) for θ in θpleft ]
    βI = [ P.βI(θ) for θ in θpleft ]
    p = P.p
    ε = P.ε
    ζ = P.ζ
    k = [ P.k(θ) for θ in θpleft ]
    q = [ P.q(θ) for θ in θpleft ]
    ξ = P.ξ(1)
    χ = [ P.χ(θ) for θ in θpleft ]
    γA = P.γA(1)
    γI = P.γI(1)

    # Initial values for time = 0
    TIV_S = IC.S
    TIV_V = IC.V
    TIV_e(θ) = IC.e
    TIV_a(θ) = IC.a
    TIV_i(θ) = IC.i

    # Initial values for age = 0
    AIV_e(S,V,e,a,i) = (S + (1 - ε) * V ) * dt * sum( βA .* a .+ βI .* i )
    AIV_a(S,V,e,a,i) = dt * sum( k .* q .* e )
    AIV_i(S,V,e,a,i) = dt *sum( k .* (1 .- q) .* e .+ χ .* (1 .- ξ) .* a )

    # Reaction functions
    f_S(S,V,e,a,i) = μ * N0 - ( p + dt * sum( βA .* a .+ βI .* i ) + μ) * S
    f_V(S,V,e,a,i) = p * S - ( ζ*ε + dt * sum( βA .* a .+ βI .* i ) .* (1 .- ε) .+ μ) * V
    f_e(S,V,e,a,i) = @. - (k + μ) * e
    f_a(S,V,e,a,i) = @. - (γA * ξ + χ*(1 - ξ) + μ) * a
    f_i(S,V,e,a,i) = @. - (γI + μ) * i

    # Initialization of variables
    S = TIV_S
    V = TIV_V
    e = [ TIV_e(θ) for θ in θp ]
    a = [ TIV_a(θ) for θ in θp ]
    i = [ TIV_i(θ) for θ in θp ]
    R = N0 - S - V - dt * sum(e) - dt * sum(a) - dt * sum(i)

    # Initialization of left vectors
    Sleft = S
    Vleft = V
    eleft = e[1:end-1]
    aleft = a[1:end-1]
    ileft = i[1:end-1]
    Rleft = R

    nvis = 10 # save every nvis time-steps
    # Initilize saver
    SR = zeros( length(0:nvis:tL)  )
    VR = zeros( length(0:nvis:tL)  )
    eR = zeros( θL, length(0:nvis:tL)  )
    aR = zeros( θL, length(0:nvis:tL)  )
    iR = zeros( θL, length(0:nvis:tL)  )
    RR = zeros( length(0:nvis:tL) )


    for (ZR,Z) in zip( (eR, aR, iR), (e, a, i) )
        ZR[:,1] = Z
    end 

    for (ZR,Z) in zip( (SR, VR, RR), (S, V, R) )
        ZR[1] = Z
    end 

    j = 2 # save step counter
    for index in 2:tL
        S = Sleft + dt * f_S(Sleft, Vleft, eleft, aleft, ileft)
        V = Vleft + dt * f_V(Sleft, Vleft, eleft, aleft, ileft)
        e[2:end] = eleft .+ dt .* f_e(Sleft, Vleft, eleft, aleft, ileft)
        a[2:end] = aleft .+ dt .* f_a(Sleft, Vleft, eleft, aleft, ileft)
        i[2:end] = ileft .+ dt .* f_i(Sleft, Vleft, eleft, aleft, ileft)

        e[1] = AIV_e(Sleft, Vleft, eleft, aleft, ileft)
        a[1] = AIV_a(Sleft, Vleft, eleft, aleft, ileft)
        i[1] = AIV_i(Sleft, Vleft, eleft, aleft, ileft)

        R = N0 - S - V - dt * sum(e) - dt * sum(a) - dt * sum(i)

        # Asign the new left vectors for the next iteration
        Sleft = S
        Vleft = V
        eleft = e[1:end-1]
        aleft = a[1:end-1]
        ileft = i[1:end-1]
        Rleft = R

        if index%nvis == 0
            for (ZR,Z) in zip( (eR, aR, iR), (e, a, i) )
                ZR[:,j] = Z
            end 
            
            for (ZR,Z) in zip( (SR, VR, RR), (S, V, R) )
                ZR[j] = Z
            end 

            j += 1
        end
    end
    return SR, VR, eR, aR, iR, RR
end


