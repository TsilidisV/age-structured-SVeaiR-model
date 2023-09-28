using FileIO, JLD2

function saveToJLD2(M, S::String)
    result = Dict("S" => M[1], "V" => M[2], "R" => M[6],
            "E" => 0.05 * dropdims(sum(M[3], dims = 1), dims = 1),
            "A" => 0.05 * dropdims(sum(M[4], dims = 1), dims = 1), 
            "I" => 0.05 * dropdims(sum(M[5], dims = 1), dims = 1))
    JLD2.save("output/$(S).jld2", result; compress = true)
end

#R0 > 1
gC(x) = (16.71291/0.38)*exp(-((x-(ω*10)) / 10000 )^2)
sum(gC.(0:1:ω*90))/length(0:1:ω*90)

P1 = SVeaiR_Param(N0 = 80.0*10^6, μ = (16*10^(-3))/365, prA = 2/10, prI = 2/5, c = θ -> gC(θ), 
                    p = 10^-3, ε = 0.7, ζ = 1/14, k = θ -> kF(θ), q = θ -> qF(θ), ξ = θ -> 0.5,
                    χ = θ -> χF(θ), γA = θ -> 1/8, γI = θ -> 1/14)
R0(P1)
IC1(i,j) = SVeaiR_ICs(E=j*10.0^i)

sol1_1 = SVeaiR_Solver(P1,IC1(1,1))
saveToJLD2(sol1_1,"P1_IC1-1")

sol1_2 = SVeaiR_Solver(P1,IC1(2,1))
saveToJLD2(sol1_2,"P1_IC1-2")

sol1_3 = SVeaiR_Solver(P1,IC1(3,1))
saveToJLD2(sol1_3,"P1_IC1-3")

sol1_4 = SVeaiR_Solver(P1,IC1(4,1))
saveToJLD2(sol1_4,"P1_IC1-4")

sol1_5 = SVeaiR_Solver(P1,IC1(5,1))
saveToJLD2(sol1_5,"P1_IC1-5")

sol1_6 = SVeaiR_Solver(P1,IC1(6)) 
saveToJLD2(sol1_6,"P1_IC1-6")

sol1_6_2 = SVeaiR_Solver(P1,IC1(6,2)) 
saveToJLD2(sol1_6_2,"P1_IC1-6_2")

sol1_6_3 = SVeaiR_Solver(P1,IC1(6,3))
saveToJLD2(sol1_6_3,"P1_IC1-6_3")

sol1_6_4 = SVeaiR_Solver(P1,IC1(6,4))
saveToJLD2(sol1_6_4,"P1_IC1-6_4")

sol1_7_1 = SVeaiR_Solver(P1,IC1(7,1))
saveToJLD2(sol1_7_1,"P1_IC1-7_1")

sol1_7_2 = SVeaiR_Solver(P1,IC1(7,2))
saveToJLD2(sol1_7_2,"P1_IC1-7_2")

#R0 < 1
gC2(x) = exp(-((x-(ω*80)) / 10000 )^2)
sum(gC2.(0:1:ω*90))/length(0:1:ω*90)

P2 = SVeaiR_Param(N0 = 80.0*10^6, μ = (16*10^(-3))/365, prA = 2/10, prI = 2/5, c = θ -> gC2(θ), 
                    p = 10^-3, ε = 0.7, ζ = 1/14, k = θ -> kF(θ), q = θ -> qF(θ), ξ = θ -> 0.5,
                    χ = θ -> χF(θ), γA = θ -> 1/8, γI = θ -> 1/14)
R0(P2)
IC1(i,j) = SVeaiR_ICs(E=j*10.0^i)

sol2_1_1 = SVeaiR_Solver(P2,IC1(1,1))
saveToJLD2(sol2_1_1,"P2_IC1-1_1")

sol2_4_1 = SVeaiR_Solver(P2,IC1(4,1))
saveToJLD2(sol2_4_1,"P2_IC1-4_1")

sol2_6_1 = SVeaiR_Solver(P2,IC1(6,1))
saveToJLD2(sol2_6_1,"P2_IC1-6_1")

sol2_6_4 = SVeaiR_Solver(P2,IC1(6,4))
saveToJLD2(sol2_6_4,"P2_IC1-6_4")

sol2_7_1 = SVeaiR_Solver(P2,IC1(7,1))
saveToJLD2(sol2_7_1,"P2_IC1-7_1")


