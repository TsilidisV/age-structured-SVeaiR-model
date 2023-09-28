using CairoMakie, CodecZlib, JLD2, LaTeXStrings

xtickshort = (0:50:200, [L"0", L"25", L"50", L"75", L"100"])
xticklong = ( [0, 1400, 2800], [L"100", L"800", L"1500"] )

function powerOfTen(num)
    power = 0

    if 0 < num < 1
        num = round(num; digits = 2)
        return latexstring("$(num)")
    end

    rndnum = round(Int, num)
    while rndnum % 10 == 0 && rndnum != 0
        rndnum = rndnum / 10
        power += 1
    end

    rndnum = round(Int, rndnum)
    num = round(Int, num)

    if power <= 3 
        return latexstring("$(num)")
    else
        return latexstring(latexstring("$(rndnum)"),L"\cdot",latexstring("10^$(power)"))
    end
end

total_theme = Theme( fontsize = 12, figure_padding = 5,
    Axis = (
        titlesize = 13,
        leftspinevisible = false,
        rightspinevisible = false,
        bottomspinevisible = false,
        topspinevisible = false,
        backgroundcolor = "#E5ECF6",
        xgridcolor = :white,
        ygridcolor = :white,
        xtickcolor = "#E5ECF6",
        ytickcolor = "#E5ECF6",
        yticklabelcolor = "#2A3F5F",
        xticklabelcolor = "#2A3F5F",
        titlecolor = "#0B2D60",
        xlabelcolor = "#0B2D60",
        ylabelcolor = "#0B2D60"
    ),
    palette = (color = ["#F32EA1", "#FFD92E", "#8630BF", "#0BD976", "#FF766E"],
               linestyle = [nothing, :dash, :dot, :dashdot, :dashdotdot],
               marker = [:star5, :diamond, :diamond, :hexagon]
    ),
    Lines = (
        linewidth = 2.0,
        markersize = 7,
        cycle = Cycle([:color,:linestyle], covary = true)
    ),
    Legend = (
        labelcolor = "#2A3F5F",
        titlecolor = "#0B2D60",
        framevisible = false,
        margin = (0,0,-25,-10),
        orientation = :horizontal,
        tellheight = true,
        tellwidth = false 
    )
)

function add_axis_inset(pos=fig[1, 1]; bgcolor=:white, halign, valign, width=Relative(0.6), height=Relative(0.6), alignmode=Mixed(top=5, right=8))
    inset_box = Axis(pos; width, height, halign, valign, alignmode,
        xticklabelsize=10, yticklabelsize=10, backgroundcolor=bgcolor, xticks = xticklong,
        xgridcolor = "#E5ECF6", ygridcolor = "#E5ECF6")
    # bring content upfront
    translate!(inset_box.scene, 0, 0, 10),
    return inset_box
end

function figure_R0_more1()
    loads = ["1", "4", "6", "6_4", "7_1"]
    labels = [L"$d$ = 10", L"$d$ = 10$^4$", L"$d$ = 10$^6$", L"$d$ = 4$\cdot$10$^6$", L"$d$ = 10$^7$"]

    size_inches = (8, 11.5)
    size_pt = 72 .* size_inches
  
    fig = Figure( resolution = size_pt
             )


    axS = Axis(fig[1,1], xticks = xtickshort, titlesize = 12, xlabel = "Days", yticks = ([0, 5*10^6, 1*10^7, 1.5*10^7, 2*10^7], [L"0", L"5\cdot10^6", L"1\cdot10^7", L"1.5\cdot10^7", L"2\cdot10^7"]), title = "Susceptible"            )
    axV = Axis(fig[1,2], xticks = xtickshort, titlesize = 12, xlabel = "Days", yticks = ([0, 5*10^6, 1*10^7, 1.5*10^7, 2*10^7], [L"0", L"5\cdot10^6", L"1\cdot10^7", L"1.5\cdot10^7", L"2\cdot10^7"]), title = "Vaccinated"              )
    axE = Axis(fig[2,1], xticks = xtickshort, titlesize = 12, xlabel = "Days", yticks = ([0, 1*10^7, 2*10^7, 3*10^7], [L"0", L"1\cdot10^7", L"2\cdot10^7", L"3\cdot10^7"]), title = "Exposed"                 )
    axA = Axis(fig[2,2], xticks = xtickshort, titlesize = 12, xlabel = "Days", yticks = ([0, 5*10^6, 1*10^7], [L"0", L"5\cdot10^6", L"1\cdot10^7"]), title = "Asymptomatic infectious" )
    axI = Axis(fig[3,1], xticks = xtickshort, titlesize = 12, xlabel = "Days", yticks = ([0, 1*10^7, 2*10^7, 3*10^7], [L"0", L"1\cdot10^7", L"2\cdot10^7", L"3\cdot10^7"]), title = "Symptomatic infectious"  )
    axR = Axis(fig[3,2], xticks = xtickshort, titlesize = 12, xlabel = "Days", ytickformat = y -> powerOfTen.(y), title = "Removed"                 )

    axSs = add_axis_inset(fig[1,1]; halign=:right, valign=:top)
    axSs.ytickformat = y -> powerOfTen.(y)
    axVs = add_axis_inset(fig[1,2]; halign=:right, valign=:top)
    axVs.ytickformat = y -> powerOfTen.(y)
    axEs = add_axis_inset(fig[2,1]; halign=:right, valign=:top)
    axEs.ytickformat = y -> powerOfTen.(y)
    axAs = add_axis_inset(fig[2,2]; halign=:right, valign=:top)
    axAs.ytickformat = y -> powerOfTen.(y)
    axIs = add_axis_inset(fig[3,1]; halign=:right, valign=:top)
    axIs.ytickformat = y -> powerOfTen.(y)
    axRs = add_axis_inset(fig[3,2]; halign=:right, valign=:bottom, alignmode=Mixed(bottom=5, right=8, top = 0))
    axRs.yticks = [794 * 10^5, 796 * 10^5, 798 * 10^5]
    axRs.ytickformat = y -> powerOfTen.(y)

    for i in 1:length(loads)
        scS = lines!(axS, load("output/P1_IC1-$(loads[i]).jld2")["S"][1:200], label = labels[i] )
        scV = lines!(axV, load("output/P1_IC1-$(loads[i]).jld2")["V"][1:200] )
        scE = lines!(axE, load("output/P1_IC1-$(loads[i]).jld2")["E"][1:200] )
        scA = lines!(axA, load("output/P1_IC1-$(loads[i]).jld2")["A"][1:200] )
        scI = lines!(axI, load("output/P1_IC1-$(loads[i]).jld2")["I"][1:200] )
        scR = lines!(axR, load("output/P1_IC1-$(loads[i]).jld2")["R"][1:200] )

        scSs = lines!(axSs, load("output/P1_IC1-$(loads[i]).jld2")["S"][200:end] )
        scVs = lines!(axVs, load("output/P1_IC1-$(loads[i]).jld2")["V"][200:end] )
        scEs = lines!(axEs, load("output/P1_IC1-$(loads[i]).jld2")["E"][200:end] )
        scAs = lines!(axAs, load("output/P1_IC1-$(loads[i]).jld2")["A"][200:end] )
        scIs = lines!(axIs, load("output/P1_IC1-$(loads[i]).jld2")["I"][200:end] )
        scRs = lines!(axRs, load("output/P1_IC1-$(loads[i]).jld2")["R"][200:end] )
    end

    Legend(fig[0,:], axS,
        L"Solution of the problem for $\mathcal{R}_0 > 1$ and different initial values of $E$, $A$ and $I$", orientation = :horizontal
    )

    fig
end

function figure_R0_leq1()
    loads = ["1_1", "4_1", "6_1", "6_4", "7_1"]
    labels = [L"$d$ = 10", L"$d$ = 10$^4$", L"$d$ = 10$^6$", L"$d$ = 4$\cdot$10$^6$", L"$d$ = 10$^7$"]

    size_inches = (8, 11.5)
    size_pt = 72 .* size_inches
  
    fig = Figure( resolution = size_pt
             )

    axS = Axis(fig[1,1], title = "Susceptible",             xticks = xtickshort, titlesize = 12, xlabel = "Days", ytickformat = y -> powerOfTen.(y))
    axV = Axis(fig[1,2], title = "Vaccinated",              xticks = xtickshort, titlesize = 12, xlabel = "Days", ytickformat = y -> powerOfTen.(y))
    axE = Axis(fig[2,1], title = "Exposed",                 xticks = xtickshort, titlesize = 12, xlabel = "Days", ytickformat = y -> powerOfTen.(y))
    axA = Axis(fig[2,2], title = "Asymptomatic infectious", xticks = xtickshort, titlesize = 12, xlabel = "Days", ytickformat = y -> powerOfTen.(y))
    axI = Axis(fig[3,1], title = "Symptomatic infectious",  xticks = xtickshort, titlesize = 12, xlabel = "Days", ytickformat = y -> powerOfTen.(y))
    axR = Axis(fig[3,2], title = "Removed",                 xticks = xtickshort, titlesize = 12, xlabel = "Days", ytickformat = y -> powerOfTen.(y))

    axSs = add_axis_inset(fig[1,1]; halign=:right, valign=:top, alignmode=Mixed(right=8, top = 25))
    axSs.ytickformat = y -> powerOfTen.(y)
    axVs = add_axis_inset(fig[1,2]; halign=:right, valign=:top)
    axVs.ytickformat = y -> powerOfTen.(y)
    axEs = add_axis_inset(fig[2,1]; halign=:right, valign=:top)
    axEs.ytickformat = y -> powerOfTen.(y)
    axAs = add_axis_inset(fig[2,2]; halign=:right, valign=:top)
    axAs.ytickformat = y -> powerOfTen.(y)
    axIs = add_axis_inset(fig[3,1]; halign=:right, valign=:top, width=Relative(0.5))
    axIs.ytickformat = y -> powerOfTen.(y)
    axRs = add_axis_inset(fig[3,2]; halign=:right, valign=:bottom, alignmode=Mixed(bottom=5, right=8, top = 0))
    axRs.ytickformat = y -> powerOfTen.(y)

    for i in 1:length(loads)
        scS = lines!(axS, load("output/P2_IC1-$(loads[i]).jld2")["S"][1:200], label = labels[i] )
        scV = lines!(axV, load("output/P2_IC1-$(loads[i]).jld2")["V"][1:200] )
        scE = lines!(axE, load("output/P2_IC1-$(loads[i]).jld2")["E"][1:200] )
        scA = lines!(axA, load("output/P2_IC1-$(loads[i]).jld2")["A"][1:200] )
        scI = lines!(axI, load("output/P2_IC1-$(loads[i]).jld2")["I"][1:200] )
        scR = lines!(axR, load("output/P2_IC1-$(loads[i]).jld2")["R"][1:200] )

        scSs = lines!(axSs, load("output/P2_IC1-$(loads[i]).jld2")["S"][200:end] )
        scVs = lines!(axVs, load("output/P2_IC1-$(loads[i]).jld2")["V"][200:end] )
        scEs = lines!(axEs, load("output/P2_IC1-$(loads[i]).jld2")["E"][200:end] )
        scAs = lines!(axAs, load("output/P2_IC1-$(loads[i]).jld2")["A"][200:end] )
        scIs = lines!(axIs, load("output/P2_IC1-$(loads[i]).jld2")["I"][200:end] )
        scRs = lines!(axRs, load("output/P2_IC1-$(loads[i]).jld2")["R"][200:end] )
    end

    Legend(fig[0,:], axS,
        L"Solution of the problem for $\mathcal{R}_0 \le 1$ and different initial values of $E$, $A$ and $I$", orientation = :horizontal
    )

    fig
end

figure_R0_more1()
with_theme(figure_R0_more1, total_theme)


figure_R0_leq1()
with_theme(figure_R0_leq1, total_theme)

save("figure2.pdf", with_theme(figure_R0_leq1, total_theme), pt_per_unit = 1)




function contacts_figure()

    c_1(x) = (16.71/0.37)*exp(-((x-(360.0*10)) / 10000 )^2)
    c_2(x) = (16.71/0.38)*exp(-((x-(360*80)) / 10000 )^2)

    size_inches = (8, 2.5)
    size_pt = 72 .* size_inches
  
    fig = Figure( resolution = size_pt
             )

    ax = Axis(fig[1,1], xlabel = L"$\theta \; (\cdot 360$ days)", 
        xticks = (0:(360*10):360*90, powerOfTen.(0:10:90)),
        ylabel = L"Transmission rate $\; \left(\cdot \frac{1}{N_0} \right)$",
        ytickformat = y -> powerOfTen.(y)
        )

    θs = 0:1:360*90
    lines!(ax, θs, c_1.(θs)*(2/10), color = "#703791", label = L"\beta_{A,1}")
    lines!(ax, θs, c_1.(θs)*(2/5) , color = "#338F85", label = L"\beta_{I,1}")
    lines!(ax, θs, c_2.(θs)*(2/10), color = "#703791", linestyle = :dash, label = L"\beta_{A,2}")
    lines!(ax, θs, c_2.(θs)*(2/5) , color = "#338F85", linestyle = :dash, label = L"\beta_{I,2}")

    Legend(fig[0,:], ax)
    fig
end

save("figure2.pdf", with_theme(contacts_figure, total_theme), pt_per_unit = 1)

function asymptomatics_figure()
    size_inches = (8, 2.5)
    size_pt = 72 .* size_inches
  
    fig = Figure( resolution = size_pt
             )

    ax = Axis(fig[1,1], xlabel = L"$\theta \; (\cdot 360$ days)", 
        xticks = (0:(360*10):360*90, powerOfTen.(0:10:90)),
        ylabel = L"q(\theta) \; (%)",
        title = L"\text{Percentage of exposed individuals becoming asymptomatic infectious}",
        ytickformat = y -> powerOfTen.(100*y)
        )

    θs = 0:.5:360*90
    lines!(ax, θs, qF.(θs), color = "#F77474", label = L"\beta_{A,1}")

    fig
end

save("figure3.pdf", with_theme(asymptomatics_figure, total_theme), pt_per_unit = 1)
