# Packages
using JuMP
using LinearAlgebra
using Gurobi
using DataFrames
using Plots, Measures
using Colors

# Function to clean the terminal
function clearTerminal()

    # For Windows terminal
    if Sys.iswindows()
        Base.run(`cmd /c cls`)

    # Other terminals based on Unix
    else
        Base.run(`clear`)
    end
    
end

# https://www.esios.ree.es/es/unidades-fisicas
########## Data ##########
gTypes      = [ "Nuclear", "Carbón_1", "Carbón_2", "CC_1", "CC_2", "CC_3", "CC_4", "COGEN_1", "COGEN_2", "COGEN_3", "COGEN_4", "COGEN_5", "Fuel", "Derivados"]
N           = [         7,          2,          1,      5,     33,     18,     20,        97,        58,        86,        31,        12,      2,          52] #        - Number of units
C_m         = [        10,         66,         70,     80,     70,     68,     65,        55,        53,        50,        47,        45,    105,          98] # €/MWh  - Marginal cost
C_nl        = [        22,         18,         16,     13,     14,     16,     19,         6,         7,         8,         9,        10,      9,          12] # €/h    - No-load cost
Pg_Gen_lb   = [       700,        154,        256,   25.2,  120.2,  174.5,    187,      0.51,      1.61,      4.13,     10.32,        42,    1.8,         6.3] # MW     - Power lower bound
Pg_Gen_ub   = [    1016.7,     343.95,        570,   55.9,  267.1,  387.7,  415.3,      1.14,      3.57,      9.18,     22.94,      91.1,   3.95,       13.96] # MW     - Power upper bound
Pg_rr       = [        10,         75,        110,   22.5,   42.5,     60,     85,         2,         4,        10,        25,        50,      4,          15] # MW/h   - Maximum ramp rate for each generator of the cluster
Hg          = [         9,          6,          6,      3,      4,      4,      5,         2,         3,         3,         4,         5,      5,           5] # s      - Inertia constant
T_st        = [        48,         12,         15,      2,      2,      2,      2,         0,         0,         0,         1,         2,      1,           2] # h      - Startup time
T_mut       = [        48,          6,          6,      4,      4,      4,      4,         2,         3,         3,         3,         4,      1,           2] # h      - Minimum up time
T_mdt       = [        48,          6,          6,      4,      4,      4,      4,         2,         2,         2,         2,         2,      1,           1] # h      - Minimum down time

Rg_max      = Pg_Gen_ub.*0.05   # MW    - FR provision
Tg          = 8                 # s     - FR delivery time

# https://demanda.ree.es/visiona/peninsula/demandaau/acumulada/2025-01-20
# Hora:          00:00    01:00    02:00    03:00    04:00    05:00    06:00    07:00    08:00    09:00    10:00    11:00    12:00    13:00    14:00    15:00    16:00    17:00    18:00    19:00    20:00    21:00    22:00    23:00
# Total of power demand
Pd          = [ 15000,  15000,  15000,  15000,  15000,  15000,  15000,  15000,  15000,  15000,  15000,  15000,  15000,  15000,  15000,  15000,  15000,  15000,  15000,  15000,  15000,  15000,  15000,  15000] # MW 
# Power exportation
P_extern    = [     0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0,      0]
# Capacity factor of RES
cf_RES      = [     0,      1,      2,      3,      4,      5,      6,      7,      8,      9,     10,     11,     12,     13,     14,     15,     16,     17,     18,     19,     20,     21,     22,     23] ./ 200
# Power from Renewable Energy Source (RES)
P_RES       = 86 *10^3 # MW (Max power installed)

Rs_ub       = 0 # MW (EFR)
Ts          = 1 # s

# Max Rate of Change of Frecuency (RoCoF)
RoCoF_max = 1 # Hz/s

# Number of periods
T = length(Pd)

# Value of Lost Load (VoLL)
VoLL = 250000 # €/MWh

# Frecuency data
f_0         = 50 # Hz
Δf_max      = 0.8 # Hz
Δfss_max    = 0.2 # Hz

# Maximum power infeed
P_L_max = maximum(Pg_Gen_ub)
H_l = Hg[argmax(Pg_Gen_ub)]

# Load damping
D = 1.5 *10^-2 # %/Hz


# Check if the length of all data are the same
nTypes = length(gTypes)

fShowPlot = true
CaseStudy = 4

########## Model creation ##########
model = Model(Gurobi.Optimizer)
# set_optimizer_attribute(model, "OutputFlag", 0) # Se silencian las salidas pro defecto del solver Gurobi
set_optimizer_attribute(model, "Threads", min(Sys.CPU_THREADS, 12))
set_optimizer_attribute(model, "MIPFocus", 1)       # Prioriza encontrar factibles buenas
set_optimizer_attribute(model, "Heuristics", 0.2)
set_optimizer_attribute(model, "Presolve", 2)       # Agresivo
set_optimizer_attribute(model, "NumericFocus", 1)   # Robustez numérico
set_optimizer_attribute(model, "NonConvex", 2)


########## Variables ##########
@variable(model, Ng[1:nTypes, 1:T] >= 0, Int) # Number of generators of each nTypes
@variable(model, Ng_sd[1:nTypes, 1:T] >= 0, Int) # Number of generators of each nTypes - Shutdown
@variable(model, Ng_sg[1:nTypes, 1:T] >= 0, Int) # Number of generators of each nTypes - Start generating
@variable(model, Ng_st[1:nTypes, 1:T] >= 0, Int) # Number of generators set to Startup
@variable(model, Pg[1:nTypes, 1:T] >= 0) # Power from each cluster of generators
@variable(model, Rg[1:nTypes, 1:T] >= 0) # PRF from each generators
@variable(model, P_curt[1:T] >= 0) # RES Power curtailment
@variable(model, Rs[1:T] >= 0) # EFR from BESS (Battery Energy Storage Systems)
@variable(model, P_L[1:T] >= 0) # Power infeed
@variable(model, Dshed[1:T] >= 0) # Demand shedding


########## Objective ##########
@objective(model, Min, sum(Ng[g, t] * C_nl[g] + Pg[g, t] * C_m[g] for g in 1:nTypes, t in 1:T) + sum(VoLL * Dshed[t] for t in 1:T))

# System inertia
@expression(model, H[t=1:T], sum(Hg[g] * Pg_Gen_ub[g] * Ng[g, t] for g in 1:nTypes) - P_L[t] * H_l)

########## Constraints ##########
@constraint(model, [t in 1:T], sum(Pg[g, t] for g in 1:nTypes) + P_RES * cf_RES[t] - P_curt[t] + P_extern[t] == Pd[t] - Dshed[t])

# Number of active generators of each group constraint
@constraint(model, [g in 1:nTypes, t in 1:T], Ng[g, t] <= N[g])
@constraint(model, [t in 1:T], Ng[1, t] == N[1])

# Power from each generator constraints
@constraint(model, [g in 1:nTypes, t in 1:T], Pg[g, t] >= Pg_Gen_lb[g] * Ng[g, t])
@constraint(model, [g in 1:nTypes, t in 1:T], Pg[g, t] <= Pg_Gen_ub[g] * Ng[g, t])

# Largest power infeed constraints
@constraint(model, [g in 1:nTypes, t in 1:T], P_L[t] * Ng[g, t] >= Pg[g, t])

# PFR provision from g constraints
@constraint(model, [g in 1:nTypes, t in 1:T], Rg[g, t] <= Rg_max[g] * Ng[g, t])
@constraint(model, [g in 1:nTypes, t in 1:T], Rg[g, t] <= Pg_Gen_ub[g] * Ng[g, t] - Pg[g, t])

# EFR provision from RES constraints
@constraint(model, [t in 1:T], P_curt[t] <= P_RES * cf_RES[t])
@constraint(model, [t in 1:T], Rs[t] >= 0) # Rs (storage)
@constraint(model, [t in 1:T], Rs[t] <= Rs_ub)

if (CaseStudy == 3) || (CaseStudy == 4)
    # RoCoF constraint
    @constraint(model, [t in 1:T], P_L[t] * f_0 <= RoCoF_max * 2 * H[t])

    # Quasi-steady-state security constraint
    @constraint(model, [t in 1:T], Rs[t] + sum(Rg[g, t] for g in 1:nTypes) - P_L[t] <= Δfss_max * D * Pd[t])
    @constraint(model, [t in 1:T], -(Rs[t] + sum(Rg[g, t] for g in 1:nTypes) - P_L[t]) <= Δfss_max * D * Pd[t])


    # RotatedSecondOrderCone()
    # 2x₁x₂ ≥ ∣∣x₃₋ₙ∣∣²
    # @constraint(model, [x₁, x₂, x₃, ..., xₙ] in RotatedSecondOrderCone())

    # Auxiliar variables
    @variable(model, x[1:5, 1:T])
    @constraint(model, [t in 1:T], x[1, t] == (H[t]/f_0) - (Rs[t] * Ts) / (4*Δf_max))
    @constraint(model, [t in 1:T], x[2, t] == sum(Rg[g, t] for g in 1:nTypes) / Tg)
    @constraint(model, [t in 1:T], x[3, t] == (P_L[t] - Rs[t]) / sqrt(4*Δf_max))
    @constraint(model, [t in 1:T], x[4, t] == D * Pd[t])
    @constraint(model, [t in 1:T], x[5, t] == (P_L[t] - Rs[t]) / 4)
end

if CaseStudy == 3
    # Nadir constraint without load damping
    # @constraint(model, [t in 1:T], (H[t] / f_0 - Rs[t] * Ts / (4*Δf_max)) * sum(Rg[g, t] for g in 1:nTypes) >= (P_L[t] - Rs[t])^2 * Tg / (4*Δf_max))
    
    # Nadir constraint without load damping using RotatedSecondOrderCone()
    # x1*x2 >= x3^2
    @constraint(model, [t in 1:T], [x[1,t]/2, x[2,t], x[3,t]] in RotatedSecondOrderCone())
end


if CaseStudy == 4
    # Nadir constraint with load damping
    # @constraint(model, [t in 1:T], (H[t] / f_0 - Rs[t] * Ts / (4*Δf_max)) * sum(Rg[g, t] for g in 1:nTypes) >= (P_L[t] - Rs[t])^2 * Tg / (4*Δf_max) - (P_L[t] - Rs[t]) * Tg * D * Pd[t] / 4)
    # Nadir constraint with load damping using RotatedSecondOrderCone()
    @variable(model, s_aux[1:T] >= 0)
    @variable(model, t_aux[1:T] >= 0)
    # x1*x2 >= x3^2 - x5*x4
    # x1*x2 = s^2
    @constraint(model, [t in 1:T], [x[1, t]/2, x[2, t], s_aux[t]] in RotatedSecondOrderCone())
    # x5*x4 = t^2
    @constraint(model, [t in 1:T], [x[5, t]/2, x[4, t], t_aux[t]] in RotatedSecondOrderCone())
    # s^2 >= x3^2 - t^2 -----> s^2 + t^2 >= x3^2
    @variable(model, y_aux[1:T] >= 0)
    @variable(model, z_aux[1:T] >= 0)
    # y[t] ≥ s_aux[t]^2  ⇔  (y[t], 1/2, s_aux[t])
    @constraint(model, [t in 1:T], [y_aux[t], 0.5, s_aux[t]] in RotatedSecondOrderCone())
    # z[t] ≥ t_aux[t]^2  ⇔  (z[t], 1/2, t_aux[t])
    @constraint(model, [t in 1:T], [z_aux[t], 0.5, t_aux[t]] in RotatedSecondOrderCone())
    # y[t] + z[t] ≥ x[3,t]^2  ⇔  ((y[t] + z[t]), 1/2, x[3,t])
    @constraint(model, [t in 1:T], [(y_aux[t] + z_aux[t]), 0.5, x[3,t]] in RotatedSecondOrderCone())
end

clearTerminal()
########## Optimization ##########
optimize!(model)



########## Solution display ##########
if termination_status(model) == OPTIMAL || termination_status(model) == LOCALLY_SOLVED
    
    println("\n\n##### $(termination_status(model)) solution found #####")

    Rg_values = [round(sum(value(Rg[g, t]) for g in 1:nTypes), digits=2) for t in 1:T]

    println("Rg values = ", Rg_values, " GW")

    ########## Solution plotting ##########
    if fShowPlot
        plotColors = Dict(
            "Nuclear"   => colorant"#454394",  # purple

            # Combined Cycle (orange)
            "CC_1"      => colorant"#956312",
            "CC_2"      => colorant"#a68702",
            "CC_3"      => colorant"#b89502",
            "CC_4"      => colorant"#c9a302",

            # Carbón (marrones)
            "Carbón_1"  => colorant"#934e2c",
            "Carbón_2"  => colorant"#c56a3e",

            # Cogeneración y residuos (pink)
            "COGEN_1"   => colorant"#bc7db7",
            "COGEN_2"   => colorant"#c690c1",
            "COGEN_3"   => colorant"#cfa2cb",
            "COGEN_4"   => colorant"#d8b4d5",
            "COGEN_5"   => colorant"#e2c7df",

            # Otros
            "Fuel"      => colorant"#003366",  # Turbina de gas (darkblue)
            "Derivados" => colorant"#a08000"   # Turbina de vapor (ocre)
        )

        x = 0:T-1
        Pg_values = [value(Pg[g, t]/1000) for g in 1:nTypes, t in 1:T]
        y = Array{Float64}(undef, nTypes, T)
        y[1, :] .= Pg_values[1, :]
        P_curt_values = [value(P_curt[t])/1000 for t in 1:T]
        plot(x, y[1, :],
            size=(1200, 700),
            dpi=160,
            legend=:outerright,
            left_margin=10mm,
            bottom_margin=10mm,
            xticks=(x, round.(Int, cf_RES*P_RES)),
            xrotation=45,
            yticks=[(5*floor((minimum(-P_curt_values))/5)):5:maximum(Pd/1000); 20],
            ylims=(minimum(-P_curt_values), maximum((Pd-P_extern)/1000) + 0.2),
            fillrange=0,
            lw=0.5,
            linecolor=plotColors[gTypes[1]],
            label=gTypes[1],
            color=plotColors[gTypes[1]],
            xlabel="Potencia de fuentes renovables (MW)",
            ylabel="Demanda (GW)")

        for g in 2:nTypes
            y[g, :] .= y[g-1, :] .+ Pg_values[g, :]
            plot!(x, y[g, :], fillrange=y[g-1, :], lw=0.5, linecolor=plotColors[gTypes[g]], label=gTypes[g], color=plotColors[gTypes[g]])
        end
        y_RES = [value(P_RES * cf_RES[t] - P_curt[t])/1000 for t in 1:T]
        plot!(x, y_RES + y[nTypes, :], fillrange=y[nTypes, :], lw=0, label="RES power", color="green")
        plot!(x, Pd/1000, label="Demand", color="red", lw=3)
        plot!(x, P_extern/1000, label="Import&Export", fillrange=0)
        plot!(x, -P_curt_values, label="P_curt", fillrange=0)

        display(current())
    end
    println(solve_time(model))

else
    println("ERROR: ", termination_status(model))

end