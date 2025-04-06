sources = string(pwd(), "\\heatgrid\\src\\")
push!(LOAD_PATH, sources)

using Revise
using Pkg
if isfile("Project.toml") && isfile("Manifest.toml")
    Pkg.activate(".")
end

# using heatGrid, Tools, Interpolations, XLSX, MAT, Random, Dates, Plots, SciMLBase, OrdinaryDiffEq, DifferentialEquations, CSV
using heatGrid, Tools, Interpolations

using XLSX, MAT, Random, Dates, Plots, JSON, CSV

using DelimitedFiles, Statistics, DataFrames

using SparseArrays: Vector
using NLsolve
using SparseArrays
using JuMP, GLPK, Gurobi

Random.seed!(1111)

# Creating the MPC problem and the topology of the Verbier DHN
topology_file = "heatgrid\\satom_case_study\\topology_satom_extension_1.xlsx"
details_file = "heatgrid\\satom_case_study\\details_satom_extension_1.json"
measurements_file = "heatgrid\\satom_case_study\\case_study_satom_extension_1_measurements.xlsx"
hydraulic_file = "heatgrid\\satom_case_study\\hydraulic_results_extension_1_v5.mat"
thermal_file = "heatgrid\\satom_case_study\\thermal_results_extension_1_v5.mat"
pydhn_mw_file = "heatgrid\\satom_case_study\\pipes_mw_from_pydhn.mat"

function get_init_temp_file(name)
    base = "heatgrid\\satom_case_study\\"
    return string(base, "init_$(name).xlsx")
end
    


# Some constant values
rho_pipe = 940.0
cp_pipe = 2000.0

# Simulation control values
dt_dynamic = 1.0/60.0 # 1min
dt = 1.0 # 15 min
n_steps = Int(dt/dt_dynamic)
windowSpan = 3.0
alpha_losses = 0.2
t0 = 1.0
# tf = 2.0 * 7.0 # only the present data
tf = 2.0 * 30.0 * 24.0
# n_dynamic should be 18691 = (1247 - 1) * 15 + 1
sliding_window_size = 18
parameters = Dict(
    "dt" => dt,
    "dt_dynamic" => dt_dynamic,
    "windowSpan" => windowSpan,
    "alpha_losses" => alpha_losses,
    "t0" => t0,
    "tf" => tf,
    "sliding_window_size" => sliding_window_size,
    "sources_indices" => [],
    "source_slack_index" => 1,
    "init_temp" => 50.0,
)

mpc_dx = Tools.create_MPC_Problem(topology_file, 
                                        parameters["t0"], 
                                        parameters["tf"],
                                        parameters["dt"], 
                                        parameters["dt_dynamic"], 
                                        parameters["windowSpan"], 
                                        [], 
                                        0.0, # dx
                                        20, # nvols (has to be 0 if dx != 0)
                                        parameters["init_temp"],
                                        true,
                                        false)

details = JSON.parsefile(details_file)
sources = details["sources_nodes"]
for node_id in sources
    heatGrid.addSource(mpc_dx, node_id, 1000.0*10^3, 0.1)
    # mpc_dx.network.slackIndex = node_id
end

# MAIN is the slack wich has the node_id = 193
main_node_rel_id = details["MAIN_id"]
mpc_dx.network.slackIndex = main_node_rel_id

# ts_measures = DataFrame(XLSX.readtable(measurements_file, "ts"))
# tr_measures = DataFrame(XLSX.readtable(measurements_file, "tr"))
mc_measures = DataFrame(XLSX.readtable(measurements_file, "mc"))
loads_measures = DataFrame(XLSX.readtable(measurements_file, "loads"))
pipes_t_exterior = DataFrame(XLSX.readtable(topology_file, "pipes"))[!, "t_ext"] # All shoud be in 8 °C

init_pipes = DataFrame(XLSX.readtable(get_init_temp_file("pipes"), "Sheet1"))
init_nodes = DataFrame(XLSX.readtable(get_init_temp_file("nodes"), "Sheet1"))

tsin_pipes_init = init_pipes[!, "tsin"]
tsout_pipes_init = init_pipes[!, "tsout"]
ts_pipes_init = init_pipes[!, "ts"]
trin_pipes_init = init_pipes[!, "trin"]
trout_pipes_init = init_pipes[!, "trout"]
tr_pipes_init = init_pipes[!, "tr"]

fd_pipes_found = init_pipes[!, "fd_s"]

ts_nodes_init = init_nodes[!, "ts"]
tr_nodes_init = init_nodes[!, "tr"]

for node_id in 1:mpc_dx.network.Nnodes
    if string(node_id) in names(loads_measures)
        load = loads_measures[!, string(node_id)]
        n_c = min(size(mpc_dx.load)[2], length(load))
        mpc_dx.load[node_id, 1:n_c] = load[1:n_c]
    end
end

# Load are in W, reconvert them back to kW 
mpc_dx.load .*= 0.001
mpc_dx.load[sources, :] .= 0.0

mpc_dx.mc .= 0.0
# sets the Tss of the sources
mpc_dx.Tss .= 95.0 # TSS imposed

# Here no ms is imposed (all are compensated by the source node)
# Set the mc of the consumers imposed by measures
for node_id in 1:mpc_dx.network.Nnodes
    if string(node_id) in names(mc_measures)
        mc = mc_measures[!, string(node_id)]
        n_c = min(size(mpc_dx.mc)[2], length(mc))
        if !(node_id in sources) 
            mpc_dx.mc[node_id, 1:n_c] = mc[1:n_c]
        end
    end
end

# Sets the external temperatures, used only for the thermal model but do it here 
mpc_dx.network.pipes.Tex .= pipes_t_exterior
# mpc_dx.mc[sources, :] .= 0.0 # verify this
dynamic_pressure_source = 1000000.0 # in Pa 
dynamic_head_source = dynamic_pressure_source / heatGrid.RHO / 9.81

function fmassSatomV2!(F,x,netM,mbase)
    
    netM.nodes.H[netM.producers.indexes] .= dynamic_head_source
    Nnodes = netM.Nnodes
    Npipes = netM.Npipes
    # println(mbase)
    netM.pipes.massflows = x[1:Npipes] .* mbase
    for j=1:Nnodes
        if j != netM.slackIndex
            netM.nodes.H[j] = x[Npipes+j] .* dynamic_head_source
        end
    end

    heatGrid.update_pipes_k_pydhn(netM) #!! Change here!
    heatGrid.updateSlackMassFlow!(netM)

    F[Npipes+1:Npipes+Nnodes] = heatGrid.massBalance(netM) ./ mbase
    F[1:Npipes] = heatGrid.dH(netM) ./ dynamic_head_source
end

function fmassSatom(F,x,netM)
    Nnodes = netM.Nnodes
    Npipes = netM.Npipes

    if (netM.slackIndex > heatGrid.ONE_INDEX)
        netM.nodes.H[1:netM.slackIndex-1] = x[1:netM.slackIndex-1]*netM.Hbase
    end
    netM.nodes.H[netM.slackIndex+1: Nnodes] = x[netM.slackIndex:Nnodes-1]*netM.Hbase
    netM.pipes.massflows = x[Nnodes:Nnodes+Npipes-1]*netM.mbase

    heatGrid.update_pipes_k_pydhn(netM) #!! Change here!

    heatGrid.updateSlackMassFlow!(netM)
    F[1:Nnodes] = heatGrid.massBalance(netM) ./ netM.mbase
    F[netM.slackIndex:Nnodes-1] = F[netM.slackIndex+1:Nnodes]
    F[Nnodes:Nnodes+Npipes-1] = heatGrid.dH(netM) ./netM.Hbase
end

function solveMassSatom(mpc, i)

    network = mpc.network
    Npipes = network.Npipes
    Nnodes = network.Nnodes
    x0 = zeros(1:Npipes+Nnodes)
    previous_mdot = zeros(size(mpc.network.pipes.massflows))
    if i == 1
        # previous_mdot = (mpc.network.pipes.D.^2 .* 0.25 .* heatGrid.RHO .* 1.0)
        previous_mdot .= 100.0
    else
        previous_mdot = mpc.mw[:,i-1] .* 1.0 # with 20% margin
        # x0[1:Npipes] = mpc.mw[:,i-1] ./ previous_mdot
        # dynamic_pressure_source = 1300.0 * 1000.0 # in Pa 
        # dynamic_head_source = dynamic_pressure_source / heatGrid.RHO / 9.81
        # x0[Npipes+1:Npipes+Nnodes] = 1 ./ dynamic_head_source
    end

    mpc.network.pipes.massflows .= 1.0 .* previous_mdot
    mpc.network.nodes.H .= 100.0
    # mpc.network.nodes.H .= 0.0

    Fmass!(F,x) = fmassSatomV2!(F,x,network,maximum(abs.(previous_mdot)))
    sol = NLsolve.nlsolve(Fmass!, x0,method = :trust_region, show_trace=true, iterations = 2000, xtol=1e-8, ftol=1e-10)
    return sol
end

function get_mass_balance(mpc, i)
    mpc.network.nodes.mc .= mpc.mc[:,i]
    mpc.network.nodes.ms[mpc.network.producers.indexes] .= mpc.ms[:,i]
    mpc.network.pipes.massflows .= mpc.mw[:,i]
    mass_balance = heatGrid.massBalance(mpc.network)
    return mass_balance
end

function perform_hydraulic_solve(mpc, i)
    mpc.network.nodes.mc .= mpc.mc[:,i]
    mpc.network.nodes.ms[mpc.network.producers.indexes] .= mpc.ms[:,i]
    mpc.network.mbase = sum(mpc.mc[:,i]) # utiliser peut etre la somme sur les mc
    # mpc.network.mbase = 1.0 # utiliser peut etre la somme sur les mc
    mpc.network.Hbase = dynamic_head_source
    sol = solveMassSatom(mpc, i)
    return sol
end

# Initialize temperatures
mpc_dx.Ts .= ts_nodes_init
mpc_dx.network.nodes.Ts .= ts_nodes_init
mpc_dx.Tr .= tr_nodes_init
mpc_dx.network.nodes.Tr .= tr_nodes_init
mpc_dx.TsDynamic .= ts_nodes_init
mpc_dx.TrDynamic .= tr_nodes_init

mpc_dx.Tsin_Dynamic .= tsin_pipes_init
mpc_dx.Tsout_Dynamic .= tsout_pipes_init
mpc_dx.Trin_Dynamic .= trin_pipes_init
mpc_dx.Trout_Dynamic .= trout_pipes_init
mpc_dx.network.pipes.Tsin .= tsin_pipes_init
mpc_dx.network.pipes.Tsout .= tsout_pipes_init
mpc_dx.network.pipes.Trin .= trin_pipes_init
mpc_dx.network.pipes.Trout .= trout_pipes_init

for (i, p) in enumerate(mpc_dx.pipesS)
    p.T .= tsout_pipes_init[i]
    # p.T[end] .= tsout_pipes_init[i]
end
for (i, p) in enumerate(mpc_dx.pipesS_stored)
    p.T .= tsout_pipes_init[i]
end
for (i, p) in enumerate(mpc_dx.pipesR)
    p.T .= trout_pipes_init[i]
end
for (i, p) in enumerate(mpc_dx.pipesR_stored)
    p.T .= trout_pipes_init[i]
end

heatGrid.storeDescritizedPipes(mpc_dx)



# nt = 6 * 24 * 4
nt = size(mpc_dx.load, 2)
dH_found = zeros(size(mpc_dx.mw))
dP_found = zeros(size(mpc_dx.mw))
h_found = zeros(size(mpc_dx.Ts))
k_pipes_found = zeros(size(mpc_dx.mw))

n_save = nt
# n_save = 5
PERFORM_HYDRAULIC = true
if PERFORM_HYDRAULIC
    # nt_total = size(mpc_dx.load, 2)
    for i = 13:n_save
        println("Hydraulic step $(i)/$(nt)")
        sol = perform_hydraulic_solve(mpc_dx, i)
        mpc_dx.mc[:,i] .= mpc_dx.network.nodes.mc
        mpc_dx.ms[:,i] .= mpc_dx.network.nodes.ms[mpc_dx.network.producers.indexes]
        mpc_dx.mw[:,i] .= mpc_dx.network.pipes.massflows
        k_pipes_found[:,i] .= mpc_dx.network.pipes.k
        dh_found[:,i] = mpc_dx.network.incidenceT * mpc_dx.network.nodes.H 
        h_found[:,i] = mpc_dx.network.nodes.H
        if !converged(sol)
            println("--> step $(i) not converged! review")
            push!(not_converged_steps, i)
        end
    end

    # n_save = nt

    # Saving hydraulic results
    file = matopen(hydraulic_file, "w")
    write(file, "tss", mpc_dx.Tss[:,1:n_save])
    write(file, "mc", mpc_dx.mc[:,1:n_save])
    write(file, "ms", mpc_dx.ms[:,1:n_save])
    write(file, "pipes_mw", mpc_dx.mw[:,1:n_save])
    write(file, "pipes_dh", dh_found[:,1:n_save])
    write(file, "pipes_k", k_pipes_found[:,1:n_save])
    write(file, "nodes_h", h_found[:,1:n_save])
    close(file)

else 
    file = matopen(hydraulic_file, "r")
    mpc_dx.Tss[:, 1:n_save] .= read(file, "tss")[:,1:n_save]
    mpc_dx.ms[:, 1:n_save] .= read(file, "ms")[:,1:n_save]
    mpc_dx.mc[:, 1:n_save] .= read(file, "mc")[:,1:n_save]
    mpc_dx.mw[:, 1:n_save] .= read(file, "pipes_mw")[:,1:n_save]
    dh_found[:, 1:n_save] .= read(file, "pipes_dh")[:,1:n_save]
    k_pipes_found[:, 1:n_save] .= read(file, "pipes_k")[:,1:n_save]
    close(file)
end


# Thermal model
function custom_integration_dt(MPC, p, Tsin, Tex, i)
    dT = zeros(p.N)
    diameter = MPC.network.pipes.D[i]
    outer_diameter = diameter + 2*MPC.network.pipes.TT[i]
    # length = diameter + MPC.network.pipes.L[i]
    area_pipe = pi * (outer_diameter/2)^2
    area_water = pi * (diameter/2)^2 
    rho_water = heatGrid.get_rho_pydhn(Tsin) # can be ts can be tr
    cp_water = heatGrid.get_cp_pydhn(Tsin)

    # println(rho_water)
    # println(cp_water)

    pipe_mass_per_unit_length = (area_pipe - area_water) * rho_pipe # Kg/m
    water_mass_per_unit_length = area_water * rho_water # Kg/m
    cp_effective = cp_water + (pipe_mass_per_unit_length/water_mass_per_unit_length) * cp_pipe
    coeff_area = area_water/area_pipe
    coeff_cp = cp_water/cp_effective
    loss_coeff = (MPC.network.pipes.H[i] * pi * diameter) / (rho_water * area_pipe * cp_effective)
    # println(loss_coeff)
    # println((rho_water * area_pipe * cp_effective))
    # println((area_pipe))
    # println((coeff_cp))
    # println((coeff_area))
    if p.v >= 0.0
        dT[1] = p.v*coeff_area*coeff_cp*(-p.T[1]+Tsin)/p.dx - loss_coeff * (p.T[1]-Tex)
        for j = 2:p.N
            dT[j] = p.v*coeff_area*coeff_cp*(-p.T[j]+p.T[j-1])/p.dx - loss_coeff*(p.T[j]-Tex)
        end
    else
        dT[p.N] = p.v*coeff_area*coeff_cp*(p.T[p.N]-Tsin)/p.dx - loss_coeff*(p.T[p.N]-Tex)
        for j = 1:p.N-1
            dT[p.N-j] = p.v*coeff_area*coeff_cp*(p.T[p.N-j]-p.T[p.N-j+1])/p.dx - loss_coeff*(p.T[p.N-j]-Tex)
        end
    end
    return dT
end

# Maybe review CP and RHO to make it identical to pydhn rho, cp
# rho = 986 # (from PyDHN)
# cp = 4190 # (from PyDHN)
# Verify the load (there are already in W)
#!!!! mpc_dx.network.pipes.H .*= pi .* 0.263

mpc_dx.network.pipes.Tex .= 8.0
heatGrid.updateDescritizedPipes(mpc_dx)
tic = Dates.now()
heatGrid.dynamicSimulationWithCustomIntegrationSATOM(mpc_dx, 1, n_save, custom_integration_dt)
toc = Dates.now()
# Saving hydraulic results
file = matopen(thermal_file, "w")
dynamic_step = n_save*n_steps
write(file, "tsDynamic", mpc_dx.TsDynamic[:,1:dynamic_step])
write(file, "ts", mpc_dx.Ts[:,1:nt])
write(file, "trDynamic", mpc_dx.TrDynamic[:,1:dynamic_step])
write(file, "tr", mpc_dx.Tr[:,1:nt])
write(file, "tsin", mpc_dx.Tsin_Dynamic[:, 1:dynamic_step])
write(file, "trin", mpc_dx.Trin_Dynamic[:, 1:dynamic_step])
write(file, "tsout", mpc_dx.Tsout_Dynamic[:, 1:dynamic_step])
write(file, "trout", mpc_dx.Trout_Dynamic[:, 1:dynamic_step])
write(file, "tss", mpc_dx.Tss[:,1:nt])
write(file, "Ps", mpc_dx.Ps[:,1:nt])
write(file, "real_Pc", mpc_dx.Cons_Dynamic[:, 1:dynamic_step])
write(file, "load", mpc_dx.Load_Dynamic[:, 1:dynamic_step])
write(file, "pipes_mw_dyn", mpc_dx.Mw_Dynamic[:, 1:dynamic_step])
close(file)

# Total simulation file used for clustering/ML
total_simulation_file = "heatgrid\\satom_case_study\\satom_simulation_results.mat"
file = matopen(total_simulation_file, "w")
dynamic_step = n_save*n_steps
write(file, "tsDynamic", mpc_dx.TsDynamic[:,1:dynamic_step])
write(file, "ts", mpc_dx.Ts[:,1:nt])
write(file, "trDynamic", mpc_dx.TrDynamic[:,1:dynamic_step])
write(file, "tr", mpc_dx.Tr[:,1:nt])
write(file, "tsin", mpc_dx.Tsin_Dynamic[:, 1:dynamic_step])
write(file, "trin", mpc_dx.Trin_Dynamic[:, 1:dynamic_step])
write(file, "tsout", mpc_dx.Tsout_Dynamic[:, 1:dynamic_step])
write(file, "trout", mpc_dx.Trout_Dynamic[:, 1:dynamic_step])
write(file, "tss", mpc_dx.Tss[:,1:nt])
write(file, "Ps", mpc_dx.Ps[:,1:nt])
write(file, "mc", mpc_dx.mc[:,1:nt])
write(file, "real_Pc", mpc_dx.Cons_Dynamic[:, 1:dynamic_step])
write(file, "load", mpc_dx.Load_Dynamic[:, 1:dynamic_step])
write(file, "pipes_mw_dyn", mpc_dx.Mw_Dynamic[:, 1:dynamic_step])
close(file)



# Compute and verify mass balances for all nodes and all steps
mw_heatgrid = deepcopy(mpc_dx.mw)
file = matopen(pydhn_mw_file, "r")
pydhn_mw = read(file, "pipes_mw")
close(file) 
# mpc_dx.mw[:,1:size(pydhn_mw)[2]] .= pydhn_mw
mpc_dx.mw[:,1:size(pydhn_mw)[2]] .= mw_heatgrid[:,1:size(pydhn_mw)[2]]
is_mass_balance_violated = []
for t=1:n_save
    mass_balance_nodes = get_mass_balance(mpc_dx, t) 
    for n=1:mpc_dx.network.Nnodes
        mass_b = mass_balance_nodes[n]
        if abs(mass_b) >= 0.01
            println("Mass balance violated at node $(n) for step $(t)")
            if !(n in is_mass_balance_violated)
                push!(is_mass_balance_violated, n)
            end
        end
    end
end



# fd = k_pipes_found .* (4 .* pi^2 .* heatGrid.RHO .* (0.263/2)^5) ./ mpc_dx.network.pipes.L