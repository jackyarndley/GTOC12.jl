using Printf
using FiniteDiff
using SparseArrays
using Statistics
using JuMP
using MosekTools
using Gurobi
using DiffEqCallbacks
using DifferentialEquations
using ForwardDiff
using CairoMakie
using GLMakie
using GraphMakie
using Graphs
using ColorSchemes
using SwarmMakie

CairoMakie.activate!() 

include("assorted.jl")
include("configuration.jl")
include("integration.jl")
include("lambert.jl")
include("scp.jl")
include("mip.jl")
include("transforms.jl")
include("plotting.jl")
include("data_handling.jl")

# Scaling constants
μ = 1.32712440018e11 # km^3 / s^2
r_scale = 1.49597870691e8 # km / LU
v_scale = sqrt(μ/r_scale) # km/s / LU/TU
t_scale = r_scale/v_scale # s / TU
a_scale = v_scale/t_scale
m_scale = 3000 # kg / MU
f_scale = m_scale*a_scale # N

# Derived scaling constants
day_scale = 24*60*60/t_scale

# Global variable for number of SCP iterations to take
scp_iterations = 10

# Acceleration limit for low-thrust system assuming full mass
maximum_mass = 3000.0 / m_scale
dry_mass = 500.0 / m_scale
miner_mass = 40.0 / m_scale
mining_rate = 10.0 / m_scale / 365.25 / day_scale

# Scaling nondimensional
μ = 1.0
μ_planets = [3.24858592000e5, 3.98600435436e5, 4.28283752140e4] .* t_scale^2 / r_scale^3
rp_min_planets_actual = ([6351.0, 6678.0, 3689.0]) ./ r_scale 
rp_min_planets = ([6351.0, 6678.0, 3689.0] .+ 1.0) ./ r_scale # Add 1km extra
g0 = 9.80665e-3/(r_scale/t_scale^2) # LU/TU^2
isp = 4000/t_scale # TU
g0_isp = g0*isp

thrust_scale = 1.0
thrust = thrust_scale*0.6e-3/(m_scale*r_scale/t_scale^2)

# Maximum time
mjd_start = 64328
mjd_end = 69807
maximum_time = (mjd_end - mjd_start)*day_scale

# For the spacing will space nodes out until the last node which will be shorter if needed
node_time_spacing = 5.0 * day_scale

# Asteroids
asteroids_classical = read_file_classical_data("data/GTOC12_Asteroids_Data.txt")
asteroids_cartesian = classical_to_cartesian(asteroids_classical)
number_asteroids = size(asteroids_classical, 2)

# Planets
planets_classical = read_file_classical_data("data/GTOC12_Planets_Data.txt")
planets_cartesian = classical_to_cartesian(planets_classical)

bad_asteroids = (asteroids_classical[2, :] .>= 0.4) .|  (asteroids_classical[1, :] .>= 4.0) .| (asteroids_classical[3, :] .>= deg2rad(15))

all_ids = collect(1:60000)

pruned_ids = all_ids[.!bad_asteroids];