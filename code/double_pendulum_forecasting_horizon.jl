using Distributions
using Distances
using DataFrames,  DataFramesMeta
using DynamicalSystems, DifferentialEquations
using DynamicalSystems.Systems: double_pendulum
using Plots
using CSV
using ProgressMeter
"""
Types for representing different types of double pendulums


"""
abstract type DoublePendulum end
struct SameLengthSameMass <: DoublePendulum
    mass::Float64
    length::Float64
end
SameLengthSameMass() = SameLengthSameMass(1.,1.)


struct SameLengthDifferentMass <: DoublePendulum
    massInner::Number
    massOuter::Number
    length::Number
    SameLengthDifferentMass() = SameLengthDifferentMass(0.5, 1, 1)
end

struct DifferentLengthSameMass <: DoublePendulum
    mass::Number
    lengthInner::Number
    lengthOuter::Number
    DifferentLengthSameMass() = DifferentLengthSameMass(1, 0.5, 1)
end

struct DifferentLengthDifferentMass <: DoublePendulum
    massInner::Number
    massOuter::Number
    lengthInner::Number
    lengthOuter::Number
    DifferentLengthDifferentMass() = DifferentLengthDifferentMass(0.5, 1, 0.5, 1)
end

"""
    Types for representing the position of a double pendulum in both polar and cartesian, with converters
"""
abstract type DoublePendulumState end
struct DoublePendulumStatePolar <: DoublePendulumState
    lengthInner::Number
    lengthOuter::Number
    θInner::Number
    θOuter::Number
end
DoublePendulumStatePolar(pend::SameLengthSameMass, θ1, θ2) = DoublePendulumStatePolar(pend.length, pend.length, θ1, θ2)
DoublePendulumStatePolar(pend::SameLengthDifferentMass, θ1, θ2) = DoublePendulumStatePolar(pend.length, pend.length, θ1, θ2)
DoublePendulumStatePolar(pend::DifferentLengthDifferentMass, θ1, θ2) = DoublePendulumStatePolar(pend.length1, pend.length2, θ1, θ2)
DoublePendulumStatePolar(pend::DifferentLengthSameMass, θ1, θ2) = DoublePendulumStatePolar(pend.length1, pend.length2, θ1, θ2)

struct DoublePendulumStateCartesian <: DoublePendulumState
    xInner::Number
    yInner::Number
    xOuter::Number
    yOuter::Number
end
DoublePendulumStateCartesian(polar_pt::DoublePendulumStatePolar) = begin
    xInner = polar_pt.lengthInner * cos(polar_pt.θInner)
    yInner = polar_pt.lengthInner * sin(polar_pt.θInner)
    xOuter = polar_pt.lengthOuter * cos(polar_pt.θOuter)
    yOuter = polar_pt.lengthOuter * sin(polar_pt.θOuter)
    return DoublePendulumStateCartesian(xInner,yInner,xOuter,yOuter)
end


"""
    Convert a dataframe of a trajectory in polar coordinates (θ-space) to a dataframe of the cartesian trajectory
"""
function convert_theta_space_to_cartesian_space(theta_trajectory::DataFrame, pend::DoublePendulum)
    cartesian_trajectory = DataFrame(time=zeros(nrow(theta_trajectory)), xInner=zeros(nrow(theta_trajectory)), yInner=zeros(nrow(theta_trajectory)), xOuter=zeros(nrow(theta_trajectory)), yOuter=zeros(nrow(theta_trajectory)))
    for (i,row) in enumerate(eachrow(theta_trajectory))
        cartesian_state = DoublePendulumStateCartesian(DoublePendulumStatePolar(pend, row.θInner, row.θOuter))
        cartesian_trajectory.time[i] = i
        cartesian_trajectory.xInner[i] = cartesian_state.xInner
        cartesian_trajectory.yInner[i] = cartesian_state.yInner
        cartesian_trajectory.xOuter[i] = cartesian_state.xInner + cartesian_state.xOuter
        cartesian_trajectory.yOuter[i] = cartesian_state.yInner + cartesian_state.yOuter
    end
    return cartesian_trajectory
end


"""
    Functions to make a double pendulum instance from a pandulum object   
"""

# the L1 and L2 keyword args might be flipped in the double_pendulum constructor 
make_initial_ODE(pend::SameLengthSameMass, initial_theta::Vector{Float64}) = double_pendulum([initial_theta[1], 0, initial_theta[2], 0], G=10.0, L1 = pend.length, L2 = pend.length, M1 = pend.mass, M2 = pend.mass)
make_initial_ODE(pend::SameLengthDifferentMass, initial_theta::Vector{Float64}) = double_pendulum([initial_theta[1], 0, initial_theta[2], 0], G=10.0, L1 = pend.length, L2 = pend.length, M1 = pend.massInner, M2 = pend.massOuter)
make_initial_ODE(pend::DifferentLengthSameMass, initial_theta::Vector{Float64}) = double_pendulum([initial_theta[1], 0, initial_theta[2], 0], G=10.0, L1 = pend.lengthInner, L2 = pend.lengthOuter, M1 = pend.mass, M2 = pend.mass)
make_initial_ODE(pend::DifferentLengthDifferentMass, initial_theta::Vector{Float64}) = double_pendulum([initial_theta[1], 0, initial_theta[2], 0], G=10.0, L1 = pend.lengthInner, L2 = pend.lengthOuter, M1 = pend.massInner, M2 = pend.massOuter)

"""
    runs a double pendulum with both a true initial condition and a 
    and returns true and predicted trajectories in dataframes 
"""
function run_double_pendulum(;pend::DoublePendulum=SameLengthSameMass(), measurement_error_std_dev::Number=0.001, output_num_timepoints::Int=5000)
    true_initial_theta = rand(Uniform(0.5*pi, 1.5*pi),2)
    trueODE = make_initial_ODE(pend, true_initial_theta)
    num_timepoints = output_num_timepoints / 100  # the dynamical system gives 100x the num_timepoints passed because delta_t = 0.01. 
    true_trajacetory_in_theta_space = DataFrame(Matrix(trajectory(trueODE, num_timepoints)), [:θInner, :dθInner, :θOuter, :dθOuter])
    true_trajacetory_in_cartesian = convert_theta_space_to_cartesian_space(true_trajacetory_in_theta_space, pend)
    names!(true_trajacetory_in_cartesian, [:time, :xInner_True, :yInner_True, :xOuter_True, :yOuter_True])


    measured_initial_theta = true_initial_theta + rand(Normal(0, measurement_error_std_dev), 2)

    measuredODE = make_initial_ODE(pend, measured_initial_theta)
    forecasted_trajectory_in_theta_space = DataFrame(Matrix(trajectory(measuredODE, num_timepoints)), [:θInner, :dθInner, :θOuter, :dθOuter])
    forecasted_trajectory_in_cartesian = convert_theta_space_to_cartesian_space(forecasted_trajectory_in_theta_space, pend)
    names!(forecasted_trajectory_in_cartesian, [:time, :xInner_Forecasted, :yInner_Forecasted, :xOuter_Forecasted, :yOuter_Forecasted])

    df = innerjoin(true_trajacetory_in_cartesian, forecasted_trajectory_in_cartesian, on = :time, makeunique=true)

    df = @transform(df, 
                    outerDifference = sqrt.((:xOuter_True .- :xOuter_Forecasted).^2 + (:yOuter_True .- :yOuter_Forecasted).^2),
                    innerDifference = sqrt.((:xInner_True .- :xInner_Forecasted).^2 + (:yInner_True .- :yInner_Forecasted).^2)
                    )

    return(df)
end

function time_until_forecast_fails(df; ϵ=0.01 )
    min_time_inner = 0
    min_time_outer = 0

    for row in eachrow(df)
        if row.innerDifference > ϵ && min_time_inner == 0
            min_time_inner = row.time
        end
        if row.outerDifference > ϵ && min_time_outer == 0
            min_time_outer = row.time
        end
    end

    return (min_time_inner, min_time_outer)
end

function get_forecasting_horzon(;)
    σ =  [10.0^(-1*i) for i in 12:-2:1]   # what is the measurement error sd?
    Nᵣ = 5000     # number of replicates
    ϵ = [0.000001]    # how far away counts as a correct prediction?
    pend = SameLengthSameMass()

    mean_diff_df = DataFrame(sigma=[], epsilon=[], time=[], outerDifference=[], innerDifference=[])

    nsteps = 3000
    df = DataFrame(sigma=[], epsilon=[], replicate=[],inner_time_until_fail=[], outer_time_until_fail=[])
    @showprogress for e ∈ ϵ
        for s ∈ σ    
            # 
            in_distance_matrix::Matrix{Float64} = zeros(nsteps+1, Nᵣ)
            out_distance_matrix::Matrix{Float64} = zeros(nsteps+1, Nᵣ)
            for r ∈ 1:Nᵣ
                replicate = run_double_pendulum(pend=pend, measurement_error_std_dev=s, output_num_timepoints=nsteps)
                (in_fail, out_fail) = time_until_forecast_fails(replicate, ϵ=e)
                push!(df.sigma, s)
                push!(df.epsilon, e)
                push!(df.replicate, r)
                push!(df.inner_time_until_fail, in_fail)
                push!(df.outer_time_until_fail, out_fail)
                for (t,row) in enumerate(eachrow(replicate))
                    in_distance_matrix[t,r] = row.innerDifference
                    out_distance_matrix[t,r] = row.outerDifference
                end
            end

            for t ∈ 1:nsteps
                push!(mean_diff_df.time, t)
                push!(mean_diff_df.outerDifference, mean(out_distance_matrix[t,:]) )
                push!(mean_diff_df.innerDifference, mean(in_distance_matrix[t,:]) )
                push!(mean_diff_df.sigma, s)
                push!(mean_diff_df.epsilon, e)
            end
        end
    end
    CSV.write("forecasting_horizon.csv", df)
  #  CSV.write("trajectories.csv", mean_diff_df)
end

get_forecasting_horzon()






anim = @animate for r ∈ 1:nrow(forecasted_trajectory_in_cartesian)
    x = [forecasted_trajectory_in_cartesian.y1[r], forecasted_trajectory_in_cartesian.y2[r]]
    y =  [-1*forecasted_trajectory_in_cartesian.x1[r],  -1*forecasted_trajectory_in_cartesian.x2[r]]
    scatter([0], [0], lims=((-2,2)), color="black", markersize=5)
    scatter!(x,y,markersize=10)
end

gif(anim, "outtest.gif")