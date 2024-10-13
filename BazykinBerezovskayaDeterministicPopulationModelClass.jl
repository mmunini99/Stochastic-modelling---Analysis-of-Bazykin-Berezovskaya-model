module BazykinBerezovskayaDeterministicPopulationModelClass


export BazykinBerezovskayaDeterministicPopulationModel, simulate_deterministic_traj

using Random

struct BazykinBerezovskayaDeterministicPopulationModel
    r::Float64
    l::Float64
    k::Float64
    m::Float64


    function BazykinBerezovskayaDeterministicPopulationModel(; r::Float64, l::Float64, k::Float64, m::Float64)
        return new(r, l, k, m)
    end

end


function check_parameters(process::BazykinBerezovskayaDeterministicPopulationModel)
    # Check initial parameters
    if process.r <= 0
        throw(ArgumentError("The intrinsic growth rate must be greater than 0"))
    end

    if process.l <= 0
        throw(ArgumentError("The prey survival threshold must be greater than 0"))
    end

    if process.k <= 0
        throw(ArgumentError("The carrying capacity rate must be greater than 0"))
    end

    if process.m <= 0
        throw(ArgumentError("The mortality rate of the predators must be greater than 0"))
    end

    if process.l >= process.k
        throw(ArgumentError("The carrying capacity rate must be greater than prey survival threshold"))
    end

end

function check_inputs(process::BazykinBerezovskayaDeterministicPopulationModel, x0::Float64, y0::Float64, T::Float64, N::Int)
    # Check x0
    if x0 < 0
        throw(ArgumentError("Initial prey population must be positive"))
    end

    # Check y0
    if y0 < 0
        throw(ArgumentError("Initial predator population must be positive"))
    end

    # Check T
    if T < 0
        throw(ArgumentError("Time interval must have a positive length"))
    end

    # Check N
    if N <= 1
        throw(ArgumentError("The simulation must have at least two steps"))
    end
end

function rk4(process::BazykinBerezovskayaDeterministicPopulationModel, t_n::Float64, x_n::Float64, y_n::Float64, h::Float64)
    


    if h <= 0
        throw(ArgumentError("Given step is negative or 0"))
    end

    fun_prey = (x,y) -> process.r * x * (x - process.l) *(process.k - x) - x * y
    fun_pred = (x,y) -> y*(x - process.m)

    

    k_1_prey, k_1_pred = fun_prey(x_n,y_n), fun_pred(x_n,y_n)
    k_2_prey, k_2_pred = fun_prey((x_n + h * k_1_prey / 2), (y_n + h * k_1_pred / 2)), fun_pred((x_n + h * k_1_prey / 2), (y_n + h * k_1_pred / 2))
    k_3_prey, k_3_pred = fun_prey((x_n + h * k_2_prey / 2), (y_n + h * k_2_pred / 2)), fun_pred((x_n + h * k_2_prey / 2), (y_n + h * k_2_pred / 2))
    k_4_prey, k_4_pred = fun_prey((x_n + h * k_3_prey), (y_n + h * k_3_pred)), fun_pred((x_n + h * k_3_prey), (y_n + h * k_3_pred))

    return h * (k_1_prey + 2 * k_2_prey + 2 * k_3_prey + k_4_prey) / 6, h * (k_1_pred + 2 * k_2_pred + 2 * k_3_pred + k_4_pred) / 6
end


function simulate_deterministic_traj(process::BazykinBerezovskayaDeterministicPopulationModel, x0::Float64, y0::Float64, T::Float64, N::Int)
    check_parameters(process)
    # Check inputs
    check_inputs(process, x0, y0, T, N)

    if T == 0
        return [x0], [y0], [0.0]
    end

    # Setup step length and trajectory array
    h = T / N
    traj_prey = zeros(Float64, N + 1)
    traj_prey[1] = x0
    traj_pred = zeros(Float64, N + 1)
    traj_pred[1] = y0
    time = zeros(Float64, N + 1)
    time[1] = 0.0


    for i in 2:(N + 1)
        time[i] = (i - 1) * h
        dx, dy = rk4(process, time[i - 1], traj_prey[i - 1], traj_pred[i - 1], h)  # Deterministic part
        
        traj_prey[i] = traj_prey[i - 1] + dx
        traj_pred[i] = traj_pred[i - 1] + dy
        # lower bound condition of prey and predator density
        if traj_prey[i] < 0
            traj_prey[i] = 0.0
        end
        if traj_pred[i] < 0
            traj_pred[i] = 0.0
        end
    end

    return traj_prey, traj_pred, time
end




end