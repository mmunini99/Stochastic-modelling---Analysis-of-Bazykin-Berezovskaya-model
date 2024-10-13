module PhasePortraits


include("BazykinBerezovskayaDeterministicPopulationModelClass.jl")
include("BazykinBerezovskayaStochasticPopulationModelClass.jl")
import .BazykinBerezovskayaDeterministicPopulationModelClass as BBDPM
import .BazykinBerezovskayaStochasticPopulationModelClass as BBSPM

export StorePhaseDeterministic, StorePhaseStochastic



function StorePhaseDeterministic(param_r::Float64,
                                 param_l::Float64, 
                                 param_k::Float64, 
                                 param_m::Float64, 
                                 param_T::Float64, 
                                 param_N::Int, 
                                 array_prey_dens::Vector{Float64},
                                 array_pred_dens::Vector{Float64})


    @assert length(array_prey_dens) == length(array_pred_dens) "the sequence of initial value must be the same for prey and predator"


    tensor = zeros(Float64, Int(param_N+1), 3, Int(length(array_prey_dens)))

    system = BBDPM.BazykinBerezovskayaDeterministicPopulationModel(r=param_r, l=param_l, k=param_k, m=param_m);

    counter = Int(1)

    for ind in range(1, length(array_prey_dens))

        tr_prey, tr_pred, time_ite = BBDPM.simulate_deterministic_traj(system,  array_prey_dens[ind], array_pred_dens[ind], param_T, param_N);
        

       
        tensor[:, 1, counter] = time_ite
        tensor[:, 2, counter] = tr_prey
        tensor[:, 3, counter] = tr_pred

        counter += 1

    end



    return tensor

end



function StorePhaseStochastic(param_r::Float64,
                             param_l::Float64, 
                             param_k::Float64, 
                             param_m::Float64,
                             param_a::Float64,
                             param_b::Float64, 
                             param_T::Float64, 
                             param_N::Int, 
                             array_prey_dens::Vector{Float64},
                             array_pred_dens::Vector{Float64})


    @assert length(array_prey_dens) == length(array_pred_dens) "the sequence of initial value must be the same for prey and predator"


    tensor = zeros(Float64, Int(param_N+1), 3, Int(length(array_prey_dens)))

    system = BBSPM.BazykinBerezovskayaStochasticPopulationModel(r=param_r, l=param_l, k=param_k, m=param_m, a=param_a, b=param_b);

    counter = Int(1)

    for ind in range(1, length(array_prey_dens))

        tr_prey, tr_pred, time_ite = BBSPM.simulate_stochastic_traj(system,  array_prey_dens[ind], array_pred_dens[ind], param_T, param_N);



        tensor[:, 1, counter] = time_ite
        tensor[:, 2, counter] = tr_prey
        tensor[:, 3, counter] = tr_pred

        counter += 1

    end



    return tensor

end






end