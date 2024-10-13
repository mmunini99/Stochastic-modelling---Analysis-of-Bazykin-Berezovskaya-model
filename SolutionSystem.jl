module SimulationStochasticEnding



include("C:\\Users\\matte\\Documents\\STO MOD AND SIM\\Progetto Esame\\BazykinBerezovskayaStochasticPopulationModelClass.jl")
include("C:\\Users\\matte\\Documents\\STO MOD AND SIM\\Progetto Esame\\BazykinBerezovskayaDeterministicPopulationModelClass.jl")
import .BazykinBerezovskayaStochasticPopulationModelClass as BBSPM
import .BazykinBerezovskayaDeterministicPopulationModelClass as BBDPM
using Random

export GenerateSolution, normalized_simulative_pdf,  probability_of_exiting, link_analysis_lm_det



function probability_of_exiting(param_r::Float64,
                                param_l::Float64, 
                                param_k::Float64, 
                                param_m::Float64,
                                param_a::Vector{Float64}, 
                                param_b::Vector{Float64}, 
                                param_T::Float64, 
                                param_N::Int, 
                                initial_prey::Float64,
                                initial_predator::Float64,
                                number_of_sim::Int)

    matrix_prob = zeros(Float64, length(param_a), length(param_b))

    for ind_a in 1:length(param_a)
        for ind_b in 1:length(param_b)

            system = BBSPM.BazykinBerezovskayaStochasticPopulationModel(r = param_r, l = param_l, k = param_k, m = param_m, a = param_a[ind_a], b = param_b[ind_b])

            counter = 0
        

            for _ in 1:number_of_sim
    
                
    
                tr_prey, tr_pred, time_ite = BBSPM.simulate_stochastic_traj(system,  initial_prey, initial_predator, param_T, param_N);

                if tr_prey[end] <= 0.005 && tr_pred[end] <= 0.005

                    counter += 1

                end                                                    
    
               
            end

            matrix_prob[ind_a, ind_b] += counter/number_of_sim

        end

    end


    return matrix_prob



    

end


function GenerateSolution(param_r::Vector{Float64},
                          param_l::Vector{Float64}, 
                          param_k::Vector{Float64}, 
                          param_m::Vector{Float64},
                          param_a::Vector{Float64}, 
                          param_b::Vector{Float64}, 
                          param_T::Vector{Float64}, 
                          param_N::Int, 
                          initial_prey::Float64,
                          initial_predator::Float64,
                          choice::String,
                          number_of_sim::Int)



    config_sim_params = Dict(
    "r" => param_r,
    "l" => param_l,
    "k" => param_k,
    "m" => param_m,
    "a" => param_a,
    "b" => param_b)

    
    @assert choice in ["r", "l", "k", "m", "a", "b"] "You have to choose a simulation over one of these parameters: r, l, k, m, a and b"
    @assert all(key -> key == choice || length(config_sim_params[key]) == 1, keys(config_sim_params)) "A vector of length 1 for the fixed parameters"


    tensor = zeros(Float64, number_of_sim * length(param_T), 2, Int(length(config_sim_params[choice]))) #store the end of trajectories, 3rd D for the parameter to move along

    for param_choice_ind in 1:length(config_sim_params[choice])

        system = BBSPM.BazykinBerezovskayaStochasticPopulationModel(
                                                                    r=choice == "r" ? Float64(config_sim_params[choice][param_choice_ind]) : Float64(param_r[1]),
                                                                    l=choice == "l" ? Float64(config_sim_params[choice][param_choice_ind]) : Float64(param_l[1]),
                                                                    k=choice == "k" ? Float64(config_sim_params[choice][param_choice_ind]) : Float64(param_k[1]),
                                                                    m=choice == "m" ? Float64(config_sim_params[choice][param_choice_ind]) : Float64(param_m[1]),
                                                                    a=choice == "a" ? Float64(config_sim_params[choice][param_choice_ind]) : Float64(param_a[1]),
                                                                    b=choice == "b" ? Float64(config_sim_params[choice][param_choice_ind]) : Float64(param_b[1]))

        counter = Int(1)
        

        for n_sim in 1:number_of_sim

            for t_sim in param_T

                tr_prey, tr_pred, time_ite = BBSPM.simulate_stochastic_traj(system,  initial_prey, initial_predator, t_sim, param_N);


                tensor[counter, 1, param_choice_ind] = tr_prey[end]
                tensor[counter, 2, param_choice_ind] = tr_pred[end]
        
                counter += 1
                                                                

            end
        end

    end


    return tensor


end


function total_var_dist(P::Array{Float64, 2}, Q::Array{Float64, 2})

    @assert size(P)[1] == size(Q)[1] "The length of the matrixes of the two distributions must be the same"
    @assert size(P)[2] == size(Q)[2] "The length of the matrixes of the two distributions must be the same"

    dist = sum(abs.(P .- Q))

    return dist / 2.0

end


function normalized_simulative_pdf(param_r::Float64,
                                   param_l::Float64, 
                                   param_k::Float64, 
                                   param_m::Float64, 
                                   param_T::Float64, 
                                   param_N::Int,
                                   omega::Vector{Float64},
                                   n_sim::Int,
                                   n_bins::Int)

      
    bins = zeros(Float64, n_bins, n_bins) #matrix for the pdf

    m = (omega[2]-omega[1])/n_bins #unique since x and y bounded between 0 and 1

    for _ in 1:n_sim

        system = BBDPM.BazykinBerezovskayaDeterministicPopulationModel(r=param_r, l=param_l, k=param_k, m=param_m)

        rng = MersenneTwister()

        sampled_prey = rand(rng) * (omega[2] - omega[1]) + omega[1]
        sampled_predator = 1 - sampled_prey

        tr_prey, tr_pred, time_ite = BBDPM.simulate_deterministic_traj(system,  sampled_prey, sampled_predator, param_T, param_N);

        position_i = 1
        position_j = 1

        for i in 1:n_bins

            if tr_prey[end] < (m * (i ) + omega[1])
                position_i += i
                break
            end
        end

        for j in 1:n_bins

            if tr_pred[end] < (m * (j ) + omega[1])
                position_i += j
                break
            end
        end

        if position_i > n_bins
            position_i = n_bins
        end

        if position_j > n_bins
            position_j = n_bins
        end
     

        bins[position_i, position_j] += 1.0
        


    end

    bin_widths = fill(m, n_bins)
    # Calculate the total area
    # Calculate total area using a comprehension
    total_area = sum(bins[i, j] * m * m for i in 1:n_bins, j in 1:n_bins)


    # Normalize the bins
    t_bins = bins ./ total_area  # Element-wise division


    return t_bins

end


function link_analysis_lm_det(param_r::Float64,
                              param_l::Vector{Float64}, 
                              param_k::Float64, 
                              param_m::Vector{Float64}, 
                              param_T::Float64, 
                              param_N::Int,
                              omega::Vector{Float64},
                              n_sim::Int,
                              n_bins::Int)

    matrix_diff = zeros(Float64, length(param_l)^2, length(param_m)^2)
    matrix_pdf = Array{Array{Float64, 2}, 2}(undef, length(param_l), length(param_m))
    sets_names = String[]
    

    for ind_l in 1:length(param_l)
        for ind_m in 1:length(param_m)
            push!(sets_names, "$ind_l-$ind_m")
            matrix_pdf[ind_l, ind_m] = normalized_simulative_pdf(param_r, param_l[ind_l], param_k, param_m[ind_m], param_T, param_N, omega, n_sim, n_bins)
        end
    end
    

    for ind_l in 1:length(param_l)^2
        for ind_m in 1:length(param_m)^2

            if ind_m >= ind_l+1

                matrix_diff[ind_l,ind_m] += total_var_dist(
                    matrix_pdf[div(ind_l-1, length(param_l))+1, (ind_l-1) % length(param_l) + 1], 
                    matrix_pdf[div(ind_m-1, length(param_m))+1, (ind_m-1) % length(param_m) + 1]
                )

            end


        end

    end

    return matrix_diff

end


end