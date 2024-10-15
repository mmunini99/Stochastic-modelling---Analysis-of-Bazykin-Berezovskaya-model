module grid_search_tuning


include("C:\\Users\\matte\\Documents\\STO MOD AND SIM\\Progetto Esame\\BazykinBerezovskayaDeterministicPopulationModelClass.jl")
import .BazykinBerezovskayaDeterministicPopulationModelClass as BBDPM
using Statistics


export grid_search_parameter_BB_model


function grid_search_parameter_BB_model(train_df, list_r, list_l, list_k, list_m)

    loss_min = 1000000
    param_best = [0.01, 0.01, 0.01, 0.01]

    for r_tent in list_r
        for l_tent in list_l
            for k_tent in list_k
                for m_tent in list_m

                    if l_tent < k_tent

                        det_system = BBDPM.BazykinBerezovskayaDeterministicPopulationModel(r=r_tent, l=l_tent, k=k_tent, m=m_tent);

                        tr_x, tr_y, time_sim = BBDPM.simulate_deterministic_traj(det_system, train_df[1, :prey], train_df[1, :pred], 49.0, 49000 );

                        tr_x = tr_x[collect(range(1, stop=length(tr_x), step=1000))]  #365*1000
                        tr_y = tr_y[collect(range(1, stop=length(tr_y), step=1000))]

                        loss = sqrt(mean((tr_x .- train_df[:, "prey"]).^2)) + sqrt(mean((tr_y .- train_df[:, "pred"]).^2))

                        if loss < loss_min

                            loss_min = loss

                            param_best = [r_tent, l_tent, k_tent, m_tent]

                        end
                    end

                end
            end
        end
    end

    return param_best

end

end