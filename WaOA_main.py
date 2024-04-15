import random
from  SBMA_MatingFunction import *
from  SB_SupportFunction import *
from evaluate_all_sols import *
from SBMA_Fucntion import *
import  itertools
import time
value_heavy = 40
random.seed(2222)



def walrus_optimization_algorithm(name_path_input, num_gen, pop_size):
    df_item_pool, df_item_sas_random = read_input(name_path_input)
    num_item = df_item_pool.shape[0]

    E_all = list(itertools.permutations(range(num_item), 2))
    sub_E_list = [[arc for arc in E_all if arc[0] == item or arc[1] == item] for item in range(num_item)]

    # สร้าง set ของ item ที่ถือว่าเป็น item หนัก
    num_item = len(df_item_pool)
    heavy_item_set = set(df_item_pool[df_item_pool['weight'] >= value_heavy].index)

    # Initialize gbest and pbest values for males and females
    gbest_value, gbest_sol, gbest_arc_sol_cut = 100000, [], []
    pbest_value = [100000] * pop_size
    pbest_sol= [[] for _ in range(pop_size)]
    pbest_arc_sols = [[] for _ in range(pop_size)]
    pbest_arc_sols_cut= [[] for _ in range(pop_size)]

    evaluations = []
    cur_sols = []
    cur_sols_value = []
    cur_arc_sols = []
    arc_sols_cut = []
    velocity_dict = []

    walrus_population = [[i for i in range(num_item)] for _ in range(pop_size)]
    for walrus in walrus_population:
        random.shuffle(walrus)
        evaluation = evaluate_all_sols_check(walrus, df_item_pool, heavy_item_set, name_path_input)
        evaluations.append(evaluation)
        cur_sols.append(extract_and_flatten(evaluation))
        cur_sols_value.append(evaluation[1])
        cur_arc_sols.append(sol_from_list_to_arc(cur_sols[-1]))
        arc_sols_cut.append(cut_arc_sol(cur_arc_sols[-1]))
        velocity_dict.append(init_velocity_sol(arc_sols_cut[-1]))

    gbest_each_gen = []
    for gen in range(num_gen):
        a = random.uniform(0.1, 0.3)
        b = random.uniform(0.7, 2)
    # PHASE 1 FEEDING STRATEGY (EXPLORATION)
        # 1.1. Update the strongest walrus in the population
        strongest_walrus = cur_sols_value.index(min(cur_sols_value))
        if cur_sols_value[strongest_walrus] <= gbest_value:
            gbest_value = cur_sols_value[strongest_walrus]
            gbest_sol = cur_sols[strongest_walrus]
            gbest_arc_sol_cut = arc_sols_cut[strongest_walrus]

        # 1.2. Update the personal bests of all walruses
        for i in range(pop_size):
            if cur_sols_value[i] <= pbest_value[i]:
                pbest_value[i] = cur_sols_value[i]
                pbest_sol[i] = cur_sols[i]
                pbest_arc_sols[i] = cur_arc_sols[i]
                pbest_arc_sols_cut[i] = arc_sols_cut[i]

        # 1.3 Feeding strategy for each walrus in the population (EXPLORATION)
        for i in range(pop_size):
            walrus_phase1_velocity = check_velocity_inconsistency(add_velocity(velocity_dict[i],(coef_times_position(random.uniform(a, b),position_minus_position(pbest_arc_sols_cut[i], arc_sols_cut[i])))))
            walrus_phase1_position = sol_position_update(creat_cut_set(walrus_phase1_velocity,alpha_calculation(gen,num_gen)),arc_sols_cut[i],sub_E_list,cur_sols[i][0],pbest_sol[i][0],gbest_sol[0])[0]

            # 1.3.1 Update the walrus if the new position is better
            cur_sols[i] = walrus_phase1_position
            cur_arc_sols[i] = sol_from_list_to_arc(cur_sols[i])
            arc_sols_cut[i] = cut_arc_sol(cur_arc_sols[i])
            velocity_dict[i] = walrus_phase1_velocity

    # PHASE 2 MIGRATION
        # 2.1 Choose a random walrus as the migrant walrus and a random walrus as the receiver walrus in the population for migration process (MIGRATION)
            migrant_walrus = random.choice(range(pop_size))
            receiver_walrus = random.choice(range(pop_size))

            # 2.2 Migrate the migrant walrus to the receiver walrus
            walrus_phase2_velocity = check_velocity_inconsistency(add_velocity(velocity_dict[migrant_walrus],(coef_times_position(random.uniform(a,b),position_minus_position(arc_sols_cut[receiver_walrus],arc_sols_cut[migrant_walrus])))))
            walrus_phase2_position = sol_position_update(creat_cut_set(walrus_phase2_velocity,alpha_calculation(gen,num_gen)),arc_sols_cut[migrant_walrus],sub_E_list,cur_sols[migrant_walrus][0],pbest_sol[migrant_walrus][0],gbest_sol[0])[0]

            # 2.3 Update the receiver walrus if the new position is better
            cur_sols[receiver_walrus] =  walrus_phase2_position
            cur_arc_sols[receiver_walrus] = sol_from_list_to_arc(cur_sols[receiver_walrus])
            arc_sols_cut[receiver_walrus] = cut_arc_sol(cur_arc_sols[receiver_walrus])
            velocity_dict[receiver_walrus] = walrus_phase2_velocity

        # PHASE 3 PREDATOR-PREY INTERACTION (EXPLOITATION)
            # 3.1 Choose a random walrus as the predator
            predator = random.choice(range(pop_size))

            # 3.2 Escape from the predator if the predator is stronger than the current walrus
            if cur_sols_value[predator] <= cur_sols_value[strongest_walrus]:
                walrus_phase3_velocity = check_velocity_inconsistency(add_velocity(velocity_dict[strongest_walrus] ,(coef_times_position(random.uniform(a, b),position_minus_position(pbest_arc_sols_cut[predator], arc_sols_cut[strongest_walrus])))))
                walrus_phase3_position = sol_position_update(creat_cut_set(walrus_phase3_velocity,alpha_calculation(gen,num_gen)),arc_sols_cut[strongest_walrus],sub_E_list,cur_sols[strongest_walrus][0],pbest_sol[predator][0],gbest_sol[0])[0]

                # 3.2.1 Update the walrus if the new position is better
                cur_sols[strongest_walrus] = walrus_phase3_position
                cur_arc_sols[strongest_walrus] = sol_from_list_to_arc(cur_sols[strongest_walrus])
                arc_sols_cut[strongest_walrus] = cut_arc_sol(cur_arc_sols[strongest_walrus])
                velocity_dict[strongest_walrus] = walrus_phase3_velocity

            # 3.3 Fight against the predator if the predator is weaker than the current walrus
            else:
                walrus_phase3_velocity = check_velocity_inconsistency(add_velocity(velocity_dict[predator], (coef_times_position(random.uniform(a, b),position_minus_position(pbest_arc_sols_cut[strongest_walrus], arc_sols_cut[predator])))))
                walrus_phase3_position = sol_position_update(creat_cut_set(walrus_phase3_velocity,alpha_calculation(gen,num_gen)),arc_sols_cut[predator],sub_E_list,cur_sols[predator][0],pbest_sol[strongest_walrus][0],gbest_sol[0])[0]

                # 3.3.1 Update the walrus if the new position is better
                cur_sols[predator] = walrus_phase3_position
                cur_arc_sols[predator] = sol_from_list_to_arc(cur_sols[predator])
                arc_sols_cut[predator] = cut_arc_sol(cur_arc_sols[predator])
                velocity_dict[predator] = walrus_phase3_velocity

        # PHASE 4 COMMUNICATION STRATEGY
            # 4.1 Choose a random walrus as the sender and a random walrus as the receiver in the population for communication process
            sender = random.choice(range(pop_size))
            receiver = random.choice(range(pop_size))

            # 4.2 Communicate the sender walrus to the receiver walrus
            walrus_phase4_velocity = check_velocity_inconsistency(add_velocity(velocity_dict[receiver],(coef_times_position(random.uniform(a, b),position_minus_position(arc_sols_cut[sender],arc_sols_cut[receiver])))))
            walrus_phase4_position = sol_position_update(creat_cut_set(walrus_phase4_velocity,alpha_calculation(gen,num_gen)),arc_sols_cut[receiver],sub_E_list,cur_sols[receiver][0],pbest_sol[sender][0],gbest_sol[0])[0]

            # 4.3 Update the receiver walrus if the new position is better
            cur_sols[receiver] = walrus_phase4_position
            cur_arc_sols[receiver] = sol_from_list_to_arc(cur_sols[receiver])
            arc_sols_cut[receiver] = cut_arc_sol(cur_arc_sols[receiver])
            velocity_dict[receiver] = walrus_phase4_velocity

            # Evaluate the new solutions
            evaluation = evaluate_all_sols_check(cur_sols[i], df_item_pool, heavy_item_set, name_path_input)
            evaluations[i] = evaluation
            cur_sols_value[i] = evaluation[1]
            cur_arc_sols[i] = sol_from_list_to_arc(cur_sols[i])
            arc_sols_cut[i] = cut_arc_sol(cur_arc_sols[i])

            # Update the best solution
            if cur_sols_value[i] <= gbest_value:
                gbest_value = cur_sols_value[i]
                gbest_sol = cur_sols[i]
                gbest_arc_sol_cut = arc_sols_cut[i]

            # Update the personal bests
            if cur_sols_value[i] <= pbest_value[i]:
                pbest_value[i] = cur_sols_value[i]
                pbest_sol[i] = cur_sols[i]
                pbest_arc_sols[i] = cur_arc_sols[i]
                pbest_arc_sols_cut[i] = arc_sols_cut[i]


        gbest_each_gen.append(gbest_value)

    # return gbest_value

        print('Generation:', gen, 'Best Solution:', gbest_value)
    batch_solution = evaluate_all_sols_check(gbest_sol, df_item_pool, heavy_item_set, name_path_input)
    print('Best solution:', batch_solution)



num_gen = 100
pop_size = 50
start_time = time.time()
name_path_input = '1R-20I-150C-2P'
walrus_optimization_algorithm(name_path_input, num_gen, pop_size)
end_time = time.time()
duration = end_time - start_time
print('Duration:', round(duration,2))