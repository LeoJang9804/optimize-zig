# optimize-zig
indeter optimization to maximize sensitivity 

genotype A = [a_1', a_2', ..., a_20']
a_i' = a_i/10
a_i = distace between axis and surface at z = 10*i/20 mm

for gen in num_gen:
    cand_genotypes = [genotype_1, genotype_2, ..., genotype_n]

    for genType in cand_genotypes:
        In abaqus python: 
            make zig from information in genType
            perform 3 simulation for each parameters (ref, ref+delta_p, ref-delta_p)
        calculate reward (fitness)
            sensitivity of each parameter p : S_p = argmax_{i}(|F_i(theta+delta_p)-F_i(theta)|)
            reward = mean(S_p)/(1+std(S_p))
        delete relevant files for memory
    
    select parents which have top k fitness
    delete some genTypes by random selection
    make new children genTypes from parents by cross and mutation algorithm 
    
save each genType and fitness information 

- [] create function in abaqus python that make and execute simulation and control zig shape
- [] create function that can execute several simulation to calculate reward 
- [] create abaqus python file that can calculate rewards from simulation results
- [] wrap above three function in one function 
- [] perform genetic algorithm

concerns: optimization time 
This concerns can be mitigated by multiprocessing. refer below url to learn how to use multiprocessing in PyGAD
https://hackernoon.com/how-genetic-algorithms-can-compete-with-gradient-descent-and-backprop-9m9t33bq