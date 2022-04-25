import os
from braket.jobs import save_job_result
from braket.jobs.metrics import log_metric
import json
import networkx as nx
import numpy as np
import pandas as pd
import math
import glob
import minorminer
from dwave_qbsolv import QBSolv
from dwave.system.composites import FixedEmbeddingComposite
from braket.aws import AwsDevice
from braket.ocean_plugin import BraketSampler, BraketDWaveSampler
from dwave.system.composites import EmbeddingComposite
import dimod
import random
import time
from qiskit.algorithms.optimizers import SPSA
from IPython.display import display, clear_output

def get_structures(pks, size, input_dir):
    subdirectory = input_dir +'/data'+'/'+pks+'/'+size
    bprna = []
    fasta = [f for f in os.listdir(subdirectory) if f.endswith('.fasta.txt')]
    for f in fasta:
        bprna.append(subdirectory+"/"+f.split(".")[0])
    return bprna

def actual_stems(seq_ss, seq_ps):
    
    with open(seq_ss) as file:
        ss_lines = file.readlines()
    
    with open(seq_ps) as file:
        ps_lines = file.readlines()
    
    rna = ps_lines[1]
    
    stems_actual = []

    sip = False                       # stem in progress?
    sl = 0                            # stem length
    last_line = [0, 0, 0, 0, 0, 0]    # initiate last line

    for i in range(0, len(ss_lines)):
        line = ss_lines[i].strip().split()
        
        if (int(line[4]) != 0 and sip == False):
            sip = True
            temp = [int(line[0]), int(line[4])]
            sl += 1
        elif (int(line[4]) != 0 and sip == True and (int(last_line[4])-int(line[4]) == 1)):
            sl += 1
        elif (int(line[4]) == 0 and sip == True):
            sip = False
            temp.append(sl)
            if temp[1] > temp[0]:
                stems_actual.append(temp)
            sl = 0
        elif ((int(last_line[4])-int(line[4]) != 1) and int(last_line[4]) != 0  and sip == True):
            temp.append(sl)
            if temp[1] > temp[0]:
                stems_actual.append(temp)
            temp = [int(line[0]), int(line[4])]
            sl = 1
        
        last_line = line
        
    return stems_actual

def potential_stems(seq_ps):
    
    with open(seq_ps) as file:
        lines = file.readlines()
    
    rna = lines[1]
    
    matrix = np.zeros((len(rna),len(rna)))
    for diag in range(0, len(matrix)):
        for row in range(0, len(matrix)-diag):
            col = row + diag
            base1 = rna[row]
            base2 = rna[col]
            if row != col:
                if ((base1 == ("A" or "a")) and (base2 == ("U" or "u"))) or ((base1 == ("U" or "u")) and (base2 == ("A" or "a"))) or ((base1 == ("G" or "g")) and (base2 == ("U" or "u"))) or ((base1 == ("U" or "u")) and (base2 == ("G" or "g"))) or ((base1 == ("G" or "g")) and (base2 == ("C" or "c"))) or ((base1 == ("C" or "c")) and (base2 == ("G" or "g"))):
                    matrix[row][col] = 1
    
    stems_potential = []
    mu = 0

    for row in range(0, len(matrix)):
        for col in range (row, len(matrix)):
            if row != col:
                if matrix[row][col] == 1:
                    temp_row = row
                    temp_col = col
                    stem = [row+1,col+1,0]
                    length = 0
                    while (matrix[temp_row][temp_col] != 0) and (temp_row != temp_col):
                        length+=1
                        temp_row+=1
                        temp_col-=1
                        if length >= 2 and col-row-2*length >= 3:
                            stem[2] = length
                            stems_potential.append(stem.copy())
                    if length > mu:
                        mu = length
    
    return [stems_potential, mu, rna, len(rna)]

def potential_pseudoknots(stems_potential):

    pseudoknots_potential = []

    for i in range(len(stems_potential)):
        for j in range(i + 1, len(stems_potential)):
            
            stem1 = stems_potential[i]
            stem2 = stems_potential[j]
    
            i_a = stem1[0]
            j_a = stem1[1]
            i_b = stem2[0]
            j_b = stem2[1]
    
            pseudoknot = [i, j, 0]
    
            if (i_a < i_b and i_b < j_a and j_a < j_b) or (i_b < i_a and i_a < j_b and j_b < j_a):
        
                pseudoknot[2] = 1
    
            pseudoknots_potential.append(pseudoknot)
            
    return pseudoknots_potential

def potential_overlaps(stems_potential):
    
    overlaps_potential = []
    overlap_penalty = 1e6

    for i in range(len(stems_potential)):
        for j in range(i+1, len(stems_potential)):
    
            stem1 = stems_potential[i]
            stem2 = stems_potential[j]
    
            overlap = [i, j, 0]
    
            stem1_cspan1 = set(range(stem1[1]-stem1[2]+1, stem1[1]+1))
            stem2_cspan1 = set(range(stem2[1]-stem2[2]+1, stem2[1]+1))
            
            stem1_cspan2 = set(range(stem1[0], stem1[0]+stem1[2]))
            stem2_cspan2 = set(range(stem2[0], stem2[0]+stem2[2]))
    
            if (len(stem1_cspan1 & stem2_cspan1) != 0) or (len(stem1_cspan2 & stem2_cspan2) != 0)  or (len(stem1_cspan1 & stem2_cspan2) != 0) or (len(stem1_cspan2 & stem2_cspan1) != 0):
        
                overlap[2] = overlap_penalty
        
            overlaps_potential.append(overlap)
            
    return overlaps_potential


def model(stems_potential, overlaps_potential, pseudoknots_potential, cl, cb, delta):
    
    L = {}
    Q = {}
    k = 0
    
    for i in range(0, len(stems_potential[0])):
        L[str(i)] = cl*((stems_potential[0][i][2]-stems_potential[1])**2)-cb*(stems_potential[0][i][2]**2)
        for j in range(i+1, len(stems_potential[0])):
            Q[(str(i), str(j))] = -(2*cb*stems_potential[0][i][2]*stems_potential[0][j][2]*(delta*pseudoknots_potential[k][2] + (1-pseudoknots_potential[k][2]))) + overlaps_potential[k][2]
            k += 1
            
    return L, Q

def energy(stems_actual, mu, cl, cb, delta):
    
    k = 0
    
    pseudoknots_actual = potential_pseudoknots(stems_actual)
    cost = 0

    for i in range(0, len(stems_actual)):
        cost += cl*((stems_actual[i][2]-mu)**2)-cb*(stems_actual[i][2]**2)
        for j in range(i+1, len(stems_actual)):
            cost += -2*cb*stems_actual[i][2]*stems_actual[j][2]*(delta*pseudoknots_actual[k][2] + (1-pseudoknots_actual[k][2]))
            k += 1
    
    return cost

def evaluation_1(stems_actual, stems_potential):
    
    bp_actual = []
    bp_predicted = []

    for i in range(0, len(stems_actual)):
        for j in range(0, stems_actual[i][2]):
            bp_actual.append((stems_actual[i][0]+j, stems_actual[i][1]-j))
        
    for i in range(0, len(stems_potential)):
        for j in range(0, stems_potential[i][2]):
            bp_predicted.append((stems_potential[i][0]+j, stems_potential[i][1]-j))
            
    TP = 0    # number of correctly identified base pairs
    FP = 0    # number of the predicted base pairs missing from the known structure
    FN = 0    # number of non-predicted base pairs present in the known structure

    for i in range(0, len(bp_predicted)):
        if bp_predicted[i] in bp_actual:
            TP += 1
        else:
            FP += 1

    for i in range(0, len(bp_actual)):
        if bp_actual[i] not in bp_predicted:
            FN += 1
            
    if TP+FP != 0:
        ppv = TP/(TP+FP)
    else:
        ppv = 0
    if TP+FN != 0:
        sensitivity = TP/(TP+FN)
    else:
        sensitivity = 0
    
    return [ppv, sensitivity]

# function to compare actual and predicted structure based on comparison of bases involved in pairing:

def evaluation_2(stems_actual, stems_predicted):
    
    pb_actual = []
    pb_predicted = []

    for i in range(0, len(stems_actual)):
        for j in range(0, stems_actual[i][2]):
            pb_actual.append(stems_actual[i][0]+j)
            pb_actual.append(stems_actual[i][1]-j)
        
    for i in range(0, len(stems_predicted)):
        for j in range(0, stems_predicted[i][2]):
            pb_predicted.append(stems_predicted[i][0]+j)
            pb_predicted.append(stems_predicted[i][1]-j)
            
    TP = 0    # number of correctly identified bases that are paired
    FP = 0    # number of the predicted paired bases missing from the known structure
    FN = 0    # number of non-predicted paired bases present in the known structure

    for i in range(0, len(pb_predicted)):
        if pb_predicted[i] in pb_actual:
            TP += 1
        else:
            FP += 1

    for i in range(0, len(pb_actual)):
        if pb_actual[i] not in pb_predicted:
            FN += 1
            
    if TP+FP != 0:
        ppv = TP/(TP+FP)
    else:
        ppv = 0
    if TP+FN != 0:
        sensitivity = TP/(TP+FN)
    else:
        sensitivity = 0
    
    return [ppv, sensitivity]

def MCC(stems_actual, stems_predicted, stems_potential):
    
    bp_actual = []
    bp_predicted = []
    bp_potential = []

    for i in range(0, len(stems_actual)):
        for j in range(0, stems_actual[i][2]):
            bp_actual.append((stems_actual[i][0]+j, stems_actual[i][1]-j))
        
    for i in range(0, len(stems_predicted)):
        for j in range(0, stems_predicted[i][2]):
            bp_predicted.append((stems_predicted[i][0]+j, stems_predicted[i][1]-j))
            
    for i in range(0, len(stems_potential)):
        for j in range(0, stems_potential[i][2]):
            bp_potential.append((stems_potential[i][0]+j, stems_potential[i][1]-j))
            
    TP = 0    # number of correctly identified base pairs
    TN = 0    # number of correctly excluded base pairs
    FP = 0    # number of the predicted base pairs missing from the known structure
    FN = 0    # number of non-predicted base pairs present in the known structure

    for i in range(0, len(bp_predicted)):
        if bp_predicted[i] in bp_actual:
            TP += 1
        else:
            FP += 1

    for i in range(0, len(bp_actual)):
        if bp_actual[i] not in bp_predicted:
            FN += 1
            
    for i in range(0, len(bp_potential)):
        if (bp_potential[i] not in bp_predicted) and (bp_potential[i] not in bp_actual):
            TN += 1
    
    if TP + FP != 0 and TP + FN != 0 and TN + FP != 0 and TN + FN != 0:
        score = (TP*TN-FP*FN)/np.sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
    else:
        score = -1
        
    return score

def calculate_cost_function(expectation_values, target_values):
    product_zt = expectation_values*target_values
    all_costs = ((1-product_zt)/2)**2
    return all_costs

def get_random_split(bpRNAs, seed):
    random.Random(seed).shuffle(bpRNAs)
    return bpRNAs[:15], bpRNAs[15:]

def spsa_optimizer_callback(nb_fct_eval, params, fct_value, stepsize, step_accepted, train_history):
    print("In callback")
    train_history.append((nb_fct_eval,params,fct_value))
    clear_output(wait=True)
    display(f'evaluations : {nb_fct_eval} loss: {fct_value:0.4f}')

def optimize_params(optimizer, hyper_params, inital_point, full_bprna, solver, num_repeats, solver_limit, num_reads, seed):
    past_params = {}
    target_value = 1
    bp_metric = []
    a_stems = {}
    p_stems = {}
    p_psudoknots = {}
    p_overlaps = {}
    stems_f = {}
    for bprna in full_bprna:
        bprna_id = bprna.split("/")[8]
        print(bprna_id)
        fasta_file = bprna + ".fasta.txt"
        ct_file = bprna + ".ct.txt" 
        p_stems[bprna_id] = potential_stems(fasta_file)
        p_psudoknots[bprna_id] = potential_pseudoknots(p_stems[bprna_id][0])
        p_overlaps[bprna_id] = potential_overlaps(p_stems[bprna_id][0])
        a_stems[bprna_id] = actual_stems(ct_file, fasta_file)
    print("Finished pre-processing")
    
    def cost_function(hyper_params):
        print(hyper_params)
        if str(hyper_params) in past_params:
            print("Previous value")
            return past_params[str(hyper_params)]
        else:
            cl = hyper_params[0]
            cb = hyper_params[1]
            delta = hyper_params[2]
            problems = {}
            a_energies = {}
            for bprna in full_bprna:
                print(bprna)
                bprna_id = bprna.split("/")[8]
                a_energies[bprna_id] = energy(a_stems[bprna_id], p_stems[bprna_id][1], cl, cb, delta)
                md = model(p_stems[bprna_id], p_overlaps[bprna_id], p_psudoknots[bprna_id], cl, cb, delta)
                problems[bprna_id] = dimod.BinaryQuadraticModel(md[0], md[1], vartype = 'BINARY', offset = 0.0) 
            print("Finished creating the models")

            for key, value in problems.items():
                t0 = time.time()
                sampleset = QBSolv().sample_qubo(
                    value.to_qubo()[0], solver=solver, num_repeats=num_repeats,solver_limit=solver_limit, num_reads=num_reads, seed=seed
                    )
                print(f"Model: {key} spent {round(time.time()-t0, 2)} seconds for {num_repeats} repeats")
                for datum in sampleset.data(['sample', 'energy', 'num_occurrences']):
                    results_hybrid = datum.sample
                    predicted_energy = datum.energy
        
                f_stems = []

                for j in range(0, len(results_hybrid)):
                    if results_hybrid[str(j)] == 1:
                        f_stems.append(p_stems[key][0][j])
                    
                stems_f[key] = ([f_stems, predicted_energy])
                bp_metric.append(MCC(a_stems[key], stems_f[key][0], p_stems[key][0]))
                print("Metric:",bp_metric[len(bp_metric)-1])
            print("Finished running the models")
            all_costs = calculate_cost_function(np.array(bp_metric), target_value)
            cost = np.sum(all_costs)/len(all_costs)
            print(cost)
            past_params[str(hyper_params)] = cost
            return cost
    model_values, loss, nfev = optimizer.optimize(len(hyper_params), cost_function, initial_point=inital_point)
    return model_values, loss, nfev 

def main():
    print(os.environ)
    input_dir = os.environ["AMZN_BRAKET_INPUT_DIR"]
    hp_file = os.environ["AMZN_BRAKET_HP_FILE"]
    job_name = os.environ["AMZN_BRAKET_JOB_NAME"]
    s3_bucket = os.environ["AMZN_BRAKET_OUT_S3_BUCKET"]
    device_arn = os.environ["AMZN_BRAKET_DEVICE_ARN"]

    with open(hp_file, "r") as f:
        hyperparams = json.load(f)
    print(hyperparams)

    solver_limit = int(hyperparams["solver_limit"])
    num_repeats = int(hyperparams["num_repeats"])
    num_reads = int(hyperparams["num_reads"])
    seed = int(hyperparams["seed"])


    s3_task_prefix = f"jobs/{job_name}/tasks"
    s3_folder = (s3_bucket, s3_task_prefix)

    system = BraketDWaveSampler(s3_folder, device_arn)
    G_sub = nx.complete_graph(solver_limit)
    embedding = minorminer.find_embedding(G_sub.edges, system.edgelist)
    solver = FixedEmbeddingComposite(system, embedding)

    size = "s"  
    bprna_wPKs = get_structures("wPKs", size, input_dir)
    bprna_woutPKs = get_structures("woutPKs", size, input_dir)
    training_data = bprna_wPKs + bprna_woutPKs
    print(len(training_data))
    #testing_data = ve_wPKs + ve_woutPKs
    print("Finished data set creation")
    initial_point = [1, 1, 1]
    model_params = [1, 1, 1]
    train_history = []
    optimizer = SPSA(maxiter=30, callback=lambda n, p, v, ss, sa: spsa_optimizer_callback(n, p, v, ss, sa, train_history))
    model_values, loss, nfev = optimize_params(optimizer, model_params, initial_point, training_data, solver, num_repeats, solver_limit, num_reads, seed)
    print("Final values:",model_values)
    print('Save results')
    save_job_result({"final_values": str(model_values), "loss":str(loss), "hyperparams": str(hyperparams)})

if __name__ == "__main__":
    main()