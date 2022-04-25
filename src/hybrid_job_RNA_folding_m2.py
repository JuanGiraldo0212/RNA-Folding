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

def SQ_energy(sq):
    sqe = 0
    if len(sq) > 1:
        for i in range(1, len(sq)):
            if sq[i] == "AU":
                if sq[i-1] == "AU": 
                    sqe += 0.9
                if sq[i-1] == "CG":
                    sqe += 2.2
                if sq[i-1] == "GC":
                    sqe += 2.1
                if sq[i-1] == "UA":
                    sqe += 1.1
                if sq[i-1] == "GU":
                    sqe += 0.6
                if sq[i-1] == "UG":
                    sqe += 1.4
            if sq[i] == "CG":
                if sq[i-1] == "AU": 
                    sqe += 2.1
                if sq[i-1] == "CG":
                    sqe += 3.3
                if sq[i-1] == "GC":
                    sqe += 2.4
                if sq[i-1] == "UA":
                    sqe += 2.1
                if sq[i-1] == "GU":
                    sqe += 1.4
                if sq[i-1] == "UG":
                    sqe += 2.1
            if sq[i] == "GC":
                if sq[i-1] == "AU": 
                    sqe += 2.4
                if sq[i-1] == "CG":
                    sqe += 3.4
                if sq[i-1] == "GC":
                    sqe += 3.3
                if sq[i-1] == "UA":
                    sqe += 2.2
                if sq[i-1] == "GU":
                    sqe += 1.5
                if sq[i-1] == "UG":
                    sqe += 2.5
            if sq[i] == "UA":
                if sq[i-1] == "AU": 
                    sqe += 1.3
                if sq[i-1] == "CG":
                    sqe += 2.4
                if sq[i-1] == "GC":
                    sqe += 2.1
                if sq[i-1] == "UA":
                    sqe += 0.9
                if sq[i-1] == "GU":
                    sqe += 1.0
                if sq[i-1] == "UG":
                    sqe += 1.3
            if sq[i] == "GU":
                if sq[i-1] == "AU": 
                    sqe += 1.3
                if sq[i-1] == "CG":
                    sqe += 2.5
                if sq[i-1] == "GC":
                    sqe += 2.1
                if sq[i-1] == "UA":
                    sqe += 1.4
                if sq[i-1] == "GU":
                    sqe += 0.5
                if sq[i-1] == "UG":
                    sqe += -1.3
            if sq[i] == "UG":
                if sq[i-1] == "AU": 
                    sqe += 1.0
                if sq[i-1] == "CG":
                    sqe += 1.5
                if sq[i-1] == "GC":
                    sqe += 1.4
                if sq[i-1] == "UA":
                    sqe += 0.6
                if sq[i-1] == "GU":
                    sqe += -0.3
                if sq[i-1] == "UG":
                    sqe += 0.5
    return sqe

def actual_SQs(seq_ss, seq_ps):    # seq_ss: secondary structure, seq_ps: primary structure (sequence)
    
    with open(seq_ss) as file:
        ss_lines = file.readlines()
    
    with open(seq_ps) as file:
        ps_lines = file.readlines()
    
    rna = ps_lines[1]
    
    SQs_actual = []

    sip = False                       # SQ in progress?
    sl = 0                            # SQ length
    sp = []                           # SQ pairs
    last_line = [0, 0, 0, 0, 0, 0]    # initiate last line

    for i in range(0, len(ss_lines)):
        line = ss_lines[i].strip().split()
        
        if (int(line[4]) != 0 and sip == False):
            sip = True
            temp = [int(line[0]), int(line[4])]
            if (rna[i] == ('G' or 'g') and rna[int(line[4])-1] == ('C' or 'c')):
                sp.append("GC")
            elif (rna[i] == ('C' or 'c') and rna[int(line[4])-1] == ('G' or 'g')):
                sp.append("CG")
            elif (rna[i] == ('G' or 'g') and rna[int(line[4])-1] == ('U' or 'u')):
                sp.append("GU")
            elif (rna[i] == ('U' or 'u') and rna[int(line[4])-1] == ('G' or 'g')):
                sp.append("UG")
            elif (rna[i] == ('A' or 'a') and rna[int(line[4])-1] == ('U' or 'u')):
                sp.append("AU")
            elif (rna[i] == ('U' or 'u') and rna[int(line[4])-1] == ('A' or 'a')):
                sp.append("UA")
            else: 
                sp.append("noncanonical")
            sl += 1
            
        elif (int(line[4]) != 0 and sip == True and (int(last_line[4])-int(line[4]) == 1)):
            if (rna[i] == ('G' or 'g') and rna[int(line[4])-1] == ('C' or 'c')):
                sp.append("GC")
            elif (rna[i] == ('C' or 'c') and rna[int(line[4])-1] == ('G' or 'g')):
                sp.append("CG")
            elif (rna[i] == ('G' or 'g') and rna[int(line[4])-1] == ('U' or 'u')):
                sp.append("GU")
            elif (rna[i] == ('U' or 'u') and rna[int(line[4])-1] == ('G' or 'g')):
                sp.append("UG")
            elif (rna[i] == ('A' or 'a') and rna[int(line[4])-1] == ('U' or 'u')):
                sp.append("AU")
            elif (rna[i] == ('U' or 'u') and rna[int(line[4])-1] == ('A' or 'a')):
                sp.append("UA")
            else: 
                sp.append("noncanonical")
            sl += 1
            if sl == 2:
                if temp[1] > temp[0]:
                    temp.append(SQ_energy(sp[-2:]))
                    SQs_actual.append(temp)
                temp = [int(line[0]), int(line[4])]
                sl = 1
            
        elif (int(line[4]) == 0 and sip == True):
            sip = False
            if sl == 2:
                if temp[1] > temp[0]:
                    temp.append(SQ_energy(sp[-2:]))
                    SQs_actual.append(temp)
            elif sl == 1:
                sp.pop()
            sl = 0
            
        elif ((int(last_line[4])-int(line[4]) != 1) and int(last_line[4]) != 0  and sip == True):

            if sl == 2:
                if temp[1] > temp[0]:
                    temp.append(SQ_energy(sp[-2:]))
                    SQs_actual.append(temp)
            elif sl == 1:
                sp.pop()
            temp = [int(line[0]), int(line[4])]
            sl = 0
            
            if (rna[i] == ('G' or 'g') and rna[int(line[4])-1] == ('C' or 'c')):
                sp.append("GC")
            elif (rna[i] == ('C' or 'c') and rna[int(line[4])-1] == ('G' or 'g')):
                sp.append("CG")
            elif (rna[i] == ('G' or 'g') and rna[int(line[4])-1] == ('U' or 'u')):
                sp.append("GU")
            elif (rna[i] == ('U' or 'u') and rna[int(line[4])-1] == ('G' or 'g')):
                sp.append("UG")
            elif (rna[i] == ('A' or 'a') and rna[int(line[4])-1] == ('U' or 'u')):
                sp.append("AU")
            elif (rna[i] == ('U' or 'u') and rna[int(line[4])-1] == ('A' or 'a')):
                sp.append("UA")
            else: 
                sp.append("noncanonical")
            sl += 1
        
        last_line = line
        
    return SQs_actual

def potential_SQs(seq_ps):
    
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
    
    SQs_potential = []

    for row in range(0, len(matrix)):
        for col in range (row, len(matrix)):
            if row != col:
                if matrix[row][col] != 0:
                    SQp = []                 # stacked quartet pairs
                    temp_row = row
                    temp_col = col
                    SQ = [row+1, col+1, 0] # [SQ start, SQ end, SQ energy]
                    length = 0
                    while (matrix[temp_row][temp_col] != 0) and (temp_row != temp_col):
                        base1 = rna[temp_row]
                        base2 = rna[temp_col]
                        if (base1 == ('G' or 'g') and base2 == ('C' or 'c')):
                            SQp.append("GC")
                        if (base1 == ('C' or 'c') and base2 == ('G' or 'g')):
                            SQp.append("CG")
                        if (base1 == ('G' or 'g') and base2 == ('U' or 'u')):
                            SQp.append("GU")
                        if (base1 == ('U' or 'u') and base2 == ('G' or 'g')):
                            SQp.append("UG")
                        if (base1 == ('A' or 'a') and base2 == ('U' or 'u')):
                            SQp.append("AU")
                        if (base1 == ('U' or 'u') and base2 == ('A' or 'a')):
                            SQp.append("UA")
                        length += 1
                        temp_row += 1
                        temp_col -= 1
                        if length == 2 and col-row-2*length >= 3:
                            SQ[2] = SQ_energy(SQp)
                            SQs_potential.append(SQ.copy())
                            break
    
    return [SQs_potential, rna, len(rna)]

def potential_couplings(SQs_potential):
    
    overlaps_potential = []
    nestings_potential = []
    pseudoknots_potential = []
    
    overlap_penalty = 1e6

    for i in range(len(SQs_potential)):
        for j in range(i+1, len(SQs_potential)):
    
            SQ1 = SQs_potential[i]
            SQ2 = SQs_potential[j]
            
            i_a = SQ1[0]
            j_a = SQ1[1]
            i_b = SQ2[0]
            j_b = SQ2[1]
    
            overlaps = [i, j, 0]    
            nestings = [i, j, 0]
            pseudoknots = [i, j, 0]
    
            SQ1_cspan1 = set(range(SQ1[1]-2+1, SQ1[1]+1))
            SQ2_cspan1 = set(range(SQ2[1]-2+1, SQ2[1]+1))
            
            SQ1_cspan2 = set(range(SQ1[0], SQ1[0]+2))
            SQ2_cspan2 = set(range(SQ2[0], SQ2[0]+2))
    
            if (len(SQ1_cspan1 & SQ2_cspan1) != 0) or (len(SQ1_cspan2 & SQ2_cspan2) != 0)  or (len(SQ1_cspan1 & SQ2_cspan2) != 0) or (len(SQ1_cspan2 & SQ2_cspan1) != 0):
                
                if (SQ1[0] == SQ2[0]+1 and SQ1[1] == SQ2[1]-1) or (SQ2[0] == SQ1[0]+1 and SQ2[1] == SQ1[1]-1):
                    nestings[2] = 1
                else:
                    overlaps[2] = overlap_penalty
            elif (i_a < i_b and i_b < j_a and j_a < j_b) or (i_b < i_a and i_a < j_b and j_b < j_a):
                pseudoknots[2] = 1
            
            overlaps_potential.append(overlaps)
            nestings_potential.append(nestings)
            pseudoknots_potential.append(pseudoknots)
            
    return (overlaps_potential, nestings_potential, pseudoknots_potential)

def model(SQs_potential, overlaps_potential, nestings_potential, pseudoknots_potential, mplus, mminus):
    
    L = {}
    Q = {}
    k = 0
    
    for i in range(0, len(SQs_potential)):
        L[str(i)] = -SQs_potential[i][2]
        for j in range(i+1, len(SQs_potential)):
            Q[(str(i), str(j))] = overlaps_potential[k][2] - mplus*nestings_potential[k][2] - mminus*pseudoknots_potential[k][2]
            k += 1
            
    return L, Q

def energy(SQs_actual, mplus, mminus):
    
    k = 0
    couplings_actual = potential_couplings(SQs_actual)
    nestings_actual = couplings_actual[1]
    pseudoknots_actual = couplings_actual[2]
    cost = 0
        
    for i in range(0, len(SQs_actual)):
        cost -= SQs_actual[i][2]
        for j in range(i+1, len(SQs_actual)):
            cost -= mplus*nestings_actual[k][2] + mminus*pseudoknots_actual[k][2]
            k += 1
    
    return cost

def evaluation_1(SQs_actual, SQs_potential):
    
    bp_actual = []
    bp_predicted = []

    for i in range(0, len(SQs_actual)):
        for j in range(0, 2):
            bp_actual.append((SQs_actual[i][0]+j, SQs_actual[i][1]-j))
        
    for i in range(0, len(SQs_potential)):
        for j in range(0, 2):
            bp_predicted.append((SQs_potential[i][0]+j, SQs_potential[i][1]-j))
            
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

def evaluation_2(SQs_actual, SQs_predicted):
    
    pb_actual = []
    pb_predicted = []

    for i in range(0, len(SQs_actual)):
        for j in range(0, 2):
            pb_actual.append(SQs_actual[i][0]+j)
            pb_actual.append(SQs_actual[i][1]-j)
        
    for i in range(0, len(SQs_predicted)):
        for j in range(0, 2):
            pb_predicted.append(SQs_predicted[i][0]+j)
            pb_predicted.append(SQs_predicted[i][1]-j)
            
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

def MCC(SQs_actual, SQs_predicted, SQs_potential):
    
    bp_actual = []
    bp_predicted = []
    bp_potential = []

    for i in range(0, len(SQs_actual)):
        for j in range(0, 2):
            bp_actual.append((SQs_actual[i][0]+j, SQs_actual[i][1]-j))
        
    for i in range(0, len(SQs_predicted)):
        for j in range(0, 2):
            bp_predicted.append((SQs_predicted[i][0]+j, SQs_predicted[i][1]-j))
            
    for i in range(0, len(SQs_potential)):
        for j in range(0, 2):
            bp_potential.append((SQs_potential[i][0]+j, SQs_potential[i][1]-j))
            
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
    a_SQs = {}
    p_SQs = {}
    p_couplings = {}
    prediction = {}
    for bprna in full_bprna:
        bprna_id = bprna.split("/")[8]
        print(bprna_id)
        fasta_file = bprna + ".fasta.txt"
        ct_file = bprna + ".ct.txt" 
        p_SQs[bprna_id] = potential_SQs(fasta_file)
        p_couplings[bprna_id] = potential_couplings(p_SQs[bprna_id][0])
        a_SQs[bprna_id] = actual_SQs(ct_file, fasta_file)
    print("Finished pre-processing")
    
    def cost_function(hyper_params):
        print(hyper_params)
        if str(hyper_params) in past_params:
            print("Previous value")
            return past_params[str(hyper_params)]
        else:
            mplus = hyper_params[0]
            mminus = hyper_params[1]
            problems = {}
            a_energies = {}
            for bprna in full_bprna:
                print(bprna)
                bprna_id = bprna.split("/")[8]
                a_energies[bprna_id] = energy(a_SQs[bprna_id], mplus, mminus)
                md = model(p_SQs[bprna_id][0], p_couplings[bprna_id][0], p_couplings[bprna_id][1], p_couplings[bprna_id][2], mplus, mminus)
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
        
                SQs_found = []

                for j in range(0, len(results_hybrid)):
                    if results_hybrid[str(j)] == 1:
                        SQs_found.append(p_SQs[key][0][j])
                    
                prediction[key] = ([SQs_found, predicted_energy])
                bp_metric.append(MCC(a_SQs[bprna_id], SQs_found, p_SQs[key][0]))
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

    print("m2", len(training_data))    
    print("Finished data set creation")
    initial_point = [1, 1]
    model_params = [1, 1]
    train_history = []
    optimizer = SPSA(maxiter=30, callback=lambda n, p, v, ss, sa: spsa_optimizer_callback(n, p, v, ss, sa, train_history))
    model_values, loss, nfev = optimize_params(optimizer, model_params, initial_point, training_data, solver, num_repeats, solver_limit, num_reads, seed)
    print("Final values:",model_values)
    print('Save results')
    save_job_result({"final_values": str(model_values), "loss":str(loss), "hyperparams": str(hyperparams)})

if __name__ == "__main__":
    main()