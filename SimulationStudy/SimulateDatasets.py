from IGCexpansion.CodonSimulator import CodonSimulator
from IGCexpansion.IndCodonGeneconv import IndCodonGeneconv
from IGCexpansion.JSGeneconv import JSGeneconv
from IGCexpansion.Simulator import Simulator
import argparse, os
import numpy as np

if __name__ == '__main__':

##### Simulate according to MG94 Model
########## Check MLE
##    paralog = ['EDN', 'ECP']
##    gene_to_orlg_file = '../EDN_ECP_GeneToOrlg.txt'
##    alignment_file = '../EDN_ECP_Cleaned.fasta'
##    newicktree = '../input_tree.newick'
##    DupLosList = '../EDN_ECP_DupLost.txt'
##    Force = None
##    terminal_node_list = ['Tamarin', 'Macaque', 'Orangutan', 'Gorilla', 'Chimpanzee']
##    node_to_pos = {'D1':0}
##    seq_index_file = '../' + '_'.join(paralog) +'_seq_index.txt'
##
##    Ind_MG94_IGC = IndCodonGeneconv( newicktree, alignment_file, paralog, Model = 'MG94', Force = Force, clock = None, save_path = './save/',
##                             save_name = './save/EDN_ECP_Ind_MG94_IGC_save.txt')
##    #Ind_MG94_IGC.get_mle(True, True, 0, 'BFGS')

######## Now simulate datasets
##    display = False
##    mean_tract_list = [3.0, 10.0, 50.0, 100.0, 200.0, 300.0, 400.0, 500.0]
##    sim_num = 1
##    seed_number_start = 27606
##    
##    for geo in mean_tract_list[:]:
##        for sim_num in range(1, 101):
##
##            seed_number = seed_number_start + sim_num
##
##            if not os.path.isdir('./MG94_Tract_' + str(geo)):
##                os.mkdir('./MG94_Tract_' + str(geo))
##            if not os.path.isdir('./MG94_Tract_' + str(geo) + '/sim_' + str(sim_num)):
##                os.mkdir('./MG94_Tract_' + str(geo) + '/sim_' + str(sim_num))
##            seq_file = './MG94_Tract_' + str(geo) + '/sim_' + str(sim_num) + '/EDN_ECP_sim_' + str(sim_num) + '.fasta'
##            IGC_log_file = './MG94_Tract_' + str(geo) + '/sim_' + str(sim_num) + '/EDN_ECP_sim_' + str(sim_num) + '_IGC.log'
##            PM_log_file = './MG94_Tract_' + str(geo) + '/sim_' + str(sim_num) + '/EDN_ECP_sim_' + str(sim_num) + '_PM.log'
##
##            pm_model_name = 'MG94'
##            #x_pm = np.log([0.4, 0.5, 0.2, 9.2, 1.0])
##            x_pm = Ind_MG94_IGC.x_process[:-1]
##            rate_variation = False
##
##            x_IGC = [Ind_MG94_IGC.tau * 3.0 / geo, 3.0 / geo]
##            init_pm = 'One rate'
##            tract_pm = 'One rate'
##            pm_IGC = [init_pm, tract_pm]
##
##        ##            x_rates = [-4.170654939766711422e+00,
##        ##                       -5.674236262981605883e+00,
##        ##                       -4.140979602575983520e+00,
##        ##                       -4.344239699023852097e+00,
##        ##                       -6.496123290482403334e+00,
##        ##                       -6.063647134296714647e+00,
##        ##                       -6.043806966727234276e+00,
##        ##                       -5.111657692573940537e+00,
##        ##                       -6.404488905061815451e+00,
##        ##                       -5.467996717925044159e+00,
##        ##                       -5.460686727891754799e+00,
##        ##                       -6.459940982759793116e+00]
##            
##            x_rates = Ind_MG94_IGC.x_rates
##            
##            test = CodonSimulator(pm_model_name, x_pm, rate_variation,
##                             x_IGC, pm_IGC, newicktree, DupLosList, x_rates,
##                             terminal_node_list, node_to_pos, gene_to_orlg_file, seq_file, IGC_log_file, PM_log_file, seed_number, seq_index_file)
##
##            self = test
##            print test
##            
##            test.sim_root()
##            #edge = ('N0', 'kluyveri')
##            #test.sim_one_branch(edge, True)
##
##            test.sim(display = display)
##            test.output_seq()
##
##            test.seq_file = './MG94_Tract_' + str(geo) + '/sim_' + str(sim_num) + '/EDN_ECP_sim_' + str(sim_num) + '_newformat.fasta'
##            test.output_seq(True)


# Simulate according to HKY+rv Model
###### Check MLE
    paralog = ['EDN', 'ECP']
    gene_to_orlg_file = '../EDN_ECP_GeneToOrlg.txt'
    newicktree = '../input_tree.newick'
    DupLosList = '../EDN_ECP_DupLost.txt'
    Force = None
    terminal_node_list = ['Tamarin', 'Macaque', 'Orangutan', 'Gorilla', 'Chimpanzee']
    node_to_pos = {'D1':0}
    seq_index_file = '../' + '_'.join(paralog) +'_seq_index.txt'
    alignment_file = '../EDN_ECP_Cleaned_input.fasta'

    pm_model = 'HKY'
    IGC_pm = 'One rate'
    rate_variation = True
    save_file = './save/EDN_ECP_HKY_rv_IGC_nonclock_save.txt'

    x_js = np.log([ 0.5, 0.5, 0.5,  4.35588244, 0.5, 5.0, 0.3])
    test_JS = JSGeneconv(alignment_file, gene_to_orlg_file, True, newicktree, DupLosList, x_js, pm_model, IGC_pm,
                         rate_variation, node_to_pos, terminal_node_list, save_file)

    #test_JS.get_mle()

############ Now simulate datasets
##    display = False
##    mean_tract_list = [3.0, 10.0, 50.0, 100.0, 200.0, 300.0, 400.0, 500.0]
##    sim_num = 1
##    seed_number_start = 27606
##    
##    for geo in mean_tract_list[:]:
##        for sim_num in range(1, 101):
####    geo = 10.0
####    sim_num = 1
##
##            seed_number = seed_number_start + sim_num
##
##            if not os.path.isdir('./Tract_' + str(geo) + '_HKY'):
##                os.mkdir('./Tract_' + str(geo) + '_HKY')
##            if not os.path.isdir('./Tract_' + str(geo) + '_HKY' + '/sim_' + str(sim_num)):
##                os.mkdir('./Tract_' + str(geo) + '_HKY/sim_' + str(sim_num))
##            seq_file = './Tract_' + str(geo) + '_HKY/sim_' + str(sim_num) + '/EDN_ECP_sim_' + str(sim_num) + '.fasta'
##            IGC_log_file = './Tract_' + str(geo) + '_HKY/sim_' + str(sim_num) + '/EDN_ECP_sim_' + str(sim_num) + '_IGC.log'
##            PM_log_file = './Tract_' + str(geo) + '_HKY/sim_' + str(sim_num) + '/EDN_ECP_sim_' + str(sim_num) + '_PM.log'
##
##            pm_model_name = 'HKY'
##            #x_pm = np.log([0.4, 0.5, 0.2, 9.2, 1.0])
##            x_pm = test_JS.jsmodel.x_pm
##            rate_variation = True
##
##            x_IGC = [test_JS.jsmodel.IGCModel.parameters['Tau'] * 1.0 / geo, 1.0 / geo]
##            init_pm = 'One rate'
##            tract_pm = 'One rate'
##            pm_IGC = [init_pm, tract_pm]
##
##        ##            x_rates = [-4.170654939766711422e+00,
##        ##                       -5.674236262981605883e+00,
##        ##                       -4.140979602575983520e+00,
##        ##                       -4.344239699023852097e+00,
##        ##                       -6.496123290482403334e+00,
##        ##                       -6.063647134296714647e+00,
##        ##                       -6.043806966727234276e+00,
##        ##                       -5.111657692573940537e+00,
##        ##                       -6.404488905061815451e+00,
##        ##                       -5.467996717925044159e+00,
##        ##                       -5.460686727891754799e+00,
##        ##                       -6.459940982759793116e+00]
##            
##            if test_JS.root_by_dup:
##                x_rates = test_JS.x[-len(test_JS.tree.edge_list):]
##            else:
##                x_rates = test_JS.x[-(len(test_JS.tree.edge_list) - 1):]
##            
##            test = Simulator(pm_model_name, x_pm, rate_variation,
##                             x_IGC, pm_IGC, newicktree, DupLosList, x_rates,
##                             terminal_node_list, node_to_pos, gene_to_orlg_file, seq_file, IGC_log_file, PM_log_file, seed_number, seq_index_file)
##
##            self = test
##            print test
##            
##            test.sim_root()
##            #edge = ('N0', 'kluyveri')
##            #test.sim_one_branch(edge, True)
##
##            test.sim(display = display)
##            test.output_seq()
##
##            test.seq_file = './Tract_' + str(geo) + '_HKY/sim_' + str(sim_num) + '/EDN_ECP_sim_' + str(sim_num) + '_newformat.fasta'
##            test.output_seq(True)

############ Now simulate 1/10 Tau datasets
##    display = False
##    mean_tract_list = [3.0, 10.0, 50.0, 100.0, 200.0, 300.0, 400.0, 500.0]
##    sim_num = 1
##    seed_number_start = 27606
##    
##    for geo in mean_tract_list[:]:
##        for sim_num in range(1, 101):
####    geo = 10.0
####    sim_num = 1
##
##            seed_number = seed_number_start + sim_num
##
##            if not os.path.isdir('./TenthTau/Tract_' + str(geo) + '_HKY'):
##                os.mkdir('./TenthTau/Tract_' + str(geo) + '_HKY')
##            if not os.path.isdir('./TenthTau/Tract_' + str(geo) + '_HKY' + '/sim_' + str(sim_num)):
##                os.mkdir('./TenthTau/Tract_' + str(geo) + '_HKY/sim_' + str(sim_num))
##            seq_file = './TenthTau/Tract_' + str(geo) + '_HKY/sim_' + str(sim_num) + '/EDN_ECP_sim_' + str(sim_num) + '.fasta'
##            IGC_log_file = './TenthTau/Tract_' + str(geo) + '_HKY/sim_' + str(sim_num) + '/EDN_ECP_sim_' + str(sim_num) + '_IGC.log'
##            PM_log_file = './TenthTau/Tract_' + str(geo) + '_HKY/sim_' + str(sim_num) + '/EDN_ECP_sim_' + str(sim_num) + '_PM.log'
##
##            pm_model_name = 'HKY'
##            #x_pm = np.log([0.4, 0.5, 0.2, 9.2, 1.0])
##            x_pm = test_JS.jsmodel.x_pm
##            rate_variation = True
##
##            x_IGC = [test_JS.jsmodel.IGCModel.parameters['Tau'] * 0.1 / geo, 1.0 / geo]
##            init_pm = 'One rate'
##            tract_pm = 'One rate'
##            pm_IGC = [init_pm, tract_pm]
##
##        ##            x_rates = [-4.170654939766711422e+00,
##        ##                       -5.674236262981605883e+00,
##        ##                       -4.140979602575983520e+00,
##        ##                       -4.344239699023852097e+00,
##        ##                       -6.496123290482403334e+00,
##        ##                       -6.063647134296714647e+00,
##        ##                       -6.043806966727234276e+00,
##        ##                       -5.111657692573940537e+00,
##        ##                       -6.404488905061815451e+00,
##        ##                       -5.467996717925044159e+00,
##        ##                       -5.460686727891754799e+00,
##        ##                       -6.459940982759793116e+00]
##            
##            if test_JS.root_by_dup:
##                x_rates = test_JS.x[-len(test_JS.tree.edge_list):]
##            else:
##                x_rates = test_JS.x[-(len(test_JS.tree.edge_list) - 1):]
##            
##            test = Simulator(pm_model_name, x_pm, rate_variation,
##                             x_IGC, pm_IGC, newicktree, DupLosList, x_rates,
##                             terminal_node_list, node_to_pos, gene_to_orlg_file, seq_file, IGC_log_file, PM_log_file, seed_number, seq_index_file)
##
##            self = test
##            print test
##            
##            test.sim_root()
##            #edge = ('N0', 'kluyveri')
##            #test.sim_one_branch(edge, True)
##
##            test.sim(display = display)
##            test.output_seq()
##
##            test.seq_file = './TenthTau/Tract_' + str(geo) + '_HKY/sim_' + str(sim_num) + '/EDN_ECP_sim_' + str(sim_num) + '_newformat.fasta'
##            test.output_seq(True)

######## Now simulate 1/2 Tau datasets
##    display = False
##    mean_tract_list = [3.0, 10.0, 50.0, 100.0, 200.0, 300.0, 400.0, 500.0]
##    sim_num = 1
##    seed_number_start = 27606
##    
##    for geo in mean_tract_list[:]:
##        for sim_num in range(1, 101):
####    geo = 10.0
####    sim_num = 1
##
##            seed_number = seed_number_start + sim_num
##
##            if not os.path.isdir('./HalfTau/Tract_' + str(geo) + '_HKY'):
##                os.mkdir('./HalfTau/Tract_' + str(geo) + '_HKY')
##            if not os.path.isdir('./HalfTau/Tract_' + str(geo) + '_HKY' + '/sim_' + str(sim_num)):
##                os.mkdir('./HalfTau/Tract_' + str(geo) + '_HKY/sim_' + str(sim_num))
##            seq_file = './HalfTau/Tract_' + str(geo) + '_HKY/sim_' + str(sim_num) + '/EDN_ECP_sim_' + str(sim_num) + '.fasta'
##            IGC_log_file = './HalfTau/Tract_' + str(geo) + '_HKY/sim_' + str(sim_num) + '/EDN_ECP_sim_' + str(sim_num) + '_IGC.log'
##            PM_log_file = './HalfTau/Tract_' + str(geo) + '_HKY/sim_' + str(sim_num) + '/EDN_ECP_sim_' + str(sim_num) + '_PM.log'
##
##            pm_model_name = 'HKY'
##            #x_pm = np.log([0.4, 0.5, 0.2, 9.2, 1.0])
##            x_pm = test_JS.jsmodel.x_pm
##            rate_variation = True
##
##            x_IGC = [test_JS.jsmodel.IGCModel.parameters['Tau'] * 0.5 / geo, 1.0 / geo]
##            init_pm = 'One rate'
##            tract_pm = 'One rate'
##            pm_IGC = [init_pm, tract_pm]
##
##        ##            x_rates = [-4.170654939766711422e+00,
##        ##                       -5.674236262981605883e+00,
##        ##                       -4.140979602575983520e+00,
##        ##                       -4.344239699023852097e+00,
##        ##                       -6.496123290482403334e+00,
##        ##                       -6.063647134296714647e+00,
##        ##                       -6.043806966727234276e+00,
##        ##                       -5.111657692573940537e+00,
##        ##                       -6.404488905061815451e+00,
##        ##                       -5.467996717925044159e+00,
##        ##                       -5.460686727891754799e+00,
##        ##                       -6.459940982759793116e+00]
##            
##            if test_JS.root_by_dup:
##                x_rates = test_JS.x[-len(test_JS.tree.edge_list):]
##            else:
##                x_rates = test_JS.x[-(len(test_JS.tree.edge_list) - 1):]
##            
##            test = Simulator(pm_model_name, x_pm, rate_variation,
##                             x_IGC, pm_IGC, newicktree, DupLosList, x_rates,
##                             terminal_node_list, node_to_pos, gene_to_orlg_file, seq_file, IGC_log_file, PM_log_file, seed_number, seq_index_file)
##
##            self = test
##            print test
##            
##            test.sim_root()
##            #edge = ('N0', 'kluyveri')
##            #test.sim_one_branch(edge, True)
##
##            test.sim(display = display)
##            test.output_seq()
##
##            test.seq_file = './HalfTau/Tract_' + str(geo) + '_HKY/sim_' + str(sim_num) + '/EDN_ECP_sim_' + str(sim_num) + '_newformat.fasta'
##            test.output_seq(True)
