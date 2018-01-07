from IGCexpansion.PSJSGeneconv import *
from IGCexpansion.JSGeneconv import *
import argparse
from collections import namedtuple
import numpy as np

def main(args):
    paralog = [args.paralog1, args.paralog2]

    gene_to_orlg_file = './' + '_'.join(paralog) +'_GeneToOrlg.txt'
    alignment_file = './' + '_'.join(paralog) +'_Cleaned_input.fasta'

    tree_newick = './input_tree.newick'
    DupLosList = './EDN_ECP_DupLost.txt'
    terminal_node_list = ['Tamarin', 'Macaque', 'Orangutan', 'Gorilla', 'Chimpanzee']
    node_to_pos = {'D1':0}
    pm_model = 'HKY'
    
    IGC_pm = 'One rate'
    seq_index_file = './' + '_'.join(paralog) +'_seq_index.txt'

    if args.rate_variation:
        save_file = './save/JS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_rv_nonclock_save.txt'
        log_file = './log/JS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_rv_nonclock_log.txt'
        summary_file = './summary/JS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_rv_nonclock_summary.txt'
        x_js = np.log([ 0.5, 0.5, 0.5,  4.35588244, 0.5, 5.0, 0.3])
    else:
        save_file = './save/JS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_nonclock_save.txt'
        log_file = './log/JS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_nonclock_log.txt'
        summary_file = './summary/JS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_nonclock_summary.txt'
        x_js = np.log([ 0.5, 0.5, 0.5,  4.35588244,   0.3])



    test_JS = JSGeneconv(alignment_file, gene_to_orlg_file, args.cdna, tree_newick, DupLosList, x_js, pm_model, IGC_pm,
                         args.rate_variation, node_to_pos, terminal_node_list, save_file, log_file = log_file)
    test_JS.get_mle()
    #test_JS.get_expectedNumGeneconv()
    test_JS.get_individual_summary(summary_file)



if __name__ == '__main__':
##    parser = argparse.ArgumentParser()
##    parser.add_argument('--paralog1', required = True, help = 'Name of the 1st paralog')
##    parser.add_argument('--paralog2', required = True, help = 'Name of the 2nd paralog')
##    parser.add_argument('--L', type = float, dest = 'tract_length', default = 30.0, help = 'Initial guess tract length')
##    parser.add_argument('--D', type = int, dest = 'dim', default = 0, help = 'Dimension used in search with default value 0')
##    parser.add_argument('--heterogeneity', dest = 'rate_variation', action = 'store_true', help = 'rate heterogeneity control')
##    parser.add_argument('--homogeneity', dest = 'rate_variation', action = 'store_false', help = 'rate heterogeneity control')
##    parser.add_argument('--coding', dest = 'cdna', action = 'store_true', help = 'coding sequence control')
##    parser.add_argument('--noncoding', dest = 'cdna', action = 'store_false', help = 'coding sequence control')
##    parser.add_argument('--samecodon', dest = 'allow_same_codon', action = 'store_true', help = 'whether allow pair sites from same codon')
##    parser.add_argument('--no-samecodon', dest = 'allow_same_codon', action = 'store_false', help = 'whether allow pair sites from same codon')
##    
##    main(parser.parse_args())

  

    MyStruct = namedtuple('MyStruct', 'paralog1 paralog2 tract_length dim cdna rate_variation allow_same_codon')
    args = MyStruct(paralog1 = 'EDN', paralog2 = 'ECP', tract_length = 30.0, dim = 1, cdna = True, rate_variation = True, allow_same_codon = True)

    paralog = [args.paralog1, args.paralog2]

    gene_to_orlg_file = './' + '_'.join(paralog) +'_GeneToOrlg.txt'
    alignment_file = './' + '_'.join(paralog) +'_Cleaned_input.fasta'

    tree_newick = './input_tree.newick'
    DupLosList = './EDN_ECP_DupLost.txt'
    terminal_node_list = ['Tamarin', 'Macaque', 'Orangutan', 'Gorilla', 'Chimpanzee']
    node_to_pos = {'D1':0}
    pm_model = 'HKY'
    
    IGC_pm = 'One rate'
    seq_index_file = './' + '_'.join(paralog) +'_seq_index.txt'

    force = None
    if args.rate_variation:
        save_file = './save/JS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_rv_nonclock_save.txt'
        log_file = './log/JS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_rv_nonclock_log.txt'
        summary_file = './summary/JS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_rv_nonclock_summary.txt'
        gradient_file = './summary/JS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_rv_nonclock_gradient.txt'
        hessian_file = './summary/JS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_rv_nonclock_hessian.txt'
        godambe_file = './summary/JS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_rv_nonclock_godambe.txt'
        x_js = np.log([ 0.5, 0.5, 0.5,  4.35588244, 0.5, 5.0, 0.3])
    else:
        save_file = './save/JS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_nonclock_save.txt'
        log_file = './log/JS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_nonclock_log.txt'
        summary_file = './summary/JS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_nonclock_summary.txt'
        gradient_file = './summary/JS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_nonclock_gradient.txt'
        hessian_file = './summary/JS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_nonclock_hessian.txt'
        godambe_file = './summary/JS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_nonclock_godambe.txt'
        x_js = np.log([ 0.5, 0.5, 0.5,  4.35588244,   0.3])



    test_JS = JSGeneconv(alignment_file, gene_to_orlg_file, args.cdna, tree_newick, DupLosList, x_js, pm_model, IGC_pm,
                         args.rate_variation, node_to_pos, terminal_node_list, save_file, force = force)
    test_JS.get_mle()
    test_JS.get_expectedNumGeneconv()
    test_JS.get_expectedMutationNum()
    test_JS.get_individual_summary(summary_file)
    godambe = test_JS.get_Godambe_matrix(test_JS.x, gradient_file, hessian_file, 1e-6)
    np.savetxt(open(godambe_file, 'w+'), np.array(godambe))



    force = {6:0.0}
    if args.rate_variation:
        save_file = './save/Force_JS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_rv_nonclock_save.txt'
        log_file = './log/Force_JS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_rv_nonclock_log.txt'
        summary_file = './summary/Force_JS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_rv_nonclock_summary.txt'
        x_js = np.log([ 0.5, 0.5, 0.5,  4.35588244, 0.5, 5.0, 0.3])
    else:
        save_file = './save/Force_JS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_nonclock_save.txt'
        log_file = './log/Force_JS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_nonclock_log.txt'
        summary_file = './summary/Force_JS_HKY_'+ '_'.join(paralog) + '_' + IGC_pm.replace(' ', '_') + '_nonclock_summary.txt'
        x_js = np.log([ 0.5, 0.5, 0.5,  4.35588244,   0.3])



    test_JS_Force = JSGeneconv(alignment_file, gene_to_orlg_file, args.cdna, tree_newick, DupLosList, x_js, pm_model, IGC_pm,
                         args.rate_variation, node_to_pos, terminal_node_list, save_file, force = force)
    test_JS_Force.get_mle()
    #test_JS.get_expectedNumGeneconv()
    test_JS_Force.get_individual_summary(summary_file)
        
