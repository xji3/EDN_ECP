from IGCexpansion.CodonGeneconv import ReCodonGeneconv
from IGCexpansion.IndCodonGeneconv import IndCodonGeneconv
from IGCexpansion.HMMJSGeneconv import HMMJSGeneconv
import argparse
import numpy as np

def main(args):
    paralog = [args.paralog1, args.paralog2]
    Force = None
    alignment_file = './MafftAlignment/' + '_'.join(paralog) + '/' + '_'.join(paralog) + '_input.fasta'
    newicktree = './YeastTree.newick'
    if args.force:
        if args.model == 'MG94':
            Force = {5:0.0}
        elif args.model == 'HKY':
            Force = {4:0.0}
    else:
        Force = None
    test = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = args.model, Force = Force, clock = args.clock)
    test.get_mle(True, True, 0, 'BFGS')
    test.get_individual_summary(summary_path = './summary/')
    #test.get_SitewisePosteriorSummary(summary_path = './Summary/')
    if Force == None:
        test.get_sitewise_loglikelihood_summary('./summary/' + '_'.join(paralog) + '_MG94_nonclock_sw_lnL.txt')
    else:
        test.get_sitewise_loglikelihood_summary('./summary/Force_' + '_'.join(paralog) + '_MG94_nonclock_sw_lnL.txt')


if __name__ == '__main__':
    pair = ["EDN", "ECP"]

    paralog = pair
    Force = None
    #Force = {5:0.0}
    alignment_file = './EDN_ECP_Cleaned.fasta'
    newicktree = './input_tree.newick'

##    test = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = 'MG94', Force = Force, clock = None)
##    test.get_mle(True, True, 0, 'BFGS')
##    test.get_sitewise_loglikelihood_summary('./Summary/' + '_'.join(paralog) + '_MG94_nonclock_sw_lnL.txt')
##
##    Force = {5:0.0}
##    test_force = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = 'MG94', Force = Force, clock = None)
##    test_force.update_by_x(test.x)
##    test_force._loglikelihood2()
##    test_force.get_sitewise_loglikelihood_summary('./Summary/Force_' + '_'.join(paralog) + '_MG94_nonclock_sw_lnL.txt')
##
##    Force = None
##    test = IndCodonGeneconv( newicktree, alignment_file, paralog, Model = 'MG94', Force = Force, clock = None)
##    test.get_mle(True, True, 0, 'BFGS')
##    test.get_sitewise_loglikelihood_summary('./Summary/' + '_'.join(paralog) + '_Ind_MG94_nonclock_sw_lnL.txt')
##
##    Force = {5:0.0}
##    test_force = IndCodonGeneconv( newicktree, alignment_file, paralog, Model = 'MG94', Force = Force, clock = None)
##    test_force.update_by_x(test.x)
##    test_force._loglikelihood2()
##    test_force.get_sitewise_loglikelihood_summary('./Summary/Force_' + '_'.join(paralog) + '_Ind_MG94_nonclock_sw_lnL.txt')


    save_file = './save/HMMJS_' + '_'.join(paralog) + '_Ind_MG94_nonclock_save.txt'
    IGC_save_name = './save/HMMJS_' + '_'.join(paralog) + '_Ind_MG94_nonclock_IGC_save.txt'
    Force_save_name = './save/HMMJS_' + '_'.join(paralog) + '_Ind_MG94_nonclock_Force_save.txt'

    summary_path = './Summary/'
    IGC_sitewise_lnL_file = './Summary/' + '_'.join(paralog) + '_Ind_MG94_nonclock_sw_lnL.txt'
    Force_sitewise_lnL_file = './Summary/Force_' + '_'.join(paralog) + '_Ind_MG94_nonclock_sw_lnL.txt'
    save_path = './save/'
    seq_index_file = './' + '_'.join(paralog) + '_seq_index.txt'
    state_list = ['No IGC event (Si = 0)','At least one IGC event (Si > 0)']

    x_3 = np.array([-0.69727878,
                  -0.53710801,
                  -0.72400474,
                  0.72385788,
                  -0.18983735,
                  -2.5199719,
                  -2.01452935,
                  -1.0337633,
                  -3.29369029,
                  -1.75318807,
                  -3.25869777,
                  -2.27341043,
                  -4.20160402,
                  -4.110472,
                  -0.5])
    
    test = HMMJSGeneconv(save_file, newicktree, alignment_file, paralog, summary_path, x_3, save_path, IGC_save_name, Force_save_name,
                         IGC_sitewise_lnL_file, Force_sitewise_lnL_file,
                         state_list, seq_index_file)
    test.get_mle(display = True, two_step = False, OneDimension = True)
