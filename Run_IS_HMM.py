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

    save_name = './save/EDN_ECP_Ind_MG94_IGC_save.txt'
    summary_file = './summary/EDN_ECP_Ind_MG94_HMM_summary.txt'
    summary_file_1D = './summary/EDN_ECP_Ind_MG94_HMM_1D_summary.txt'
    Ind_MG94_IGC = IndCodonGeneconv( newicktree, alignment_file, paralog, Model = 'MG94', Force = Force, clock = None, save_path = './save/',
                             save_name = save_name)
    Ind_MG94_IGC.get_mle(True, True, 0, 'BFGS')
    x = np.concatenate((Ind_MG94_IGC.x, [0.0]))

    IGC_sitewise_lnL_file = summary_path + '_'.join(paralog) + '_Ind_MG94_nonclock_sw_lnL.txt'
    NOIGC_sitewise_lnL_file = summary_path + 'NOIGC_' + '_'.join(paralog) + '_Ind_MG94_nonclock_sw_lnL.txt'
        
    Ind_MG94_IGC.get_sitewise_loglikelihood_summary(IGC_sitewise_lnL_file, False)
    Ind_MG94_IGC.get_sitewise_loglikelihood_summary(NOIGC_sitewise_lnL_file, True)

    
    test = HMMJSGeneconv(save_file, newicktree, alignment_file, paralog, summary_path, x, save_path, IGC_sitewise_lnL_file, NOIGC_sitewise_lnL_file,                         state_list, seq_index_file)
    

    log_p_list = np.log(3.0/np.array(range(3, 1001)))
    plot_file = './plot/HMM_' + '_'.join(paralog) + '_lnL_1D_surface.txt'
    test.plot_tract_p(log_p_list, plot_file)

    test.get_mle(display = True, two_step = True, One_Dimension = True)
    test.get_summary(summary_file_1D)

    grad, hess = test.get_Hessian(True)
    hessian_summary_file = './summary/EDN_ECP_Ind_MG94_HMM_Hessian.txt'
    np.savetxt(open(hessian_summary_file, 'w+'), np.matrix([[grad], [hess]]), delimiter = ' ', footer = 'Grad Hessian')
    
    lnL_array = test.hmmtract.get_posterior()
    posterior_lnL_file = './summary/EDN_ECP_Ind_MG94_HMM_Posterior_lnL.txt'
    np.savetxt(open(posterior_lnL_file, 'w+'), np.matrix(lnL_array).T, delimiter = ' ', footer = 'Si=0 Si=1')

    Viterbi_array = test.hmmtract.Viterbi()[1]
    Viterbi_array_file = './summary/EDN_ECP_Ind_MG94_HMM_Viterbi_array.txt'
    np.savetxt(open(Viterbi_array_file, 'w+'), np.matrix(Viterbi_array).T, delimiter = ' ', footer = 'Viterbi')
