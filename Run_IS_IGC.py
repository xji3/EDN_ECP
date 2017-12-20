from IGCexpansion.CodonGeneconv import ReCodonGeneconv
from IGCexpansion.IndCodonGeneconv import IndCodonGeneconv
import argparse

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
    test = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = 'MG94', Force = Force, clock = None)
    test.get_mle(True, True, 0, 'BFGS')
    test.get_sitewise_loglikelihood_summary('./Summary/' + '_'.join(paralog) + '_MG94_nonclock_sw_lnL.txt')

    Force = {5:0.0}
    test_force = ReCodonGeneconv( newicktree, alignment_file, paralog, Model = 'MG94', Force = Force, clock = None)
    test_force.update_by_x(test.x)
    test_force._loglikelihood2()
    test_force.get_sitewise_loglikelihood_summary('./Summary/Force_' + '_'.join(paralog) + '_MG94_nonclock_sw_lnL.txt')

    Force = None
    test = IndCodonGeneconv( newicktree, alignment_file, paralog, Model = 'MG94', Force = Force, clock = None)
    test.get_mle(True, True, 0, 'BFGS')
    test.get_sitewise_loglikelihood_summary('./Summary/' + '_'.join(paralog) + '_Ind_MG94_nonclock_sw_lnL.txt')

    Force = {5:0.0}
    test_force = IndCodonGeneconv( newicktree, alignment_file, paralog, Model = 'MG94', Force = Force, clock = None)
    test_force.update_by_x(test.x)
    test_force._loglikelihood2()
    test_force.get_sitewise_loglikelihood_summary('./Summary/Force_' + '_'.join(paralog) + '_Ind_MG94_nonclock_sw_lnL.txt')

