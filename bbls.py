#!/usr/bin/env python

'''NAME
        %(progname)s

VERSION
        %(version)s

AUTHOR
        Aldo Hernandez Corchado  <aldohernandezcorchado@gmail.com>

DESCRIPTION
        Determinates the bayesian branch length score (BBLS).
        Describes the level of conservation for motif sites.
        Described in: https://academic.oup.com/bioinformatics/article/25/2/167/219205

CATEGORY
        motif evaluation?

USAGE
        %(usage)s

ARGUMENTS
  GENERAL OPTIONS
    --version             show program's version number and exit
    -h, --help            show this help message and exit
    -v #, --verbosity=#   set verbosity to level #

'''
########################################
#                                      #
# IMPORTS
#                                      #
########################################

import copy
import os
import sys
import optparse
#from pydoc import pager
# BEGIN DEPENDENCIES
import ete3
from ete3 import Tree

# END DEPENDENCIES

########################################
#                                      #
#  COMMAND LINE OPTIONS
#                                      #
########################################
sys.path.insert(1, os.path.join(sys.path[0], 'lib'))
USAGE = '''%s -l bbls [-t inputfile] [-s inputfile] [-o outputfile] [-h | --help]'''

VERSION = '0'
HEADER = ''';
; bbls
x; Describes the level of conservation for motif sites.
; %(command)s
; version                          %(version)s
; date                             %(date)s
; running time                     %(runningTime)s
; sequences                        %(sequences)s
;'''

PROG_NAME = os.path.basename(sys.argv[0])

parser = optparse.OptionParser(usage=USAGE % PROG_NAME, add_help_option=0, version=VERSION)

parser.add_option("-h", "--help", action="store_true", dest="help")
parser.add_option("-l", "--length", action="store", dest="l", type="int")
parser.add_option("-s", "--scores", action="store", dest="score", default=sys.stdin)
parser.add_option("-t", "--tree", action="store", dest="tree", default=sys.stdin)
parser.add_option("-o", "--output", action="store", dest="output", default=sys.stdout)

(options, args) = parser.parse_args()

#print(parser.print_usage())

########################################
#                                      #
#  Classes
#                                      #
########################################

class ModTree(Tree):
    '''
    Description of variables i add
    '''
    prob, complement_prob, effective_length = None, None, None

    def to_string(self):
        '''
        docstring method
        '''
        print("Probability: "                    + str(self.prob) +\
              "\nTree complement's probability: "  + str(self.complement_prob) + \
              "\nEffective length: "               + str(self.effective_length) + '\n')

    def see_state(self):
        '''
        docstring method
        '''
        print('\n- Probabilities')
        print(self.get_ascii(attributes=['name', 'prob']))
        print('\n- Tree complement probabilities')
        print(self.get_ascii(attributes=['name', 'complement_prob']))
        print('\n- Effective length')
        print(self.get_ascii(attributes=['name', 'effective_length']))


########################################
#                                      #
#  MAIN
#                                      #
########################################

def main(args, options):
    """
    Description:
    Determines the Bayesian Branch Length Score (BBLS) of a motif site.
    Describes the conservation rate of a site.

    Parameters:
    Path to newick formated tree file.
    example:
    (((aa:1,bb:3)z:10,cc:5)y:100,dd:7)R:0;

    Path to dictionary of motif sites.
    example:
    aa 0.5
    bb 0.5
    cc 0.5
    dd 0.5

    Returns:
    BBLS and Branch Length Score.
    """

    def match(tree, motifs):
        """
        Description:
        Returns boolean value. Asserts that the leaves of the tree match the motif dictionary.

        Parameters:
        ete3 tree and motif dictionary

        Returns:
        Gets dictionary from the tree & compares.
        """

        leaves = set(tree.get_leaf_names())
        motifs_keys = set(motifs.keys())

        return leaves == motifs_keys

    def open_format_input(options):
        """
        Description:
        Loads the tree and motif scores files into usable formats.

        Parameters:
        options - >  options.score path to scores.
                     options.tree path to tree in newick format.

        Returns:
        my_tree -> moded tree object base in ete3 Tree.
        motifs  -> dictionary of motif scores.
        """

        # Reads motif file into a dictionary, raises error if score > 1.

        motifs = {}
        with open(options.score, 'r') as motif_list:
            for line in motif_list:

                words = line.split(" ")
                motif_specie = str(words[0])
                motif_score = float(words[1][:-1])

                if float(motif_score) > 1:
                    raise ValueError('Motif score cannot be higher than 1')

                motifs[motif_specie] = float(motif_score)

        # Opens tree file into ModTree format.
        tree = open(str(options.tree), 'r')
        my_tree = ModTree(tree.readline(), 1)

        return my_tree, motifs
    
    def initialize_nodes(tree, motifs):
        """
        Description:
        Initializes leaves nodes to motif scores and inner tree nodes to theoretical values.
        
        Parameters:
        newick tree and dictionary motifs.

        Returns:
        None
        """

        # Raises error if the tree leaves and motif sites names do not match.
        if match(tree, motifs):

            for node in tree.traverse('levelorder'):

                if node.is_root():     # Root node
                    node.prob             = None
                    node.complement_prob  = 1 # In theory, it should be 1.
                    node.effective_length = node.dist # by theory, correct

                else:
                    if node.is_leaf(): # Leaf nodes
                        node.prob             = motifs[str(node.name)]
                        node.complement_prob  = None
                        node.effective_length = node.dist
                        
                    else:              # Middle nodes
                        node.prob             = None
                        node.complement_prob  = None
                        node.effective_length = node.dist
        else:
            raise ValueError('Tree leaves and motif species do not match')

    def tree_probabilities(tree): # Tree Probability
        """
        Description:
        
        Parameters:

        Returns:

        """

        for node in tree.traverse("postorder"):

            prob = 1
        
            if not node.is_leaf():
            
                for child_node in node.get_children():
                    prob = prob * (1 - child_node.prob)

                node.prob = 1 - prob
            
    def effective_length(tree):
        '''
        do stuff for EL
        '''
        for node in tree.traverse("postorder"):
            len_prob = 0

            if not node.is_leaf():

                for child_node in node.get_children():
                    len_prob = len_prob + child_node.prob * child_node.effective_length

            try:
                node.effective_length = node.effective_length + (len_prob / node.prob)

            except ZeroDivisionError:
                node.effective_length = 0

    def tree_complement_probabilities(tree):
        '''
        docstring method
        '''
        for node in t.traverse("preorder"):

            if node.is_root():
                node.complement_prob = 1

            if not node.is_root():
                parent_complement_prob = node.up.complement_prob
                #print parent_complement_prob
                sister_prob = 1

                for sister_node in node.get_sisters():
                    sister_prob = sister_prob * (1 - sister_node.prob)

                node.complement_prob = parent_complement_prob * sister_prob

    def bbls(tree):
        '''
        docstring method
        '''
        BBLS = 0

        for node in tree.traverse():

            if not node.is_leaf():

                node_BBLS = node.complement_prob        

                children_node_length  = 0

                for child_node in node.get_children():

                    node_BBLS = node_BBLS * child_node.prob

                    children_node_length = children_node_length + child_node.effective_length

                node_BBLS = node_BBLS * children_node_length

                BBLS = BBLS + node_BBLS

        return BBLS


    def bls(my_tree, motifs):
        '''
        docstring method
        '''

        tree_ = copy.copy(my_tree)
        motifs_ = copy.copy(motifs)

        for motif in motifs_:
            motifs_[motif] = 1

        initialize_nodes(tree_, motifs_)
        tree_probabilities(tree_)
        effective_length(tree_)
        tree_complement_probabilities(tree_)

        return bbls(tree_)


    t, motifs = open_format_input(options)    
    initialize_nodes(t, motifs)
    tree_probabilities(t)
    effective_length(t)
    tree_complement_probabilities(t)

    print('BBLS: ' + str(bbls(t)))
    print('BLS: ' + str(bls(t, motifs)))


if __name__ == '__main__':

    main(args, options)

'''
    try:
        if options.help:
            doc = globals()['__doc__'] % {'usage' : USAGE % PROG_NAME, 'version' : VERSION,\
                                          'progname' : PROG_NAME}
            #pager(doc)
            sys.exit(0)
        if not options.l:
            parser.print_usage()
            sys.exit(0)

    except KeyboardInterrupt:
        sys.stderr.write('\n')
        sys.stderr.flush()
        sys.exit(2)
    except SystemExit:
        pass

    except:
        sys.stderr.write('Fatal Error\n')
        sys.exit(2)
'''
