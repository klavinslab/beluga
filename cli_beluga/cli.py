import click
import os, sys
import errno
import logging
import parse_input
import util
import beluga_obj
#import matplotlib.pyplot as plt 
#import numpy as np
#import pylab
#import random
#from numpy import loadtxt
#from sympy import init_session, init_printing, Matrix
#from shutil import copyfile
#init_session()
#from scipy.integrate import odeint

all_colors = 'black', 'red', 'green', 'yellow', 'blue', 'magenta', \
             'cyan', 'white'

@click.command()
@click.option('--grammar', '-g', type=click.File('rb'), help='The input grammar file')
@click.option('--config', '-c', type=click.File('rb'), help='The input json file')


# def preorder(tree):
#     if tree:
#         print(tree.getRootVal())
#         preorder(tree.getLeftChild())
#         preorder(tree.getRightChild())

def main(grammar, config):
    """Automated, experimentally-driven design of genetic circuits."""
    greet = 'Hello'
    click.echo(click.style('{0}, {1}.'.format(greet, 'beluga user!'), fg='blue', blink=True))
    rd_grammar = grammar.read()
    rd_config = config.read()
    click.echo(rd_grammar)
    click.echo(rd_config)
    
    logging.basicConfig(filename='beluga.log',level=logging.DEBUG)
    logging.info('Welcome to BELUGA!')
    logging.info(rd_grammar)
    logging.info(rd_config)
    #logging.debug('This message should go to the log file')
    #logging.warning('And this, too')
    

    #set of non terminals
    N = ('<subject>', '<predicate>', '<noun phrase>', '<noun>', '<article>', '<verb>',   '<direct object>')
    #set of teminals

    T = ('the', 'boy', 'dog', 'bit')
    T = ('p0','p1','p2','g0','g1','g2')
    #productions
    P = [ ('Sigma',           ['<subject>', '<predicate>']), \
    ('<subject>',       ['<noun phrase>']),            \
    ('<predicate>',     ['<verb>']),                   \
    ('<predicate>',     ['<verb>','<direct object>']), \
    ('<noun phrase>',   ['<article>','<noun>']),       \
    ('<direct object>', ['<noun phrase>']),            \
    ('<noun>',          ['boy']),                      \
    ('<noun>',          ['dog']),                      \
    ('<article>',       ['the']),                      \
    ('<verb>',          ['bit'])                       ]

    # P = [ ('Design',           ['R']),                  \
    # ('Design',          ['C','R']),                      \
    # ('Design',          ['C','C','R']),                  \
    # ('C',               ['Prom', 'Gene']),              \
    # ('R',               ['Prom', 'g0']),                \
    # ('Prom',            ['p0']),                        \
    # ('Prom',            ['p1']),                        \
    # ('Prom',            ['p2']),                        \
    # ('Gene',            ['g1']),                        \
    # ('Gene',            ['g2'])                         ]


    P = [ ('Design',           ['R']),                  \
    ('Design',          ['C','R']),                      \
    ('C',               ['Prom', 'Gene']),              \
    ('R',               ['Prom', 'g0']),                \
    ('Prom',            ['p0']),                        \
    ('Prom',            ['p1']),                        \
    ('Gene',            ['g1'])                         ]



    input_new = parse_input.algorithm_info(P,T)
    
    #logging.info(input_new.language)    
    logging.info("Number of Designs:")
    logging.info(len(input_new.language))

    my_beluga = beluga_obj.beluga_obj(input_new.language,input_new.design_objective)
    logging.info("Design Space:")
    logging.info(my_beluga.design_space)
    my_beluga.search()

