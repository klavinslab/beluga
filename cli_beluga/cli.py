import click
import os, sys
import errno
import logging
import parse_input
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
#@click.argument('input', type=click.File('rb'), nargs=-1)
#@click.argument('output', type=click.File('wb'))
#@click.argument('name', default='world', required=False)
#@click.argument('input_file', type=click.File('rb'))
#@click.argument('out',type=click.File('w'), required=False)
def main(grammar, config):
    """Automated, experimentally-driven design of genetic circuits."""
    greet = 'Hello'
    #click.echo(out)
    click.echo(click.style('{0}, {1}.'.format(greet, 'beluga user!'), fg='blue', blink=True))
    click.echo(grammar.read())
    click.echo(config.read())
    
    logging.basicConfig(filename='beluga.log',level=logging.DEBUG)
    logging.info('Welcome to BELUGA!')
    #logging.debug('This message should go to the log file')
    #logging.warning('And this, too')
    
#    input_new = parse_input.algorithm_info(input_file)

#    my_beluga = beluga_obj.beluga_obj(...)
#    my_beluga.search(...)

