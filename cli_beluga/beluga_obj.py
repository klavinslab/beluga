import emcee
import corner
import numpy as np
from sympy import *
import os, sys
import contextlib
import errno
import util
import random
import math
import matplotlib.pyplot as pl
from matplotlib.pyplot import cm
import scipy.optimize as op
from scipy.optimize import curve_fit
from scipy.optimize import leastsq
from scipy.integrate import odeint
from scipy.spatial import distance
from scipy.signal import savgol_filter
from matplotlib.ticker import MaxNLocator
from SALib.sample import saltelli
from SALib.analyze import sobol
from SALib.test_functions import Ishigami
from IPython.core.pylabtools import figsize
import scipy.stats as stats
import pymc as pm
import pylab as P

figsize(11, 9)


sub_model = []
sim_model = []
project_path = "/Users/Leli/beluga_cli_test_dir/ABCSMC"

sims_total = 0

def fileno(file_or_fd):
    fd = getattr(file_or_fd, 'fileno', lambda: file_or_fd)()
    if not isinstance(fd, int):
        raise ValueError("Expected a file (`.fileno()`) or a file descriptor")
    return fd

@contextlib.contextmanager
def stdout_redirected(to=os.devnull, stdout=None):
    """
    http://stackoverflow.com/a/22434262/190597 (J.F. Sebastian)
    """
    if stdout is None:
       stdout = sys.stdout

    stdout_fd = fileno(stdout)
    # copy stdout_fd before it is overwritten
    #NOTE: `copied` is inheritable on Windows when duplicating a standard stream
    with os.fdopen(os.dup(stdout_fd), 'wb') as copied: 
        stdout.flush()  # flush library buffers that dup2 knows nothing about
        try:
            os.dup2(fileno(to), stdout_fd)  # $ exec >&to
        except ValueError:  # filename
            with open(to, 'wb') as to_file:
                os.dup2(to_file.fileno(), stdout_fd)  # $ exec > to
        try:
            yield stdout # allow code to be run with the redirected stdout
        finally:
            # restore stdout to its previous value
            #NOTE: dup2 makes stdout_fd inheritable unconditionally
            stdout.flush()
            os.dup2(copied.fileno(), stdout_fd)  # $ exec >&copied

def func(x, a, b, c):
	return a * np.exp(-b * x) + c

def is_sublist(a,b):
	if a == []:
		return True
	if b == []:
		return False
	return b[:len(a)] == a or is_sublist(a,b[1:])

def printMDPQueue(mdpQ):
	if mdpQ.isEmpty() != 1:
		for element in mdpQ:
			print element[0], element[1].symb_param
	return 0

#output: sums the entries in the cols of a matrix, for each row
def sum_cols(mat):
    
    sum_mat = np.zeros([len(mat)])
    
    for j in range(0,len(mat)):
        col_sum = 0
        for k in range(0,len(mat[j])):
            col_sum+=sum_mat[j]+mat[j][k]
        sum_mat[j] = col_sum       
    return sum_mat

def odeint_model(x_vec,t):
	global sub_model
	G0, G1, G2 = symbols("G0, G1, G2")
	temp_g0, temp_g1 = x_vec
	#print "This is the SUBMODEL in ode int fcn: ", sub_model
	#print "and this is temp_g0 and temp_g1 ", temp_g0, temp_g1
	new_values = list(sub_model.subs({G0: temp_g0, G1: temp_g1}))
	#print "new values", new_values
	return new_values

def gen_ODE_file(design_name, design, dir_str):
    global project_path
    
    
    #CREATE DESIGN PYTHON FILE
    if (dir_str == "modelSelection"):
        file_name = project_path + "/modelSelection" + "/design_" + design_name + ".py"
    elif (dir_str == "paramInference"):
        file_name = project_path + "/paramInference" + "/design_" + design_name + ".py"
    
    file = open(file_name, "w")

    file.write("from math import sqrt\n")
    file.write("from numpy import random\n")
    file.write("from abcsysbio.relations import *\n\n")
    
    file.write("def modelfunction((")
    
    for i in range(0,len(design.model)):
        file.write("G"+str(i)+",")
    
    file.write("),time,parameter=(")
    
    for j in range(0,len(design.params)):
        file.write("0,")
    
    file.write(")):\n")
    file.write("\tS = 100\n\n")
    
    for j in range(0,len(design.params)):
        file.write("\t"+str(design.params[j])+"=parameter["+str(j)+"]\n")
        
    for j in range(0, len(design.model)):
        file.write("\td_G"+str(j)+"="+str(design.model[j])+"\n")
        
    file.write("\treturn(")
    for j in range(0,len(design.model)):
        file.write("d_G"+str(j)+",")
    file.write(")\n\n")
    
    file.write("def rules((")
    
    for j in range(0,len(design.model)):
        file.write("G"+str(j)+",")
    file.write("),(")
    
    for j in range(0,len(design.params)):
        file.write(str(design.params[j])+",")
    file.write("),time):\n\n")
    
    file.write("\treturn((")
    for j in range(0,len(design.model)):
        file.write("G"+str(j)+",")
    file.write("),(")
    
    for j in range(0,len(design.params)):
        file.write(str(design.params[j])+",")
    file.write("))\n")
    

        
    file.close()


def genABCSMCInputFile(param_space, tester_data, model_str, models, dir_str, sim=False, 
                        data=([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])):
    
    #CREATE TEMPLATE INPUT_FILE FOR MODEL SELECTION BETWEEN THESE HYPOTHESES
    print "FILE NAME = ", model_str
    if (dir_str == "modelSelection"):
        file_name = project_path + "/modelSelection" + "/input_models"
        for j in range(0,len(models)):
            file_name += "_" + str(models[j])
        file_name += ".xml"
        cur_data = desired_data
    elif (dir_str == "paramInference"):
        file_name = project_path + "/paramInference" + "/input_models"
        file_name += "_" + model_str + ".xml"
        cur_data = data

    #print "cur_data = ", cur_data

    
    file = open(file_name, "w")

    file.write("<input>\n")
    
    if (dir_str == "modelSelection"):
        file.write("<modelnumber> 3 </modelnumber>\n")
        file.write("<restart> False </restart>\n")
        file.write("<epsilon>\n")
        file.write("<e1> 50 20 15 10 5 3 2.5 2 1.7 1.5 1.2 1 </e1>\n")
        file.write("</epsilon>\n")
    elif (dir_str == "paramInference"):
        file.write("<modelnumber> 1 </modelnumber>\n")
        file.write("<restart> False </restart>\n")
        file.write("<autoepsilon>\n")
        file.write("<finalepsilon> 0.01 </finalepsilon>\n")
        file.write("<alpha> 0.5 </alpha>")
        file.write("</autoepsilon>\n")
    
    file.write("<particles> 500 </particles>\n")
    file.write("<beta> 1 </beta>\n")
    file.write("<dt> 0.1 </dt>\n")
    file.write("<kernel> uniform </kernel>\n")
    file.write("<modelkernel> 0.7 </modelkernel>\n")
    file.write("<rtol> 1e-6 </rtol>\n")
    file.write("<atol> 1e-6 </atol>\n")
    file.write("<data>\n")
    file.write("<times>")
    
    for j in range(0, len(cur_data[0])):
        file.write(" "+str(cur_data[0][j]))
        
    file.write(" </times>\n")
    file.write("<variables>\n")
    file.write("<var1>")
    
    for j in range(0, len(cur_data[1])):
        file.write(" "+str(cur_data[1][j]))
    
    file.write(" </var1>\n")
    file.write("</variables>\n")
    file.write("</data>\n")
    file.write("################## Models\n")
    
    file.write("<models>\n")
    
    if (dir_str == "modelSelection"):
        r_val = len(models)
    elif (dir_str == "paramInference"):
        r_val = 1
    
    
    for i in range(0, r_val):
        
        if (dir_str == "modelSelection"):
            model_index = models[i]
        elif (dir_str == "paramInference"):
            model_index = models  
        
        file.write("<model"+str(i+1)+">\n")
        file.write("<name> design_"+str(model_str)+" </name>\n")
        file.write("<source> design_"+str(model_str)+".py </source>\n")
        file.write("<type> ODE </type>\n")
        file.write("<fit> species1 </fit>\n")
        file.write("<initial>\n")
        for m in range(0, len(models.model)):
            file.write("<ic"+str(m+1)+"> constant 1.0 </ic"+str(m+1)+">\n")
        file.write("</initial>\n")
        file.write("<parameters>\n")
        for j in range(0, len(models.params)):
            if sim==False:
                file.write("<"+str(models.params[j])+
                       "> uniform "+ str(param_space[models.params[j]][0]) +" "+ str(param_space[models.params[j]][1])+" </"+str(models.params[j])+">\n")
            else:
                file.write("<"+str(models.params[j])+
                       "> constant "+ str(tester_data[str(models.params[j])]) +" </"+str(models.params[j])+">\n")
        file.write("</parameters>\n")
        file.write("</model"+str(i+1)+">\n")
            
    file.write("</models>\n")
    file.write("</input>")
    
    file.close()

def lnprior(theta):
    m, b, lnf = theta
    if -5.0 < m < 0.5 and 0.0 < b < 10.0 and -10.0 < lnf < 1.0:
        return 0.0
    return -np.inf
def lnlike(theta, x, y, yerr):
    m, b, lnf = theta
    #print "\nChosen params: ", theta
    model = m * x + b
    inv_sigma2 = 1.0/(yerr**2 + model**2*np.exp(2*lnf))
    res = -0.5*(np.sum((y-model)**2*inv_sigma2 - np.log(inv_sigma2)))
    # print "model = ", model
    # print "y = ", y
    # print "Result: ", res
    return res
def lnprob(theta, x, y, yerr):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, x, y, yerr)

def mymichelis_menten(self, x_vec, t, *args):
	#print args
	#print args[0]
	#print args[1]
	global sub_model
	params = {"L00": args[0], "K00": args[1]}
	sub_model = self.design_space["['p0', 'g0']"].model.subs(params)
	G0, G1, G2 = symbols("G0, G1, G2")
	temp_g0, temp_g1 = x_vec
	new_values = list(sub_model.subs({G0: temp_g0, G1: temp_g1}))
	new_values = 0
	return new_values

def mylnprior(theta):
    L00, K00 = theta
    if 0.0 < L00 < 10.0 and 0.0 < K00 < 10.0:
        return 0.0
    return -np.inf

# def mylnlike(theta, x, g0):
#     #L00, K00 = theta
#     init_species = [1,1,1]
#     t_out = np.arange(0,25,1)
#     p = tuple(theta)
#     print "chosen params: ", p
#     params = p[0], p[1]
#     print "params again: ", params
#     sim_P = odeint(mymichelis_menten, init_species, t_out, args=params).flatten()
#     #res = distance.euclidean(np.array(sim_P.T[0]), np.array(g0))
#     res = np.sum((np.array(sim_P.T[0]) - np.array(g0))**2)
#     return res
    #return sim_P.T[0] - g0
    # model = m * x + b
    # inv_sigma2 = 1.0/(yerr**2 + model**2*np.exp(2*lnf))
    #return -0.5*(np.sum((y-model)**2*inv_sigma2 - np.log(inv_sigma2)))

def mylnprob(theta, x, g0):
    lp = mylnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + mylnlike(theta, x, g0)


class design:
	def __init__(self, des_str):
		self.id = des_str
		self.successors = []
		self.model_file = ""
		self.model = []
		self.params = []
		self.char_data = []

	def isSubDesign(self, desb):
		result = False
		cur_des_id = self.id
		test_des_id = desb.id
		cur_des_idx = 0

		# a = self.params
		# b = desb.params
		# if a == []:
		# 	return True
		# if b == []:
		# 	return False
		# return b[:len(a)] == a or a.isSubDesign(b[1:])

		# for i in range(len(test_des_id)):
		# 	#print "i = ", i
		# 	if test_des_id[i] == cur_des_id[cur_des_idx]:
		# 		#print "elements equal: ", test_des_id[i], cur_des_id[cur_des_idx]
		# 		cur_des_idx+=1
		# 	else:
		# 		cur_des_idx = 0
		# 		if test_des_id[i] == cur_des_id[cur_des_idx]:
		# 			#print "elements equal: ", test_des_id[i], cur_des_id[cur_des_idx]
		# 			cur_des_idx+=1
		# if cur_des_idx == len(cur_des_id):
		# 	result = True
		# else:
		# 	result = False

		# return result

	def isSimilar(self, desb):
		result = False
		desa_id = self.id
		desb_id = desb.id
		desa_len = len(desa_id)
		desb_len = len(desb_id)

		match_count = 0

		# if (abs(desa_len - desb_len)>1):
		# 	return False
		if desa_len != desb_len:
			return False
		else:
			for i in range(desa_len):
				if desa_id[i] == desb_id[i]:
					match_count += 1
			if abs(match_count - desa_len) > 1:
				result = False
			else:
				result = True

		return result


	def getModel(self,species,promoters,Ls,Ks,Gs):
	    params = []
	    odes = []
	    
	    #creating identity matrix for design (prom/gene only)
	    mat = np.zeros([len(species),len(promoters)])
	    #print "Species = ", species
	    for gene in species:
	    	#print gene, " in ", self.id, " = ", (gene in self.id), " at ", self.id.index(gene)
	    	if gene in self.id:
	    		#print "GENE is, ", gene
	    		indices = [d for d, x in enumerate(self.id) if x == gene]
	    		#print indices
	    		for index in indices:
		    		#rel_prom_idx = self.id.index(gene) - 1
		    		rel_prom_idx = index - 1
		    		prom = self.id[rel_prom_idx]
		    		prom_idx = promoters.index(prom)
		    		gene_idx = species.index(gene)
		    		mat[gene_idx][prom_idx] = 1

	    #print "matrix for ", self.id, " = ", mat
	    sum_mat = sum_cols(mat)

	    for n in species:
	    	odes.append(0)
	    #### FIX ME ####
	    ### Look at design 193 ... 

	    #GENERATE GENERIC ODE MODEL FROM DESIGN MATRIX
	    for a in range(0, len(species)):
	        for j in range(0, len(promoters)):
	            if (mat[a][j]!=0):
	                for i in range(0, len(species)):
	                    if (sum_mat[[i]]>0):
	                        odes[a]+= (Ls[j][i]) / (1+Gs[0][i]*Ks[i][j])
	                        Lparam_exists = Ls[j][i] in params
	                        Kparam_exists = Ks[i][j] in params
	                        #print Ls[j][i], Ks[i][j], " -> ", Lparam_exists, Kparam_exists
	                        if Lparam_exists == False:
	                            params.append(Ls[j][i])
	                        if Kparam_exists == False:
	                            params.append(Ks[i][j])

	        odes[a]+=-Gs[0][a]
	 
	    #model = Matrix([odes[0]]+[odes[1]]+[odes[2]])
	    model = Matrix([odes[0]]+[odes[1]])
	    
	    #print params

	    # global global_params2
	    # for x in params:
	    #     global_params2.update({x: [0,10]})
	    self.model = model
	    self.params = params

	    return (model,params)

	def testDesign(self, cur_param_sp, ground_truth_vals):
		param_sp = cur_param_sp
		global sub_model

		sub_model = self.model.subs(ground_truth_vals)
		init_species = [10,1]
		t_out = np.arange(0,25,1)
		sim_out = odeint(odeint_model,init_species,t_out)
		#print "SIMULATION of design ", des, ":"
		g0_res, g1_res = sim_out.T
		#print list(g0_res)

		self.char_data = ([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24],list(g0_res))
		#design_name = str(self.id).strip('[]').strip(', ')
		design_name = ""
		for elem in self.id:
			design_name+=str(elem)
		gen_ODE_file(design_name, self, "paramInference")
		genABCSMCInputFile(param_sp, ground_truth_vals, design_name, self, "paramInference", False, self.char_data)
		print "PLEASE RUN THE FOLLOWING: \n"
		print "\t cd /Users/Leli/Documents/py_try/ABC_sysbio/abc-sysbio-2.08/working_ex/beluga_12_7/paramInference \n"
		print "\t run-abc-sysbio -i input_models_"+ design_name + ".xml" + " -lc -of=results_" + design_name + " -f -sd=2 \n"
		stuff = input("Commands run? ")

		return param_sp
		#simulate with ground_truth_vals and then
		#do param inference on self design using it's cur_param_sp and using char_data from sim
		#and update the cur_param_sp
	
	def getError(self, goal):
		return 0


def getSuccessors(des,design_dict):
	successors = []
	sub_des_list = []
	exp_des_list = []
	#suc_exp = []
	new_successors = []

	#Find the build successors

	for key, nxtdes in design_dict.iteritems():
		if is_sublist(des.params,nxtdes.params) and (des.id != nxtdes.id):
			sub_des_list.append((nxtdes, len(nxtdes.params)))
		
	if len(sub_des_list) > 0:
		min_param_len = min(sub_des_list,key=lambda item:item[1])[1]
	for item in sub_des_list:
		if item[1] == min_param_len:
			successors.append(str(item[0].id))
		
	if len(successors) > 0:
		min_des_len = len(min(successors,key=lambda item:len(item)))
		for thing in successors:
			if len(thing) == min_des_len:
				new_successors.append(thing)



	#Find the explore successors
	for key2, nxtdes2 in design_dict.iteritems():
		if len(set(des.params).intersection(set(nxtdes2.params))) == 0:
			exp_des_list.append((nxtdes2, len(nxtdes2.params)))
	if len(exp_des_list) > 0:
		min_exp_len = min(exp_des_list,key=lambda item:item[1])[1]
	for item in exp_des_list:
		if item[1] == min_exp_len:
			new_successors.append(str(item[0].id))



	# #Find all designs that des is a subdesign of
	# for key, subdes in design_dict.iteritems():
	# 	if des.isSubDesign(subdes):
	# 		if des.id != subdes.id:
	# 			sub_des_list.append((subdes,len(subdes.id)))


	# #Finding successors: choose min subdes's
	# #MAY CHANGE THIS TO ACCEPT ANY DESIGN WHICH IS A SUPERSET OF THIS DES AND NOT JUST THE NEXT CLOSEST ONE
	# if len(sub_des_list) > 0:
	# 	min_des_len = min(sub_des_list,key=lambda item:item[1])[1]
	# for item in sub_des_list:
	# 	if item[1] == min_des_len:
	# 		successors.append(str(item[0].id))

	# #Now find all designs that differ only by one element from des
	# for key, otherdes in design_dict.iteritems():
	# 	if des!=otherdes and des.isSimilar(otherdes) and (str(otherdes.id) not in successors):
	# 		successors.append(str(otherdes.id))

	return new_successors

class MDP_node:
	"""
	A class for the MDP formulation of our search problem
	"""

	def __init__(self, design_indx, design_space):
		self.my_design = design_indx
		self.symb_param = {} #self.getSymbParams(design_space)
		self.actions = {} #self.getActions(design_space)
		self.successors = []
		self.percent_unknown = 1
		self.history = []
		self.policy = []
		self.terminal_reward = (0,-1)
		self.vVal = 0

	def getSymbParams(self, design_space):
		return []

	def getPossibleActions(self, design_dict):
		known_params = []
		possible_actions = {}
		#if self.perct_unknown <=0:
		if self.percent_unknown <= 0:
			return {}

		#gets a list of the known params only
		cur_known = self.symb_param
		for key, item in cur_known.iteritems():
			if item == 0:
				known_params.append(key)

		known_set = set(known_params)

		#find build actions
		successors = []
		sub_des_list = []
		exp_des_list = []
		#suc_exp = []
		actions = {}

		#Find the build successors
		for key, nxtdes in design_dict.iteritems():
			#Should this check the full history??? Or just the last action that got us here?
			if known_set.issubset(set(nxtdes.params)) and (str(nxtdes.id) not in self.history):
				#print known_set , "is a subset of ", set(nxtdes.params), "i.e design ", nxtdes.id
				sub_des_list.append((nxtdes, len(nxtdes.params)))
			# if is_sublist(known_params,nxtdes.params) and (nxtdes.id not in self.history):
			# 	sub_des_list.append((nxtdes, len(nxtdes.params)))
		
		###FIX ME###
		#change this to a percentage knowledge thing
		#instead of just the minimum design that is a superset of the set of known params 
		#say the successors are the set of designs in which we know at least 60% of the params
		#given our current knowledge ... and then this percentage can be tuned later to vary
		#simulation time vs number of experiments

		if len(sub_des_list) > 0:
			min_param_len = min(sub_des_list,key=lambda item:item[1])[1]
		for item in sub_des_list:
			#if item[1] == min_param_len:
			#we know at least 50% of the params
			#print "perct known for next considered design = ", float(float(len(known_params))/float(len(item[0].params)))
			if float(float(len(known_params))/float(len(item[0].params))) >= 0.5:
				successors.append(str(item[0].id))
		
		#then choose only smallest designs ... why do we need this??	
		if len(successors) > 0:
			min_des_len = len(min(successors,key=lambda item:len(item)))
			for thing in successors:
				#if len(known_set) 
				if len(thing) == min_des_len:
					actions.update({thing: []})

		#find explore actions
		for key2, nxtdes2 in design_dict.iteritems():
			if len(known_set.intersection(set(nxtdes2.params))) == 0:
				exp_des_list.append((nxtdes2, len(nxtdes2.params)))
		if len(exp_des_list) > 0:
			min_exp_len = min(exp_des_list,key=lambda item:item[1])[1]
		for item in exp_des_list:
			if item[1] == min_exp_len:
				actions.update({str(item[0].id): []})



		#print "actions for state ", self. symb_param, " = ", actions

		return actions

	def isTerminal(self):
		return 0
	

class beluga_obj:
    """
    A class to setup and run the beluga algorithm.
    
    """ 
    # instantiation
    def __init__(self, language, design_obj):
    	#self.species = ['g0','g1','g2']
    	self.species = design_obj.species
    	self.promoters = ['p0','p1']
    	self.goal = design_obj.output
    	self.testdata = design_obj.testdata
    	self.Ls = []
    	self.Ks = []
    	self.Gs = []
    	self.param_space = {}
    	self.design_space = self.genGraph(language)
    	self.belugaMDP = self.initMDP()
    	self.MDP_rewards = {}
    	#self.policy = self.initPolicy()
    	
    	

    def genGraph(self, language):
    	design_dict = {}
    	global sub_model
    	for elem in language:
    		str_elem = str(elem)
    		design_dict.update({str_elem: design(elem)})


    	#Generating full space of params and symbols
    	for m in range(0,len(self.species)):
			self.Ks.append(map(lambda i: Symbol('K' + str(m)+str(i)), xrange(len(self.promoters))))
			self.Ls.append(map(lambda i: Symbol('L' + str(m)+str(i)), xrange(len(self.promoters))))
			self.Gs.append(map(lambda i: Symbol('G' + str(i)), xrange(len(self.species))))
    	#finding all successors in one direction for each design in design list
    	for key, des in design_dict.iteritems():

    		#Generating ODE model and params:
    		des.getModel(self.species, self.promoters, self.Ls, self.Ks, self.Gs)
    		
    		####### TESTING ODE SIMULATION WITH SCIPY ODEINT ######
    		#######         THIS WILL NOT STAY HERE          ######
    		
    		# try:
    		# 	os.makedirs(project_path, 0777)
    		# except OSError as e:
    		# 	if e.errno != errno.EEXIST:
    		# 		raise e
    		# 	pass
    		# try:
    		# 	os.makedirs(project_path+"/paramInference", 0777)
    		# except OSError as e:
    		# 	if e.errno != errno.EEXIST:
    		# 		raise e
    		# 	pass
    		# try:
    		# 	os.makedirs(project_path+"/modelSelection", 0777)
    		# except OSError as e:
    		# 	if e.errno != errno.EEXIST:
    		# 		raise e
    		# 	pass

    		#########################^^^^##########################
    		############## THIS ALSO WON'T STAY HERE ##############
    		#######################################################

    		for x in des.params:
    			self.param_space.update({x: [0,10]})


    	for key1, des1 in design_dict.iteritems():
    		#Finding successors:
    		successors = getSuccessors(des1,design_dict)
    		for item in successors:
    			if item not in des1.successors:
    				des1.successors.append(item)
    			#updating the successor graph edges in the other direction
    			# if des.id not in design_dict[item].successors:
    			# 	design_dict[item].successors.append(des.id)

    	for key, thing in design_dict.iteritems():
    		print " key = ", key, type(key)
    		print "\n id = ", thing.id, type(thing.id)
    		print "\n model = ", thing.model
    		print "\n params = ", thing.params
    		print "\n successors = ", thing.successors
    		print "-----------------------------------\n\n"

    	print "Num Designs in Space = ", len(design_dict)
    	print "Param space: ", self.param_space
    	#print "len of intersection of [p0,g0] and [p0,g1,p0,g0] params: ", len(set(design_dict["['p0', 'g0']"].params).intersection(set(design_dict["['p1', 'g0']"].params)))
    	#print "design p0,g1|p1,g0 model: ", design_dict["['p0', 'g1', 'p1', 'g0']"].model
    	return design_dict

    def getTransitionsandProb(self,state,action):
    	"""
        Returns list of (nextState, prob) pairs
        representing the states reachable
        from 'state' by taking 'action' along
        with their transition probabilities.
        """
    	return 0

    def initMDP(self):
    	#an MDP state has associated with it:
    		#a design (indx in des_space)
    		#a symbolic set of params it could learn, 0 learned, 1 not learned
    		#actions it can take (the neighbors of the design) and their set of possible s' states
    	#start with self.design_space[0]
    	mdp_Q = util.Queue()
    	mdp_dict = {}
    	s_key = 0
    	start_node = MDP_node(None, self.design_space)
    	start_node.symb_param = dict(self.param_space) #copy all keys in global param list and put a 1 value next to it
    	for p_key, p in start_node.symb_param.iteritems():
    		start_node.symb_param.update({p_key: 1})
    	#print "Start node symb params = ", start_node.symb_param
    	#start_node.actions = {"['p0', 'g0']": [], "['p1', 'g0']": [], "['p2', 'g0']": []}
    	start_node.history.append("None")
    	start_element = ("S"+str(s_key), start_node)
    	mdp_Q.push(start_element)

    	#WHAT TO DO:
    	#figure out all the states in the MDP 
    		#for each state call getTransitionsandProb to determine successors
    		#just basically do what i was doing before but the successors will be different
    		#should probably just redo building of the inital design space state graph
    		#basically the getsuccessors needs to look at params NOT parts! But otherwise it's the same.
    	#ALSO, determine initial policies for each state by choosing simplest action

    	while mdp_Q.isEmpty()==0:
    		mdpelement = mdp_Q.pop()
    		key, mdpnode = mdpelement[0], mdpelement[1]
    		mdp_dict.update({key: mdpnode})
    		print "\n\nMDP NODE = ", key
    		#print "actions before = ", mdpnode.actions
    		
    		param_sum = 0
    		
    		for v1_key, val1 in mdpnode.symb_param.iteritems():
    			param_sum += val1
    		perct_unk = float(float(param_sum)/float(len(mdpnode.symb_param)))
    		#print "Percent unknown = ", perct_unk

    		actions = mdpnode.getPossibleActions(self.design_space)

    		for action, val in actions.iteritems():
    			new_mdpnode = MDP_node(action, self.design_space)
    			new_mdpnode.history = list(mdpnode.history)
    			new_mdpnode.history.append(action)

    			a_des_params = self.design_space[action].params
    			s_key += 1
    			new_symb_param = dict(mdpnode.symb_param)
    			for param in a_des_params:
    				new_symb_param[param] = 0
    			new_mdpnode.symb_param = new_symb_param

    			sym_param_sum = 0
    			for v_key, val in new_mdpnode.symb_param.iteritems():
    				sym_param_sum += val
    			perct_unknown = float(float(sym_param_sum)/float(len(new_symb_param)))
    			new_mdpnode.percent_unknown = perct_unknown

    			###FIX ME###
    			#This way of getting prob doesn't seem right
    			prob = 1 - (perct_unk-perct_unknown)
    			
    			actions[action].append(("S"+str(s_key),prob,(0,-1)))
    			actions[action].append((key,1-prob,(0,-1)))
    			mdpnode.actions = dict(actions)
    			new_element = ("S"+str(s_key), new_mdpnode)
    			#print "PUSHING ELEMENT", new_element[0], "TO QUEUE"
    			mdp_Q.push(new_element)
    					
    		print "actions after = ", mdpnode.actions
    		print "symb params = ", mdpnode.symb_param
    		print "history after = ", mdpnode.history

    	print "Num MDP states = ", len(mdp_dict)
    	return mdp_dict

    def initPolicy(self):
    	#start with MDP_graph[0] ... i.e. S0
    	#for every state in MDP_graph, determine leftmost action it can take (indx of first successor of design in des_space)
    	#update self.policy with {s:des_indx}
    	#actually i think i'll init to explore actions - or the simplest designs i guess
    	policy = {}
    	for key, item in self.belugaMDP.iteritems():
    		#if item.isTerminal() == False:
    		if len(item.actions) > 0:
    			policy.update({key: item.actions.keys()[len(item.actions)-1]})
    	#Start with start of MDP, follow the policy until it terminates and store and return this policy.
    	#print policy
    	return policy

    def learn(self,curstate,nxtstate):
    	total_param_range = 0
    	for p, pranges in self.param_space.iteritems():
    		total_param_range += pranges[1] - pranges[0]
    	learn = -1 * float(float(total_param_range)/float(10*len(self.param_space)))
    	print "learn val = ", learn
    	return learn

    def local_min(self,ys):
    	return [y for i, y in enumerate(ys)
            if ((i == 0) or (ys[i - 1] >= y))
            and ((i == len(ys) - 1) or (y < ys[i+1]))]

    def local_max(self,ys):
    	return [y for i, y in enumerate(ys)
            if ((i == 0) or (ys[i - 1] <= y))
            and ((i == len(ys) - 1) or (y > ys[i+1]))]

    def michelis_mentenz(self, x_vec, t, *args):
    	global sim_model
    	params = {"L00": args[0], "L01": args[1], "L10": args[2], "L11": args[3],
    	 "K00": args[4], "K01": args[5], "K10": args[6], "K11": args[7]}
    	sub_model = sim_model.subs(params)
    	#print sub_model
    	#sub_model = self.design_space["['p0', 'g1', 'p1', 'g0']"].model.subs(params)
    	sub
    	G0, G1 = symbols("G0, G1")
    	temp_g0, temp_g1 = x_vec
    	new_values = list(sub_model.subs({G0: temp_g0, G1: temp_g1}))
    	# print "actual model = ", self.design_space["['p0', 'g0']"].model
    	# print "model = ", sub_model
    	# print "sub_model = ", temp_g0, temp_g1
    	# print "the model = ", new_values
    	return new_values


    def goal_dist(self, action, iter_num):
    	#print "action = ", action
    	global sub_model
    	global sim_model
    	global sims_total
    	num_sims = 100
    	num_accepted = 0
    	epsilon_start = 70
    	sims_total = sims_total + 1
    	if epsilon_start - (iter_num*10) < 10:
    		epsilon = 10
    	else:
    		epsilon = epsilon_start - (iter_num*10)

    	sim_model = self.design_space[action].model

    	# print "\n\nIN GOAL DISTANCE. ACTION = ", action, ": "
    	# print "Remember, Goal is = ", self.goal[2]
    	print "\n\nCurrent EPSILON is = ", epsilon
    	print "ACTION IS : ", action
    	# print "And current design model is: ", sim_model
    	print "Current knowledge is: ", self.param_space
    	#print "\n"

    	names_arr = []
    	bound_arr = []
    	problem = {'num_vars': 0, 'names': [], 'bounds': []}

    	problem.update({'num_vars': len(self.design_space[action].params)})
    	for par in self.design_space[action].params:
    		names_arr.append(str(par))
    		bound_arr.append(self.param_space[par])
    	problem.update({'names': names_arr})
    	problem.update({'bounds': bound_arr})

    	#Y = np.empty(num_sims)

    	#print "Problem Is: ", problem

    	param_values = saltelli.sample(problem, 10, calc_second_order=False)
    	#print "Param values = ", len(param_values), param_values

    	#Y = np.empty([param_values.shape[0]])
    	#print "Y is now ", len(Y), Y
    	#print "working on submat ", problem['names'][0], param_values[0][0]
    	#num_sims = len(param_values)
    	sim_traj = []


    	#####
    	#Sensitivity analysis needs to be folded in here somehow. 

    	############################ PREVIOUSLY #############################
    	for i in range(0, num_sims):
    		submat = {}
    		for param, prange in self.param_space.iteritems():
    		# for j in range(0, len(param_values[i])):
    			submat.update({str(param) : random.uniform(prange[0],prange[1])})
    			# submat.update({str(problem['names'][j]) : param_values[i][j]})
    	
    		sub_model = self.design_space[action].model.subs(submat)
    		#print "current submat: ", submat
    		init_species = [10,1]
    		t_out = np.arange(0,25,1)
    		# l00_p = 0
    		# l01_p = 0
    		# l10_p = 0
    		# l11_p = 9.74
    		# k00_p = 0
    		# k01_p = 0
    		# k10_p = 0
    		# k11_p = 9.37

    		# mm_params = (l00_p, l01_p, l10_p, l11_p, k00_p, k01_p, k10_p, k11_p)
    		#Generate some synthetic data from the model
    		with stdout_redirected():
    			#sim_out = odeint(self.michelis_mentenz,init_species,t_out, args = (mm_params))
    			sim_out = odeint(odeint_model,init_species,t_out)
    			g0_res, g1_res = sim_out.T
    		sim_traj.append(g0_res)
    		dist = distance.euclidean(np.array(g0_res), np.array(self.goal[2]))

    		#np.append(Y,dist)
    		#Y[i] = dist

    		smooth_sim = savgol_filter(g0_res, 5, 3)
    		num_maxi = 0
    		global_min = min(smooth_sim)
    		local_max = self.local_max(smooth_sim)

    		for maxi in local_max:
    			if maxi > (global_min + 1):
    				num_maxi += 1
    		if num_maxi != 2:
    			dist += 10
    		#else:
    			#print "\nPulse!"

    		if dist <= epsilon:
    			num_accepted = num_accepted + 1
    		# print "SIM RESULTS: ", list(g0_res)
    		# print "GOAL: ", self.goal
    		# print "Distance = ", dist
    	goal_d = float(float(num_accepted)/float(num_sims))
    	print "\nnum_accepted = ", num_accepted, " ... num sims = ", num_sims
    	print "Goal_dist = ", float(float(num_accepted)/float(num_sims)), " for action ", action
    	#print "Y size = ", Y.size
    	#print "num vars (D) = ", problem['num_vars']
    	#print "Y % (D+2) = ", Y.size % (problem['num_vars'] + 2)
    	#Si = sobol.analyze(problem, Y, calc_second_order=False, print_to_console=False)

    	n = num_sims
    	color=iter(cm.rainbow(np.linspace(0,1,n)))
    	for i in range(n):
    		c = next(color)
    		pl.plot(t_out,sim_traj[i], c=c)
    	pl.plot(t_out, np.array(self.goal[2]),'b--')
    	pl.savefig("goal_dist_traj"+ str(iter_num)+".png")

    	# print "\nDoing param Sensitivity analysis"
    	# print problem['names']
    	# print Si['S1']
    	
    	return goal_d

    def getTerminalVal(self, state, iter_num):
    	action = self.belugaMDP[state].history[-1]
    	if iter_num == 1 :
    		if action not in self.MDP_rewards:
    			self.MDP_rewards.update({action: (iter_num, 0)})
    		elif self.MDP_rewards[action][0] != iter_num:
    			self.MDP_rewards.update({action: (iter_num, 0)})
    		return self.MDP_rewards[action][1]

    	if action not in self.MDP_rewards:
    		self.MDP_rewards.update({action: (iter_num, math.exp(iter_num+2) * self.goal_dist(action,iter_num))})
    	elif self.MDP_rewards[action][0] != iter_num:
    		self.MDP_rewards.update({action: (iter_num, math.exp(iter_num+2) * self.goal_dist(action,iter_num))})
    	
    	return self.MDP_rewards[action][1]

    def getReward(self, cur_state, action, next_state, iter_num):
    	"""
    	The reward returned by a state will be a weighted function of the number of new parameters it learned 
    	AND the average distance between a monte carlo simulation of action (design) a and the desired goal 
    	behavior (i.e the Error of action a). The weights for each element of the reward function will change 
    	as a function of time. To encourage exploration of unknown parameters early on, the weight for successfully 
    	learning new parameters will begin high and degrade over time. To encourage the agent to seek out a design 
    	which behaves to the user's defined behavior spec, the weight for the error will begin low and will 
    	increase significantly over time. I will have a funtion R(s,a,s') to represent rewards for each state.
    	"""
    	#So that we simulate all the states in the current policy only once before each policy eval/improvement
    	reward = 0
    	update = False
    	new_rew = ()
    	new_tuple = ()
    	if iter_num == 1 :
    		if action not in self.MDP_rewards:
    			self.MDP_rewards.update({action: (iter_num, self.learn(cur_state,next_state))})
    		elif self.MDP_rewards[action][0] != iter_num:
    			self.MDP_rewards.update({action: (iter_num, self.learn(cur_state,next_state))})
    		return self.MDP_rewards[action][1]

    	if action not in self.MDP_rewards:
    		self.MDP_rewards.update({action: (iter_num, self.learn(cur_state,next_state) + math.exp(iter_num) * self.goal_dist(action,iter_num))})
    	elif self.MDP_rewards[action][0] != iter_num:
    		self.MDP_rewards.update({action: (iter_num, self.learn(cur_state,next_state) + math.exp(iter_num) * self.goal_dist(action,iter_num))})
    	
    	return self.MDP_rewards[action][1]

    def findPolicy(self, policy, cur_round):
    	unchanged = False
    	cur_policy = policy
    	while unchanged==False:
    		#policy evaluation:
    		#for each state in the cur_policy, determine the value of the policy from following the cur_policy
    		#iterate until value of policy converges 
    		iterations = 100
    		discount = 0.9
    		values = util.Counter() #may change this to be read from an element in mdp node
    		#print "\nONE ITERATION: "
    		#print "Policy being evaluated: ", cur_policy
    		for k in range (0,iterations):
    			vVal_dict = util.Counter()

    			for s, action in cur_policy.iteritems():
    				qVal = 0
    				for t in self.belugaMDP[s].actions[action]:
    					t_prob = t[1]
    					nextstate = t[0]
    					#if nextstate = terminal state:
    					if len(self.belugaMDP[nextstate].actions) < 1:
    						values[nextstate] = int(self.getTerminalVal(nextstate, cur_round))
    					qVal += int(t_prob*(self.getReward(s, action, nextstate, cur_round)) + discount*values[nextstate])
    				vVal_dict[s] = qVal
    			
    			for s, action in cur_policy.iteritems():
    				values[s] = vVal_dict[s]
    		#print "Values of policy ", k , " - ", values, "\n\n"
    		#print "\nIn Find Policy - improvement: "
    		#policy improvement:
    		new_policy = dict(cur_policy)
    		#for each state in the cur_policy
    		for s, origaction in cur_policy.iteritems():
    			# val = 0
    			# maxVal = 0
    			#print self.belugaMDP[s].actions
    			for a, transitions in self.belugaMDP[s].actions.iteritems():
    				val = 0
    				for t in transitions:
    					t_prob = t[1]
    					nextstate = t[0]
    					if len(self.belugaMDP[nextstate].actions) < 1:
    						values[nextstate] = int(self.getTerminalVal(nextstate, cur_round))
    					val += int(t_prob*(self.getReward(s, a, nextstate, cur_round)) + discount*values[nextstate])

    				if val > values[s]:
    					new_policy.update({ s: a})

    		mismatch = 0
    		for s, curaction in cur_policy.iteritems():
    			if new_policy[s] != cur_policy[s]:
    				mismatch = mismatch + 1

    		if mismatch == 0:
    			unchanged = True
    		else:
    			cur_policy = dict(new_policy)

    	#print "DESE VALUES for the first policy, iteration ", k , " - ", values
    	print "CURRENT VALUES OF STATES:  ", values
    	return new_policy

    # def michelis_menten(self, x_vec, t, *args):
    # 	global sub_model
    # 	params = {"L00": args[0], "K00": args[1]}
    # 	sub_model = self.design_space["['p0', 'g0']"].model.subs(params)
    # 	G0, G1, G2 = symbols("G0, G1, G2")
    # 	temp_g0, temp_g1 = x_vec
    # 	new_values = list(sub_model.subs({G0: temp_g0, G1: temp_g1}))
    # 	print "actual model = ", self.design_space["['p0', 'g0']"].model
    # 	print "model = ", sub_model
    # 	print "sub_model = ", temp_g0, temp_g1
    # 	print "the model = ", new_values
    # 	return new_values

    def mylnprior(self, theta):
    	L00, K00 = theta
    	if 0.0 < L00 < 10.0 and 0.0 < K00 < 10.0:
        	return 0.0
    	return -np.inf

    def mylnlike(self, theta, x, g0):
    	#L00, K00 = theta
    	init_species = [10,1]
    	t_out = np.arange(0,25,1)
    	p = tuple(theta)
    	print "chosen params: ", p
    	params = p[0], p[1]
    	#print "params again: ", params
    	sim_P = odeint(self.michelis_menten, init_species, t_out, args=params)
    	#res = distance.euclidean(np.array(sim_P.T[0]), np.array(g0))
    	res = np.sum((np.array(sim_P.T[0]) - np.array(g0))**2)
    	print "sim 0 = ", sim_P.T[0]
    	print "sim = ", sim_P.T
    	print "chardata = ", g0
    	print "Res = ", res
    	return res
    	
    def mylnprob(self, theta, x, g0):
    	lp = self.mylnprior(theta)
    	if not np.isfinite(lp):
        	return -np.inf
    	return lp + self.mylnlike(theta, x, g0)

    def execute_policy1(self, policy, cur_state, cur_state_name):
    	np.random.seed(123)
    	#action = policy[cur_state_name]
    	action = "['p0', 'g0']"

    	##################### THE INTERNET EXAMPLE ##########################
    	#Choose the "true" parameters.
    	m_true = -0.9594
    	b_true = 4.294
    	f_true = 0.534
    	#Generate some synthetic data from the model
    	N = 50
    	x = np.sort(10*np.random.rand(N))
    	yerr = 0.1+0.5*np.random.rand(N)
    	y = m_true*x+b_true
    	y += np.abs(f_true*y) * np.random.randn(N)
    	y += yerr * np.random.randn(N)
    	print "this is the y : ", y
    	print "type : ", type(y)
    	
    	# global sub_model
    	# sub_model = self.design_space[action].model.subs(self.testdata)
    	# init_species = [1,1,1]
    	# t_out = np.arange(0,25,1)
    	# sim_out = odeint(odeint_model,init_species,t_out)
    	# #print "SIMULATION of design ", des, ":"
    	# g0_res, g1_res, g2_res = sim_out.T
    	print "(x,y,yerr)", (x,y,yerr)
    	print "type : ", type((x,y,yerr)[1])
    	# # Find the maximum likelihood value.
    	chi2 = lambda *args: -2 * lnlike(*args)
    	result = op.minimize(chi2, [m_true, b_true, np.log(f_true)], args=(x, y, yerr))
    	# Set up the sampler.
    	ndim, nwalkers = 3, 100
    	pos = [result["x"] + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
    	sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(x, y, yerr))
    	# Clear and run the production chain.
    	print("Running MCMC...")
    	sampler.run_mcmc(pos, 500, rstate0=np.random.get_state())
    	print("Done.")
    	# Compute the quantiles.
    	burnin = 50
    	samples = sampler.chain[:, burnin:, :].reshape((-1, ndim))
    	samples[:, 2] = np.exp(samples[:, 2])
    	m_mcmc, b_mcmc, f_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),zip(*np.percentile(samples, [16, 50, 84], axis=0)))
    	print("""MCMC result:
    		m = {0[0]} +{0[1]} -{0[2]} (truth: {1})
    		b = {2[0]} +{2[1]} -{2[2]} (truth: {3})
    		f = {4[0]} +{4[1]} -{4[2]} (truth: {5})
    	""".format(m_mcmc, m_true, b_mcmc, b_true, f_mcmc, f_true))
    	# #do param inference for design using char data
    	# #update self.param_space
    	# print "\n Doing experiment ", action, ": "
    	# print list(g0_res)

    	##################### THE PO,GO EXAMPLE ##########################
    	#Choose the "true" parameters.
    	l00_p = 2.256
    	k00_p = 5.658
    	mm_params = (l00_p, k00_p)
    	#Generate some synthetic data from the model
    	global sub_model
    	init_species = [10,1]
    	t_out = np.arange(0,25,1)
    	print "Params - ", self.design_space["['p0', 'g0']"].params
    	sim_out = odeint(self.michelis_menten,init_species,t_out, args = (mm_params))
    	g0_res, g1_res = sim_out.T
    	print g0_res
    	#create experimental data by making true sim noisy
    	exp_g0_res = g0_res + np.random.randn(len(g0_res))*0.5
    	print "args that im screwing up somehow ", exp_g0_res
    	print "type = ", type(exp_g0_res)
    	# # Find the maximum likelihood value.
    	xdata = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24]
    	mychi2 = lambda *args: -2 * self.mylnlike(*args)
    	result = op.minimize(mychi2, [l00_p, k00_p], args=(xdata,exp_g0_res))
    	# Set up the sampler.
    	ndim, nwalkers = 2, 10
    	pos = [result["x"] + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
    	sampler = emcee.EnsembleSampler(nwalkers, ndim, self.mylnprob, args=(xdata,exp_g0_res))
    	# Clear and run the production chain.
    	print("Running MCMC...")
    	sampler.run_mcmc(pos, 10, rstate0=np.random.get_state())
    	print("Done.")
    	
    	return cur_state

    def notmy_michelis_menten(self, y, t, *args):
    	Vmax = args[0]
    	km = args[1]
    	St = args[2]
    	P = y[0]
    	S = St-P

    	dP = Vmax * (S / (S+km))
    	return dP

    def execute_policy(self, policy, cur_state, cur_state_name):
    	#print cur_state
    	#action = policy[cur_state_name]
    	action = "['p0', 'g0']"
    	np.random.seed(123)
    	print "CURVE FIT TEST:"
    	xdata = np.linspace(0, 4, 50)
    	y = func(xdata, 2.5, 1.3, 0.5)
    	ydata = y + 0.2 * np.random.normal(size=len(xdata))
    	popt, pcov = curve_fit(func, xdata, ydata)
    	print popt, pcov

    	print "\nCURVE FIT ODE TEST"
    	xdata = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24]
    	global sub_model
    	sub_model = self.design_space[action].model.subs(self.testdata)
    	init_species = [10,1]
    	t_out = np.arange(0,25,1)
    	sim_out = odeint(odeint_model,init_species,t_out)
    	#print "SIMULATION of design ", des, ":"
    	g0_res, g1_res = sim_out.T
    	print g0_res
    	print "is this the same? ", sim_out.T[0]
    	#y = g0_res
    	# ydata = y + 0.2 * np.random.normal(size=len(xdata))
    	# popt, pcov = curve_fit(odeint(odeint_model,init_species,t_out).T[0], xdata, ydata)
    	# print popt, pcov
    	
    	print "\nNEW ODE WITH PARAMS AS ARGS"
    	#sub_model = self.design_space[action].model.subs(self.testdata)
    	init_species = [10,1]
    	t_out = np.arange(0,25,1)
    	l00_p = 2.256
    	k00_p = 5.658
    	mm_params = (l00_p, k00_p)
    	print "Params - ", self.design_space["['p0', 'g0']"].params
    	sim_out = odeint(self.michelis_menten,init_species,t_out, args = (mm_params))
    	#print "SIMULATION of design ", des, ":"
    	g0_res, g1_res = sim_out.T
    	print g0_res
    	#create experimental data by making true sim noisy
    	exp_g0_res = g0_res + np.random.randn(len(g0_res))*0.5
    	#ydata = g0_res + np.random.randn(len(g0_res))*0.5
    	#exp_g0_res = exp_g0_res[::5]
    	
    	def residuals(p):
    		p = tuple(p)
    		print "Params chosen = ", p
    		print "Simulating ..."
    		sim_P = odeint(self.michelis_menten, init_species,t_out, args=p) #.flatten()
    		res = sim_P.T[0] - exp_g0_res
    		print "sim res = ", sim_P.T[0]
    		print "chardata = ", exp_g0_res
    		print "residual results = ", res
    		#res = distance.euclidean(np.array(sim_P.T[0]), np.array(exp_g0_res))
    		return res
    	initial_guess = [2.0,5.0]
    	fitted_params = op.minimize(residuals, initial_guess, method = 'SLSQP', bounds = [(0,10),(0,10)])[0]
    	#fitted_params = leastsq(residuals, initial_guess)[0]
    	#x_3d = np.array([[xdata],[ydata]])
    	#popt, pcov = curve_fit(np.array(odeint(self.michelis_menten, init_species, t_out).T[0]), x_3d, x_3d[1,:], initial_guess)
    	print "Fitted Params = ", fitted_params


    	# print "Example ODE fitting problem:"
    	# Vmax = 1
    	# km = 3
    	# St = 10
    	# notmymm_params = (Vmax, km, St)
    	# #inits
    	# P_0 = 0
    	# n_steps = 100
    	# t = np.linspace(0,50,n_steps)

    	# num_P = odeint(self.notmy_michelis_menten, P_0, t, args = (notmymm_params))

    	# #create exp data
    	# exp_P = num_P + np.random.randn(len(num_P)) * 0.5
    	# print "This is expP ", exp_P
    	# #exp_P = exp_P[::5]
    	# #print "This is expP[::5] ", exp_P

    	# def residuals2(p):
    	# 	p = tuple(p)
    	# 	sim_P = odeint(self.notmy_michelis_menten, P_0, t, args = p).flatten()
    	# 	#res = sim_P[::5] - exp_P
    	# 	res = sim_P - exp_P
    	# 	return res.flatten()

    	# initial_g = [1,2,5]
    	# fitted_P = leastsq(residuals2, initial_g)[0]
    	# print "FITTETED: ", fitted_P
    	return cur_state

    def execute_fake_policy(self, policy, cur_state_name):
    	action = policy[cur_state_name]
    	print "\n\n========================================="
    	print "BUILDING/TESTING DESIGN: ", action
    	print "========================================="
    	#assume build/test/and fit just succeeds ... because i couldnt get the param fitting to work
    	cur_design_params = self.design_space[action].params
    	print cur_design_params
    	for param in cur_design_params:
    		pvalue = self.testdata[str(param)]
    		low = pvalue - random.uniform(0,2)
    		high = pvalue + random.uniform(0,2)
    		if low < 0:
    			low = 0
    		if high < 0:
    			high = 0
    		self.param_space.update({param: [low, high]})
    		#print "new param_space = ", self.param_space
    	next_state = self.belugaMDP[cur_state_name].actions[action][0][0]
    	#for t in self.belugaMDP[cur_state].actions[action]:
    	return next_state
    def updatePolicy(self, state, policy):
    	new_policy = {}
    	Q = util.Queue()
    	if len(self.belugaMDP[state].actions) != 0:
    		Q.push(state)

    	while Q.isEmpty()==0:
    		s = Q.pop()
    		print "in update policy: ", self.belugaMDP[s].actions 
    		for a, T in self.belugaMDP[s].actions.iteritems():
    			if len(self.belugaMDP[T[0][0]].actions) != 0:
    				Q.push(T[0][0])
    		new_policy.update({s: policy[s]})

    	return new_policy

    def bayes_model_inf(self, action, fit_data):
    	# for par in self.design_space[action].params:
    	# 	if str(par) != "L11" and str(par)!="K11":
    	# 		self.param_space.update({par: [0, 0.01]})
    	# #print "Param Space before optimizing: ", self.param_space
    	# # self.goal_dist(action, 3)
    	# for par in self.design_space[action].params:
    	# 	if str(par) == "L11":
    	# 		self.param_space.update({par: [8.0,10.0]})
    	# 		#self.param_space.update({par: [min(l_11_samples), max(l_11_samples)]})
    	# 	elif str(par) == "K11":
    	# 		self.param_space.update({par: [7.0,10.0]})
    	# 		#self.param_space.update({par: [min(k_11_samples), max(k_11_samples)]})
    	# 	else :
    	# 		self.param_space.update({par: [0.0 , 0.0001]})
    	# print "Param Space after optimizing: ", self.param_space

    	# self.goal_dist(action, 6)
    	#####generating char data for p0g0#######
    	global sub_model
    	# init_species = [10,1]
    	# #action = "['p0', 'g0']"
    	# tout = np.arange(0,25,1)
    	submat = {}
    	
    	# for par in self.design_space[action].params:
    	# 	submat.update({str(par): self.testdata[str(par)]})
    	# # submat.update({"L00" : 0})
    	# # submat.update({"K00" : 0})
    	# sub_model = self.design_space[action].model.subs(submat)
    	# sim_out = odeint(odeint_model,init_species,tout)
    	# g0_res, g1_res = sim_out.T
    	# #exp_g0_res = g0_res + np.random.randn(len(g0_res))*0.5
    	# count_data = np.array(g0_res)
    	# n = 2
    	# color=iter(cm.rainbow(np.linspace(0,1,n)))
    	# c=next(color)
    	# plt.plot(x,y,c=c)
    	
    	#count_data = np.array(self.goal[2])
    	count_data = np.array(fit_data)
    	print "submat is ", submat
    	# print "model is ", sub_model
    	print "Fitting to Data = ", count_data
    	
    	n_count_data = len(count_data)
    	# joke_arr = [1,2,3,4,5,4,3,2,1,2,3,4,5,4,3,2,1,2,3,4,5,4,3,2]
    	# joke_data = np.array(joke_arr)
    	# c=next(color)
    	
    	#pl.plot(np.arange(n_count_data), count_data, 'r--', np.arange(len(joke_data)), joke_data, 'b')
    	#pl.plot(np.arange(n_count_data), [count_data, joke_data], color="#348ABD")
    	# c=next(color)
    	# pl.plot(np.arange(len(joke_data)), joke_data,c=c)
    	
    	pl.xlabel("Time")
    	pl.ylabel("Fluoresence (SP)")
    	pl.title("Goal Data")
    	pl.xlim(0, n_count_data);
    	#pl.savefig("Char_data.png")
    	#########################################
    	
    	parhash = {}
    	for par in self.design_space[action].params:
    		parhash.update({str(par) : pm.Uniform(str(par),self.param_space[par][0],self.param_space[par][1], 
    			value=self.param_space[par][0])})

    	print "parhash is ", parhash
    	# l_00 = pm.Uniform('l_00', 0., 0.5, value=0.)
    	# k_00 = pm.Uniform('k_00', 0., 0.5, value=0.)
    	# l_01 = pm.Uniform('l_01', 0., 0.5, value=0.)
    	# k_10 = pm.Uniform('k_10', 0., 0.5, value=0.)
    	# l_10 = pm.Uniform('l_10', 0., 0.5, value=0.)
    	# k_01 = pm.Uniform('k_01', 0., 0.5, value=0.)
    	# l_11 = pm.Uniform('l_11', 0., 10., value=0.)
    	# k_11 = pm.Uniform('k_11', 0., 10., value=0.)
    	SI_0 = pm.Uninformative('SI_0', value=[10,1])

    	simulated_trajectories = []

    	# testhash = {}
    	# testhash.update({'l_11':pm.Uniform('l_11', 8., 10., value=9.)})
    	# testhash.update({'k_11':pm.Uniform('k_11', 7., 10., value=9.)})
    	
    	# param_arr = []
    	# for par in self.design_space[action].params:
    	# 	param_arr.append(pm.Uniform(str(par), 0., 10., value=0.))
    	# print "Param array = ", param_arr
    	# #print "Random output:", l_00.random(), l_00.random(), k_00.random()
    	# print "Random output:", param_arr[0].random(), param_arr[1].random()
    	
    	@pm.deterministic
    	#def SI(SI_0=SI_0, l_00=l_00, k_00=k_00, l_01=l_01, k_10=k_10, l_10=l_10, k_01=k_01, l_11=l_11, k_11=k_11):
    	#def SI(SI_0=SI_0, l_11=l_11, k_11=k_11):
    	def SI(SI_0=SI_0, parhash = parhash):
    	# def SI(SI_0=SI_0, parslist=(0,0)):
    		global sub_model
    		init_species = [10,1]
    		#action = "['p0', 'g1', 'p1', 'g0']"
    		t_out = np.arange(0,25,1)
    		submat = {}
    		# for key, elem in parhash.iteritems():
    		# 	submat.update({str(key): key})
    		for prm in self.design_space[action].params:
    			submat.update({str(prm) : parhash[str(prm)]})
    		# submat.update({"L00" : 0})
    		# submat.update({"K00" : 0})
    		# submat.update({"L01" : 0})
    		# submat.update({"K10" : 0})
    		# submat.update({"L10" : 0})
    		# submat.update({"K01" : 0})
    		# submat.update({"L11" : testhash['l_11']})
    		# submat.update({"K11" : testhash['k_11']})
    		print "submat in SI fcn : ", submat
    	
    		sub_model = self.design_space[action].model.subs(submat)
    		with stdout_redirected():
    			sim_out = odeint(odeint_model,init_species,t_out)
    			g0_res, g1_res = sim_out.T
    		simulated_trajectories.append(g0_res)
    		############
    		#ok we should check the goal_dist here
    		#if it's <= current epsilon then we need to save these param values 
    		#and use them to readjust the priors we'll sample from in the future

    		return sim_out.T

    	# def getpars():
    	# 	return 

    	S = pm.Lambda('S', lambda SI=SI: SI[0])
    	#simparams = pm.Lambda('simparams', lambda getpars=getpars: getpars)

    	observation = pm.Poisson("obs", mu = S, value=count_data, observed=True)
    	parhash.update({'observation': observation})
    	#model = pm.Model([observation, l_11, k_11])
    	#model = pm.Model([observation, l_00, k_00])
    	#model = pm.Model([observation, param_arr[0], param_arr[1]])
    	#testhash.update({'observation': observation})
    	
    	mcmc = pm.MCMC(parhash)
    	#mcmc = pm.MCMC(model)
    	# mcmc = pm.MCMC(set([observation, l_00, k_00, l_01, k_10, l_10, k_01, l_11, k_11]))
    	# mcmc = pm.MCMC(set([observation, l_11, k_11]))

    	###########
    	#ok we'll need to reiterate over this 
    	mcmc.sample(200, 100, 1)
    	
    	# l_00_samples = mcmc.trace('l_00')[:]
    	# k_00_samples = mcmc.trace('k_00')[:]
    	# l_01_samples = mcmc.trace('l_01')[:]
    	# k_10_samples = mcmc.trace('k_10')[:]
    	# l_10_samples = mcmc.trace('l_10')[:]
    	# k_01_samples = mcmc.trace('k_01')[:]
    	# l_11_samples = mcmc.trace('l_11')[:]
    	# k_11_samples = mcmc.trace('k_11')[:]
    	# print "\nL00 SAMPLES: ", l_00_samples
    	# print "\nK00 SAMPLES: ", k_00_samples
    	# print "\nL01 SAMPLES: ", l_01_samples
    	# print "\nK10 SAMPLES: ", k_10_samples
    	# print "\nL10 SAMPLES: ", l_10_samples
    	# print "\nK01 SAMPLES: ", k_01_samples
    	for pr in self.design_space[action].params:
    		print "\n ", str(pr), " SAMPLES: ", mcmc.trace(str(pr))[:]
    	# print "\nL11 SAMPLES: ", l_11_samples
    	# print "\nK11 SAMPLES: ", k_11_samples

    	n = len(simulated_trajectories)
    	color=iter(cm.rainbow(np.linspace(0,1,n)))
    	for i in range(n):
    		c = next(color)
    		pl.plot(np.arange(n_count_data),simulated_trajectories[i], c=c)
    	pl.plot(np.arange(n_count_data), count_data,'b--')
    	pl.savefig("sim_traj.png")

    	#ax = pl.subplot(310)
    	# P.hist(mcmc.trace('l_00')[:], normed= 1)
    	# P.savefig('hist_l00')
    	# pl.hist(mcmc.trace('l_11')[:], normed= 1)
    	# pl.hist(mcmc.trace('k_11')[:], normed= 1)
    	#pl.savefig('hist_l11')

    	############### Now let's update the current knowledge with these "learned params" #######
    	####### Even though these aren't truly learned because they're from goal fitting not char fitting ####
    	#### And then let's run goal dist and see if it gives us better results ####
    	# for par in self.design_space[action].params:
    	# 	if str(par) != "L11" and str(par)!="K11":
    	# 		self.param_space.update({par: [0, 0.5]})
    	# print "Param Space before optimizing: ", self.param_space
    	# self.goal_dist(action, 3)
    	# for par in self.design_space[action].params:
    	# 	if str(par) == "L11":
    	# 		self.param_space.update({par: [min(l_11_samples), max(l_11_samples)]})
    	# 	if str(par) == "K11":
    	# 		self.param_space.update({par: [min(k_11_samples), max(k_11_samples)]})
    	# print "Param Spacezzzzzzz after optimizing: ", self.param_space
    	# self.goal_dist(action, 3)

    	return 0

    def search(self):
    	global sims_total
    	print "IN SEARCH - here's the mdp:", self.belugaMDP
    	#cur_state = self.belugaMDP['S0']
    	cur_state_name = 'S0'
    	policy_dict = self.initPolicy() #... for each state in the full MDP state space, choose the simplest action
    	#print "Init policy is: ", policy_dict
    	alg_round = 1
    	self.bayes_model_inf("['p0', 'g0']", self.goal[2])
    	#while !cur_state.isTerminal(): ... while terminal state not experimentally reached:
    	while len(self.belugaMDP[cur_state_name].actions) > 0:
	    	policy_dict = self.findPolicy(policy_dict, alg_round) #... find optimal policy
	    	print "current optimal policy ", policy_dict
	    	#print "\ntest data = ", self.testdata
	    	#print "\nparam space = ", self.param_space
	    			#evaluate current policy
	    				#compute values (using transition and reward) for each state following current policy
	    			#improve policy
	    	next_state_name = self.execute_fake_policy(policy_dict, cur_state_name) #... do one step of current policy
	    	####X -update learned global parameters for next round of sims (should be done in execute policy)
	    	####X -update rewards
	    	cur_state_name = next_state_name #... make current step new start node
	    	print "Next state = ", cur_state_name
	    	policy_dict = self.updatePolicy(cur_state_name, policy_dict)
	    	print "Policy is = ", policy_dict
	    	alg_round = alg_round + 1
	    	print "Total number of sims = ", sims_total
    	#print "FINAL DESIGN: ", action, self.MDP_rewards[action]
    	#for state in self.belugaMDP:
    	#print "a failed design: ['p0', 'g1', 'p1', 'g0'] ", self.MDP_rewards["['p0', 'g1', 'p1', 'g0']"]
    	return 0