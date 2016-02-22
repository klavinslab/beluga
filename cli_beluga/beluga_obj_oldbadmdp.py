import numpy as np
from sympy import *
import os, sys
import errno
import util
from scipy.integrate import odeint

sub_model = []
project_path = "/Users/Leli/beluga_cli_test_dir/ABCSMC"

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
	temp_g0, temp_g1, temp_g2 = x_vec
	new_values = list(sub_model.subs({G0: temp_g0, G1: temp_g1, G2: temp_g2}))
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

		#print "DESIGNS: ", cur_des_id, test_des_id

		for i in range(len(test_des_id)):
			#print "i = ", i
			if test_des_id[i] == cur_des_id[cur_des_idx]:
				#print "elements equal: ", test_des_id[i], cur_des_id[cur_des_idx]
				cur_des_idx+=1
			else:
				cur_des_idx = 0
				if test_des_id[i] == cur_des_id[cur_des_idx]:
					#print "elements equal: ", test_des_id[i], cur_des_id[cur_des_idx]
					cur_des_idx+=1
		if cur_des_idx == len(cur_des_id):
			result = True
		else:
			result = False

		return result

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
	    		rel_prom_idx = self.id.index(gene) - 1
	    		prom = self.id[rel_prom_idx]
	    		prom_idx = promoters.index(prom)
	    		gene_idx = species.index(gene)
	    		mat[gene_idx][prom_idx] = 1

	    #print "matrix for ", self.id, " = ", mat
	    sum_mat = sum_cols(mat)

	    for n in species:
	    	odes.append(0)
	    
	    #GENERATE GENERIC ODE MODEL FROM DESIGN MATRIX
	    for a in range(0, len(species)):
	        for j in range(0, len(promoters)):
	            if (mat[a][j]!=0):
	                for i in range(0, len(species)):
	                    if (sum_mat[[i]]>0):
	                        odes[a]+= (Ls[j][i]) / (1+Gs[0][i]*Ks[i][j])
	                        Lparam_exists = Ls[j][i] in params
	                        Kparam_exists = Ks[i][j] in params
	                        if Lparam_exists == False:
	                            params.append(Ls[j][i])
	                        if Kparam_exists == False:
	                            params.append(Ks[i][j])

	        odes[a]+=-Gs[0][a]
	 
	    model = Matrix([odes[0]]+[odes[1]]+[odes[2]])


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
		init_species = [1,1,1]
		t_out = np.arange(0,25,1)
		sim_out = odeint(odeint_model,init_species,t_out)
		#print "SIMULATION of design ", des, ":"
		g0_res, g1_res, g2_res = sim_out.T
		print list(g0_res)

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

	#Find all designs that des is a subdesign of
	for key, subdes in design_dict.iteritems():
		if des.isSubDesign(subdes):
			if des.id != subdes.id:
				sub_des_list.append((subdes,len(subdes.id)))

	#Finding successors: choose min subdes's
	#MAY CHANGE THIS TO ACCEPT ANY DESIGN WHICH IS A SUPERSET OF THIS DES AND NOT JUST THE NEXT CLOSEST ONE
	if len(sub_des_list) > 0:
		min_des_len = min(sub_des_list,key=lambda item:item[1])[1]
	for item in sub_des_list:
		if item[1] == min_des_len:
			successors.append(str(item[0].id))

	#Now find all designs that differ only by one element from des
	for key, otherdes in design_dict.iteritems():
		if des!=otherdes and des.isSimilar(otherdes) and (str(otherdes.id) not in successors):
			successors.append(str(otherdes.id))

	return successors

class MDP_node:
	"""
	A class for the MDP formulation of our search problem
	"""

	def __init__(self, design_indx, design_space):
		self.my_design = design_indx
		self.symb_param = {} #self.getSymbParams(design_space)
		self.actions = {} #self.getActions(design_space)
		self.successors = []

	def getSymbParams(self, design_space):
		return []

	def getActions(self, design_space):
		return []

class beluga_obj:
    """
    A class to setup and run the beluga algorithm.
    
    """ 
    # instantiation
    def __init__(self, language, design_obj):
    	#self.species = ['g0','g1','g2']
    	self.species = design_obj.species
    	self.promoters = ['p0','p1','p2']
    	self.goal = design_obj.output
    	self.testdata = design_obj.testdata
    	self.Ls = []
    	self.Ks = []
    	self.Gs = []
    	self.param_space = {}
    	self.design_space = self.genGraph(language)
    	self.belugaMDP = self.initMDP()
    	self.policy = self.initPolicy()
    	
    	

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

    		#Finding successors:
    		successors = getSuccessors(des,design_dict)
    		for item in successors:
    			if item not in des.successors:
    				des.successors.append(item)
    			#updating the successor graph edges in the other direction
    			# if des.id not in design_dict[item].successors:
    			# 	design_dict[item].successors.append(des.id)

    	for key, thing in design_dict.iteritems():
    		print " key = ", key, type(key)
    		print "\n id = ", thing.id, type(thing.id)
    		print "\n params = ", thing.params
    		print "\n successors = ", thing.successors
    		print "-----------------------------------\n\n"

    	print "Num Designs in Space = ", len(design_dict)
    	print "Param space: ", self.param_space
    	return design_dict

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
    	start_node.symb_param = self.param_space #copy all keys in global param list and put a 1 value next to it
    	for p_key, p in start_node.symb_param.iteritems():
    		start_node.symb_param.update({p_key: 1})
    	print "Start node symb params = ", start_node.symb_param
    	start_node.actions = {"['p0', 'g0']": [], "['p1', 'g0']": [], "['p2', 'g0']": []}
    	start_element = ("S"+str(s_key), start_node)
    	mdp_Q.push(start_element)
    	#mdp_dict.update({"S"+str(s_key): start_node})
    	#only update to mdp_dict when state is popped from mdp_Q
    	#for key, mdpnode in mdp_dict.iteritems():
    	while mdp_Q.isEmpty()==0:
    		print "Q at start of while: "
    		#mdp_Q.printMDPQ()
    		mdpelement = mdp_Q.pop()
    		key, mdpnode = mdpelement[0], mdpelement[1]
    		mdp_dict.update({key: mdpnode})
    		print "MDP NODE = ", key
    		#print "Symbparam for state S",str(key), " = ", mdpnode.symb_param
    		
    		for action, s_prime in mdpnode.actions.iteritems():
    			#determine how many s_prime result states you get from this action
    				#gotta figure this part out
    				#print key, " action = ", action
    				a_des_params = self.design_space[action].params
    				exp_des_params = []

    				num_s_primes = 0
    				for param in a_des_params:
    					#print param
    					#print mdpnode.symb_param[param]
    					if mdpnode.symb_param[param] == 1:
    						num_s_primes += 1
    						exp_des_params.append(param)
    				num_s_primes = 2**num_s_primes
    				#print "Number of possible result states for action ", action, " = ", num_s_primes
    				#num_s_primes is the number of result states
    				#exp_des_params should hold the array of param symbols of which we'll expand
    				#we need to now expand
    				new_states_list = []
    				for i in range(num_s_primes-1):
    					s_key += 1
    					new_symb_param = dict(mdpnode.symb_param)
    					#take the binary representation of i
    					#for the length of exp_des_params
    					for exp_p in exp_des_params:
    						new_symb_param[exp_p] = (i >> exp_des_params.index(exp_p)) & 1

    					#new_symb_param should have the symb_param of 1 new state
    					new_mdpnode = MDP_node(action, self.design_space)
    					new_mdpnode.symb_param = new_symb_param
    					#print "Symbparam for state S",str(s_key), " = ", new_mdpnode.symb_param
    					### FIX ME ###
    					#At this point the symb_param for the new node is correct. 
    					#But when we take the node off the Q later ... it takes the orignal value of symb_param
    					#create the actions for this new mdpnode
    					#first check if it's a terminal node (i.e most params have been learned)
    					sym_param_sum = 0
    					for v_key, val in new_mdpnode.symb_param.iteritems():
    						sym_param_sum += val
    					perct_unknown = float(float(sym_param_sum)/float(len(new_symb_param)))
    					#print "Percent unknown in state S",str(s_key), " with action ", action, " = ", perct_unknown
    					if perct_unknown >= 0.3 :
	    					for succ in self.design_space[action].successors:
	    						new_mdpnode.actions.update({succ: []})
    					#go back and update the sPrime state for this action in the parent mdpnode
    					mdpnode.actions[action].append("S"+str(s_key))
    					new_element = ("S"+str(s_key), new_mdpnode)
    					new_states_list.append(new_element)
    					#print "Appending ", new_element[0], new_element[1].symb_param
    					mdp_Q.push(new_element)
    					# print "Q at end of while: "
    					# mdp_Q.printMDPQ()

    				#print "Q after 1 action expansion: "
    				#mdp_Q.printMDPQ()
    					
    		# print "actions = ", mdpnode.actions
    		# print "symb params = ", mdpnode.symb_param

    			#for each s_prime result state
    			#increment s_key
    			#create a new MDP_node for each s_prime
    			#append the current s_prime result state s_key to the parent mdpnode's 
    			#result state array - associated with this action
    			## HOPEFULLY THIS TERMINATES :/ ##
    			## DONT THINK IT DOES ##
    	print "Num MDP states = ", len(mdp_dict)

    	return mdp_dict

    def initPolicy(self):
    	#start with MDP_graph[0] ... i.e. S0
    	#for every state in MDP_graph, determine leftmost action it can take (indx of first successor of design in des_space)
    	#update self.policy with {s:des_indx}
    	return 0

    def getTransition(self):
    	#these are computed based on the current state and its successors
    	#we can initially assume something of a uniform distribution
    	return 0

    def getReward(self):
    	#this is computed by the current state's design
    	#specifically it's parameters as well as it's MC distance from the goal behavior
    	#(using current param space knowledge)
    	return 0

    def search(self):
    	#getStartNode()
    	#initPolicy()

    	#while terminal state not experimentally reached:
    		#find optimal policy
	    		#evaluate current policy
	    			#compute values (using transition and reward) for each state following current policy
	    		#improve policy
	    	#do one step of current policy
	    		#update learned global parameters for next round of sims
	    		#make current step new start node

    	return 0


    # def getDesignError(self, des_node, goal_beh):
    # 	#simulate the cur_des using it's des_node[0].model and using des_node[1] knowledge
    # 	#and determine the probability that it generated the goal behavior?  
    # 	return 1

    # def getStartState(self):
    # 	start_des = self.design_space[0]
    # 	start_knowledge = self.param_space
    # 	return (start_des,start_knowledge)

    # def isGoalState(self, des_node):
    # 	#90% of all params narrowed down? (we could eliminate params who are tied to known 0 params)
    # 	#So if 10% unknown
    # 	#OR
    # 	#design with 10% error
    # 	isgoal = False
    	
    # 	cur_des = des_node[0]
    # 	cur_param_know = des_node[1]
    # 	param_range_sum = 0
    # 	total_unknown_sum = 0

    # 	#Do parameter inference for this design and update its cur_param_know
    # 	cur_param_know = cur_des.testDesign(cur_param_know,self.testdata)

    # 	for key in cur_param_know:
    # 		param_range_sum += cur_param_know[key][1] - cur_param_know[key][0]
    # 		total_unknown_sum += 10

    # 	perc_unknown = param_range_sum / total_unknown_sum

    # 	if self.getDesignError(des_node,self.goal) < 0.10:
    # 		print "We found a design that works!"
    # 		return True

    # 	if perc_unknown < 0.10:
    # 		print "We know enough to say no such design exists in this space."
    # 		return True

    # 	return False

    # def getSearchSuccessors(self, node):
    # 	#returns successor states, actions they require, and a cost of 1
    # 	#returns all successor designs in node[0].successors and appends the current knowledge state to it and cost
    # 	return 0

    # def heuristic(self, des_node):
    # 	#distance is just how far we are from knowing all params (reducing uncertainy)
    # 	return 0

    # def search(self):
	   #  #fringe = util.PriorityQueue()
	   #  fringe = util.Queue()
	   #  start_state = self.getStartState()
	   #  #fringe.push([start_state],heuristic(start_state))
	   #  fringe.push([start_state])
	    
	   #  visited = []
	    
	   #  while fringe.isEmpty()==0:
	   #      path = fringe.pop()
	   #      node = path[-1]
	   #      if node != start_state:
	   #          node = node[0]
	        
	   #      visited.append(node)
	        
	   #      if self.isGoalState(node):
	   #          actions = getPathActions(path, problem)
	   #          return actions

	   #      for i in self.getSearchSuccessors(node):
	   #          if(i[0] not in visited):
	   #              new_path = list(path)
	   #              new_path.append(i)
	                
	   #              actions = getPathActions(new_path,problem)
	   #              fringe.push(new_path)
	   #              #fringe.push(new_path,problem.getCostOfActions(actions)+heuristic(node,problem))
	                
	   #  print "somehow fringe is empty?"
	   #  return 0
