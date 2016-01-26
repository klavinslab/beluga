import numpy as np
from sympy import *
from scipy.integrate import odeint

sub_model = []
#output: sums the entries in the cols of a matrix, for each row
def sum_cols(mat):
    
    sum_mat = np.zeros([len(mat)])
    
    for j in range(0,len(mat)):
        col_sum = 0
        for k in range(0,len(mat[j])):
            col_sum+=sum_mat[j]+mat[j][k]
        sum_mat[j] = col_sum       
    return sum_mat

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
		return param_sp
    	#simulate with ground_truth_vals and then
    	#do param inference on self design using it's cur_param_sp and using char_data from sim
    	#and update the cur_param_sp


def getSuccessors(des,design_list):
	successors = []
	sub_des_list = []

	#Find all designs that des is a subdesign of
	for subdes in design_list:
		if des.isSubDesign(subdes):
			#print des.id, "is a subdes of ", subdes.id
			if des.id != subdes.id:
				sub_des_list.append((subdes,len(subdes.id)))

	#Finding successors: choose min subdes's
	if len(sub_des_list) > 0:
		min_des_len = min(sub_des_list,key=lambda item:item[1])[1]
	for item in sub_des_list:
		if item[1] == min_des_len:
			successors.append(item[0])

	#Now find all designs that differ only by one element from des
	#similar_des_list = []
	for otherdes in design_list:
		if des!=otherdes and des.isSimilar(otherdes) and (otherdes not in successors):
			successors.append(otherdes)

	return successors

def odeint_model(x_vec,t):
	global sub_model
	G0, G1, G2 = symbols("G0, G1, G2")
	temp_g0, temp_g1, temp_g2 = x_vec
	new_values = list(sub_model.subs({G0: temp_g0, G1: temp_g1, G2: temp_g2}))
	return new_values

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
    	
    	

    def genGraph(self, language):
    	design_list = []
    	global sub_model
    	for elem in language:
    		design_list.append(design(elem))


    	#Generating full space of params and symbols
    	for m in range(0,len(self.species)):
			self.Ks.append(map(lambda i: Symbol('K' + str(m)+str(i)), xrange(len(self.promoters))))
			self.Ls.append(map(lambda i: Symbol('L' + str(m)+str(i)), xrange(len(self.promoters))))
			self.Gs.append(map(lambda i: Symbol('G' + str(i)), xrange(len(self.species))))
    	#finding all successors in one direction for each design in design list
    	for des in design_list:

    		#Generating ODE model and params:
    		des.getModel(self.species, self.promoters, self.Ls, self.Ks, self.Gs)
    		
    		####### TESTING ODE SIMULATION WITH SCIPY ODEINT ######
    		#######         THIS WILL NOT STAY HERE          ######
    		sub_model = des.model.subs(self.testdata)
    		#print "submodel: ", sub_model

    		init_species = [1,1,1]
    		t_out = np.arange(0,24,1)
    		sim_out = odeint(odeint_model,init_species,t_out)
    		#print "SIMULATION of design ", des, ":"
    		g0_res, g1_res, g2_res = sim_out.T
    		#print g0_res
    		des.char_data = g0_res



    		for x in des.params:
    			self.param_space.update({x: [0,10]})

    		#Finding successors:
    		successors = getSuccessors(des,design_list)
    		for item in successors:
    			if item not in des.successors:
    				des.successors.append(item)
    			#updating the successor graph edges in the other direction
    			if des not in item.successors:
    				item.successors.append(des)

    	for thing in design_list:
    		print thing.id, " params = \n"
    		print thing.params

    	print "Param space: ", self.param_space
    	return design_list

    def getDesignError(des_node, goal_beh):
    	#simulate the cur_des using it's des_node[0].model and using des_node[1] knowledge
    	#and determine the probability that it generated the goal behavior?  
    	return 1

    def getStartState(self):
    	start_des = self.design_space[0]
    	start_knowledge = self.param_space
    	return (start_des,start_knowledge)

    def isGoalState(self, des_node):
    	#90% of all params narrowed down? (we could eliminate params who are tied to known 0 params)
    	#So if 10% unknown
    	#OR
    	#design with 10% error
    	isgoal = False
    	
    	cur_des = des_node[0]
    	cur_param_know = des_node[1]
    	param_range_sum = 0
    	total_unknown_sum = 0

    	#Do parameter inference for this design and update its cur_param_know
    	cur_param_know = cur_des.testDesign(cur_param_know,self.testdata)

    	for key in cur_param_know:
    		param_range_sum += cur_param_know[key][1] - cur_param_know[key][0]
    		total_unknown_sum += 10

    	perc_unknown = param_range_sum / total_unknown_sum

    	if self.getDesignError(des_node,self.goal) < 0.10:
    		print "We found a design that works!"
    		return True

    	if perc_unknown < 0.10:
    		print "We know enough to say no such design exists in this space."
    		return True

    	return False

    def getSearchSuccessors(self, node):
    	#returns successor states, actions they require, and a cost of 1
    	#returns all successor designs in node[0].successors and appends the current knowledge state to it and cost
    	return 0

    def heuristic(self, des_node):
    	#distance is just how far we are from knowing all params (reducing uncertainy)
    	return 0

    def search(self):
	    fringe = util.PriorityQueue()
	    start_state = self.getStartState()
	    fringe.push([start_state],heuristic(start_state))
	    
	    visited = []
	    
	    while fringe.isEmpty()==0:
	        path = fringe.pop()
	        node = path[-1]
	        if node != start_state:
	            node = node[0]
	        
	        visited.append(node)
	        
	        if self.isGoalState(node):
	            actions = getPathActions(path, problem)
	            return actions

	        for i in self.getSearchSuccessors(node):
	            if(i[0] not in visited):
	                new_path = list(path)
	                new_path.append(i)
	                
	                actions = getPathActions(new_path,problem)
	                fringe.push(new_path,problem.getCostOfActions(actions)+heuristic(node,problem))
	                
	    print "somehow fringe is empty?"
	    return 0
