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


def getModel(cur_design):
	return([],[])

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

class beluga_obj:
    """
    A class to setup and run the beluga algorithm.
    
    """ 
    # instantiation
    def __init__(self, language):
    	self.design_space = self.genGraph(language)

    def genGraph(self, language):
    	design_list = []
    	for elem in language:
    		design_list.append(design(elem))
    	print "DESIGN: "


    	#finding all successors in one direction for each design in design list
    	for des in design_list:

    		#Generating ODE model and params:
    		des.model,des.params = getModel(des)

    		#Finding successors:
    		successors = getSuccessors(des,design_list)
    		for item in successors:
    			if item not in des.successors:
    				des.successors.append(item)
    			if des not in item.successors:
    				item.successors.append(des)

    		#updating the successor graph edges in the other direction
    		# for elem in des.successors:
    		# 	elem.successors.append(des)


    	for thing in design_list:
    		print thing.id, " successors = \n"
    		for succ in thing.successors:
    			print succ.id, "\n"


    	return design_list

    def search():
        return 0
