


class algorithm_info:
    """
    A class to parse the user-provided grammar input file and json file and return all information required to run the beluga algorithm.
    
    """ 
    
    def __init__(self, grammar_filename):
        self.grammar_file = grammar_filename

    def print_info(self):
        print "\nALGORITHM INFO"
