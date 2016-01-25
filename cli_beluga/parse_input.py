import util

def parse_config(config_file):
    return 0

def getLanguage(P,T):
    Q = util.Queue()
    Q.push(['Design'])
    found = 0
    print "START"
    language = []

    while not Q.isEmpty():
        print "One while loop done"
        # Get the next sentential form
        sf = list(Q.pop())
        sf1 = [y for y in sf]
        #sf1 = list(sf[-1])
        print "sf = ", sf
        print "sf1 = ", sf1
        expand = ''

        i = 0
        j = 0
        #expand leftmost nonterminal add all children to queue
        for j in sf1:
            if j not in T:
                print "found one ", j
                expand = j
                break
            else:
                print "this a terminal: ", j
            i+=1
        if expand == '':
            language.append(sf1)

        for production in P:   
            if production[0] == expand:
                sf_copy = list(sf)
                sf_copy[i:i+1] = [x for x in production[1]]
                #sf = list(production[1])
                print "production = ", production
                Q.push(sf_copy)
                Q.printQueue()
                #i = len(sf1)
    print "LANGUAGE = ", language
    print "size = ", len(language)
    return language


class algorithm_info:
    """
    A class to parse the user-provided grammar input file and json file and return all information required to run the beluga algorithm.
    
    """ 
    
    def __init__(self, P, T):
        self.language = self.parse_grammar(P, T)
        #self.design_objective = self.parse_config(self, config_filename)

    def parse_grammar(self, P, T):
        """
        A function that returns the design space graph derived from a grammar_file

        """
        language = getLanguage(P,T)
        return language

    def print_info(self):
        print "\nALGORITHM INFO"
