class Xms:

    def __init__(self, char):
        self.char = char
    def make_greek(self):
        if self.char not in ["r", "theta", "tau", "rho1", "rho2", "rho3"]:
            print "wrong char in Xmgrace_style class, exiting"; raise System_exit
        if self.char == "r"    :return r"r"            
        if self.char == "theta":return r"\f{Symbol}q"  
        if self.char == "tau"  :return r"\f{Symbol}t" 
        if self.char == "rho1" :return r"\f{Symbol}r\s1\N" 
        if self.char == "rho2" :return r"\f{Symbol}r\s2\N" 
        if self.char == "rho3" :return r"\f{Symbol}r\s3\N" 
