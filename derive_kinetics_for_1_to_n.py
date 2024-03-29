import sys
import time
import argparse
from itertools import combinations

class Reaction(object):
    _number=0
    _parent=[]
    _kd_number=0
    _child=[]
    def __init__(self, number, parent, kd_number, child):
        self._number=number
        self._parent=parent
        self._kd_number=kd_number
        self._child=child
    def to_string(self):
        #r1 = -a*x + kd1*a1x
        return f"r{self._number} = -{self._child}*x + kd{self._kd_number}*{self._parent}"
    def is_in(self, s):
        if s==self._child:
            return 1
        if s==self._parent:
            return -1
        return None

def find_indexes_and_sides_of_reactions_involving_symbol(reactions, symbol):
    special_positive_negative_indexes=[]
    for r in reactions:
        idx=r.is_in(symbol)
        if idx is not None:
            special_positive_negative_indexes.append(int(r._number)*idx)
    return special_positive_negative_indexes

def list_to_comma_separated_string(l):
    s=""
    for e in l:
        s+=e+", "
    return s[:-2]

def set_to_ax_symbol(tu):
    symbol="a"
    for e in tu:
        symbol+=(str(e))+"_"
    symbol+="x"
    if symbol=="ax":
        return "a"
    else:
        return symbol

def ax_symbol_to_set(symbol):
    s=set()
    symbol=symbol.replace("a","").replace("x", "")
    [s.add(int(i)) for i in symbol.split("_") if i is not ""] 
    return s
    

    

def get_parents(tree, child):
    parent_level=len(child)-1
    child=set(child)
    parents=[]
    def is_parent(parent, child):
        if len(list(parent-child))>0:
            return False
        if len(list(child-parent))==1:
            return True
        return False

    for parent in tree[parent_level]:
        if is_parent(parent, child):
            parents.append(parent)
    return parents
    
    


def generate_1_to_N_kinetic_scheme(n_sites:int, function_name:str, out_file_name:str):
    all_sites=[n for n in range(1, n_sites+1)]
    tree=[]
    for i in range(n_sites+1):
        tree.append([set(l) for l in combinations(all_sites, i)])
    print("\n\nThe tree looks like this\n--------------")
    symbols_list=[]
    for level in tree:
        for leaf in level:
            print(set_to_ax_symbol(leaf), end=",")
            symbols_list.append(set_to_ax_symbol(leaf))
        print()
    print("--------------")
    print(tree)

    output="import numpy as np\nfrom scipy.integrate import solve_ivp\n"

    output+=f"def {function_name}(a, x, {list_to_comma_separated_string(['kd'+str(i+1) for i in range(n_sites)])}, interval=(0,1)):\n"
    tabs="\t\t"
    output+=f"\tdef ode_{function_name}(concs, t, {list_to_comma_separated_string(['kd'+str(i+1) for i in range(n_sites)])}):\n{tabs}"
    
    # Add symbols
    symbols_list_with_x=symbols_list.copy()
    symbols_list_with_x.insert(1, "x")
    output+=f"{list_to_comma_separated_string(symbols_list_with_x)} = concs\n"
    print("Symbols list", symbols_list)
    reactions=[]
    
    for symbol in symbols_list:
        if symbol=='a':
            continue
        parents=get_parents(tree, ax_symbol_to_set(symbol))
        for parent in parents:
            reactions.append(Reaction(len(reactions)+1, symbol, str(list(ax_symbol_to_set(symbol)-parent)[0]), set_to_ax_symbol(parent)))

    for r in reactions:
        output+=tabs+r.to_string()+"\n"
    
    for symbol in symbols_list:

        output+=f"{tabs}d{symbol}dt = 0.0 "
        for rid in find_indexes_and_sides_of_reactions_involving_symbol(reactions, symbol):
            if rid<0:
                output+=f"-r{rid*-1} "
            else:
                output+=f"+r{rid} "
        output+="\n"
        
    output+=f"{tabs}dxdt = "
    for r in range(1,len(reactions)+1):
        output+=f"r{r}+"
    output=output[:-1]+"\n"
    
        
    
    
    
    output+=tabs+f"return [{list_to_comma_separated_string(['d'+s+'dt' for s in symbols_list_with_x])}]\n"

    tabs=tabs[:-1]
    output+=f"{tabs}res = solve_ivp(lambda t, y: ode_{function_name}(y,t,{list_to_comma_separated_string(['kd'+str(i+1) for i in range(n_sites)])}), interval, [a, x,{list_to_comma_separated_string(['0.0' for _ in range(2, len(symbols_list)+1)])}]).y[2:,-1]\n"
   
    output+=f"{tabs}return ("

    rescounter=0
    for level_i in range(1,len(tree)-1):
        output+=f"{level_i}*(sum(res[{rescounter}:{rescounter+len(tree[level_i])}]))+"
        rescounter+=len(tree[level_i])
    if(n_sites==1):
        level_i=1
        output+="0 "
    output=output[:-1]+f"+{len(tree[level_i])}*res[{rescounter}])/x\n"
    output_file=open(out_file_name, "w")
    output_file.write(output)
    output_file.close()















if __name__ == "__main__":
    #generate_1_to_N_kinetic_scheme(int(4), "multisite_1_to_3", "/tmp/out.py")
    parser=argparse.ArgumentParser(description="Derive kinetic scheme for 1:N binding")
    parser.add_argument('n_to_pick', help="How many sites on the protein")
    parser.add_argument('func_name', help="Function name, should be something like 'one_to_four'")
    parser.add_argument('outFile', help="Output file")
    args=parser.parse_args()
    generate_1_to_N_kinetic_scheme(int(args.n_to_pick), args.func_name, args.outFile)
