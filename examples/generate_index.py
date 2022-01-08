import sys, inspect
import biocrnpyler
from os import listdir
from os.path import isfile, join

# Get lists of bioCRNpyler objects of different types
species = [
    (n, o) for (n, o) in inspect.getmembers(sys.modules["biocrnpyler"]) 
    if inspect.isclass(o) and issubclass(o, biocrnpyler.Species)
]

propensities = [
    (n, o) for (n, o) in inspect.getmembers(sys.modules["biocrnpyler"]) 
    if inspect.isclass(o) and issubclass(o, biocrnpyler.Propensity)
]

components = [
    (n, o) for (n, o) in inspect.getmembers(sys.modules["biocrnpyler"]) 
    if inspect.isclass(o) and issubclass(o, biocrnpyler.Component)
]

mechanisms = [
    (n, o) for (n, o) in inspect.getmembers(sys.modules["biocrnpyler"]) 
    if inspect.isclass(o) and issubclass(o, biocrnpyler.Mechanism)
]

mixtures = [
    (n, o) for (n, o) in inspect.getmembers(sys.modules["biocrnpyler"]) 
    if inspect.isclass(o) and issubclass(o, biocrnpyler.Mixture)
]

core_objs = species+propensities+components+mechanisms+mixtures

# Find miscellanious objects
other_objs = []
for (n, o) in inspect.getmembers(sys.modules["biocrnpyler"]):
    if inspect.isclass(o) and (n, o) not in core_objs and "biocrnpyler" in str(o):
        other_objs.append((n, o))

all_objs = core_objs + other_objs

# dictionary stores the first .ipynb the object appears in
first_used = {c[0]:None for c in all_objs}

#paths to search through
paths = [".", "Specialized Tutorials"]
for path in paths:
    #find .ipynb files
    ipynb_files = [f for f in listdir(path) if isfile(join(path, f)) and f.split(".")[-1]=="ipynb"]
    for fname in ipynb_files:
        f = open(join(path, fname))
        for line in f:
            #cross references with biocrnpyler classes
            for c in first_used:
                if c in line and first_used[c] is None:
                    if path == ".":
                        first_used[c] = fname
                    else:
                        first_used[c] = path+"/"+fname
        f.close()

# Write text
txt = "BioCRNpyler Class\tObject Type\tExample Notebook\n"

# Keep track of which object names have been added to the text
written = {}

# Iterate through different object types
for n, o in species:
    if first_used[n] is not None:
        txt+=f"{n}\tSpecies\t{first_used[n]}\n"
        written[n] = True
        
for n, o in propensities:
    if first_used[n] is not None:
        txt+=f"{n}\tPropensity\t{first_used[n]}\n"
        written[n] = True
        
for n, o in components:
    if first_used[n] is not None:
        txt+=f"{n}\tComponent\t{first_used[n]}\n"
        written[n] = True
        
for n, o in mechanisms:
    if first_used[n] is not None:
        txt+=f"{n}\tMechanism\t{first_used[n]}\n"
        written[n] = True
        
for n, o in mixtures:
    if first_used[n] is not None:
        txt+=f"{n}\tMixture\t{first_used[n]}\n"
        written[n] = True
        
for n in first_used:
    if n not in written and first_used[n] is not None:
        txt+=f"{n}\tOther\t{first_used[n]}\n"
        written[n] = True

# Write the file
f = open("0. Tutorial Index.txt", 'w')
f.write(txt)
f.close()