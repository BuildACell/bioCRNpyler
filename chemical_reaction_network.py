from warnings import warn

class specie(object):
    def __init__(self, name, type = "", attributes = []):
        self.name = name
        self.type = type
        self.attributes = attributes

    def __repr__(self):
        if self.type != "complex":
            txt = self.type+"_"+self.name
        else:
            txt = self.name
        if len(self.attributes)>0 and self.attributes != None:
            for i in self.attributes:
                if i != None:
                    txt+="_"+str(i)
        return txt


    def __eq__(self, other):
        """Overrides the default implementation
           Two species are equivalent if they have the same name, type, and attributes"""
        if isinstance(other, specie) and repr(self) == repr(other):
            return True
        else:
            return False

    def __hash__(self):
        return str.__hash__(repr(self))

class reaction(object):
    def __init__(self, inputs, outputs, k, input_coefs = None, output_coefs = None, k_rev = 0):

        #Check that inputs and outputs only contain species
        for s in inputs+outputs:
            if not isinstance(s, specie):
                raise ValueError("A non-species object was used as a specie")
        self.inputs = list(set(inputs))
        self.outputs = list(set(outputs))

        #Check that rates are valid
        if k <= 0:
            raise ValueError("Reaction rate <= 0: k="+str(k))
        else:
            self.k = k
        if k_rev > 0:
            self.reversible = True
            self.k_r = k_rev
        else:
            self.k_r = 0
            self.reversible = False

        #Set input coefficients
        if input_coefs==None:
            self.input_coefs = [inputs.count(s) for s in self.inputs]
        elif input_coefs != None and len(input_coefs) == len(self.inputs):
            self.input_coefs = input_coefs
        elif len(input_coefs) == len(inputs) and len(self.inputs)!= len(inputs):
            raise ValueError("Input species and input_coefs contain contradictory counts.")
        else:
            raise ValueError("Len(input_coefs) does not match len(self.inputs)")

        #Set Output Coefs
        if output_coefs==None:
            self.output_coefs = [outputs.count(s) for s in self.outputs]
        elif output_coefs != None and len(output_coefs) == len(self.outputs):
            self.output_coefs = output_coefs
        elif len(output_coefs) == len(outputs) and len(self.outputs)!= len(outputs):
            raise ValueError("Output species and output_coefs contain contradictory counts.")
        else:
            raise ValueError("Len(output_coefs) does not match len(self.outputs)")

    def __repr__(self, **kwargs):
        txt = ""
        for i in range(len(self.inputs)):
            if self.input_coefs[i] > 1:
                txt += str(self.input_coefs[i])+"*"+str(self.inputs[i])
            else:
                txt += str(self.inputs[i])
            if i < len(self.inputs)-1:
                txt+=" + "
        if self.reversible:
            txt += " <--> "
        else:
            txt += " --> "
        for i in range(len(self.outputs)):
            if self.output_coefs[i]>1:
                txt += str(self.output_coefs[i])+"*"+str(self.outputs[i])
            else:
                txt += str(self.outputs[i])
            if i < len(self.outputs)-1:
                txt+= " + "
        tab = (" "*8)
        txt += tab
        if self.reversible:
            txt+="k_f="+str(self.k)+"\tk_r="+str(self.k_r)
        else:
            txt += "k_f=" + str(self.k)

        return txt

    def __eq__(self, other):
        """Overrides the default implementation
           Two reactions are equivalent if they have the same inputs, outputs, and rates"""
        equal = True

        #must both be reactions with the same rates and numbers of inputs and outputs
        if (not isinstance(other, reaction) or other.k != self.k or other.k_r != self.k_r or
            len(other.inputs) != len(self.inputs) or len(other.outputs) != len(self.outputs)):
            equal = False
        else:
            #check if the input species and coefficients are the same (order may be different)
            for i in range(len(self.inputs)):
                input = self.inputs[i]
                input_c = self.input_coefs[i]
                if (input not in other.inputs or input_c != other.input_coefs[other.inputs.index(input)]):
                    equal = False
            # check if the output species and coefficients are the same (order may be different)
            for i in range(len(self.outputs)):
                output = self.outputs[i]
                output_c = self.output_coefs[i]
                if (output not in other.outputs or output_c != other.output_coefs[other.outputs.index(output)]):
                    equal = False
        return equal

    def pyrepr(self):
        if self.reversible:
            return [
                ([repr(i) for i in self.inputs], self.input_coefs, [repr(i) for i in self.outputs], self.output_coefs, self.k),
                ([repr(i) for i in self.outputs], self.output_coefs, [repr(i) for i in self.inputs], self.input_coefs, self.k_r)
            ]
        else:
            return [([repr(i) for i in self.inputs], self.input_coefs, [repr(i) for i in self.outputs], self.output_coefs, self.k)]
class chemical_reaction_network(object):
    def __init__(self, species, reactions):

        #Check to make sure species are valid and only have a count of 1
        for s in species:
            if not isinstance(s, specie):
                raise ValueError("A non-species object was used as a specie: recieved "+repr(s))
            if species.count(s) > 1:
                warn("Species "+str(s)+" duplicated in CRN definition. Duplicates have been removed.")
            while species.count(s)>1:
                species.remove(s)
        self.species = species

        #Check to make sure reactions are valid meaning:
        #   only have a count of 1
        #   all species in the inputs/outputs are also in teh species list
        for r in reactions:
            if not isinstance(r, reaction):
                raise ValueError("A non-reaction object was used as a reaction")
            if reactions.count(r) > 1:
                warn("Reaction "+str(r)+" duplicated in CRN definitions. Duplicates have been removed.")
            while reactions.count(r) > 1:
                reactions.remove(r)

            for s in r.inputs:
                if s not in self.species:
                    warn("Reaction "+repr(r)+" contains a species "+repr(s)+" which is not in the CRN")

            for s in r.outputs:
                if s not in self.species:
                    warn("Reaction "+repr(r)+" contains a species "+repr(s)+" which is not in the CRN")

        self.reactions = reactions

    def __repr__(self):
        txt = "Species = "
        for s in self.species:
            txt+=repr(s)+", "
        txt = txt[:-2]+'\n'
        txt+="Reactions = ["+"\n"

        for r in self.reactions:
            txt+="\t"+repr(r)+"\n"
        txt+="]"
        return txt

    def add_species(self, new_species, debug = False):
        specie_added_dict = {}
        for s in new_species:
            if not isinstance(s, specie):
                raise ValueError("A non-species object was used as a specie")
            elif s in self.species:
                warn("CRN.add_species Warning: the specie " + str(s) + " is already contained in the CRN.")
                specie_added_dict[str(specie)] = False
            elif new_species.count(s)>1:
                warn("CRN.add_species Warning: the specie "+str(s)+" is listed multiple times.")
                if str(s) not in specie_added_dict:
                    self.species.append(s)
                    specie_added_dict[str(s)] = True
            else:
                self.species.append(s)
                specie_added_dict[str(s)] = True

        if debug:
            return self.species, specie_added_dict
        else:
            return self.species

    def pyrepr(self):
        reactions = []
        for r in self.reactions:
            reactions += r.pyrepr()
        species = [str(s) for s in self.species]
        return species, reactions

    def initial_condition_vector(self, init_cond_dict):
        x0 = [0.0 for s in self.species]
        for i in range(len(self.species)):
            s = str(self.species[i])
            if s in init_cond_dict:
                x0[i] = init_cond_dict[s]
        return x0

#BASIC CRN TESTS BELOW
"""
s1 = specie("g1", type="DNA")
s1double = specie("g1", type = "DNA")
s2 = specie("t1", type="mRNA")
s3 = specie("p1", type="protein", attributes = ["deg"])
s4 = specie("p1", type ="dimer")
r1 = reaction([s1], [s1, s2], k=1.0)
r2 = reaction([s2], [s2, s3], k=2.0)
r3 = reaction([s3, s3], [s4], k=100., k_rev=10.)
r4 = reaction([s3], [s4], k=100., k_rev = 10., input_coefs=[2], output_coefs=[1])
r5 = reaction([s3], [], k = .1)

print("s1 == s1double", s1 == s1double)
print("s1 == s2", s1 == s2)
print("r1 == r2", r1 == r2)
print("r3 == r4", r3 == r4)
CRN = chemical_reaction_network([s1, s1double, s2, s3, s4], [r1, r2, r3, r4, r5])
print(repr(CRN))
"""
