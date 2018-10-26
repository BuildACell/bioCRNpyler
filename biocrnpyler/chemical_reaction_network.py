from warnings import warn
import sbmlutil

#A formal species object for a CRN
#Species must have a name. They may also have a type (such as DNA, RNA, Protein), and a list of attributes

#A test comment for git

class specie(object):
    def __init__(self, name, type = "", attributes = []):
        self.name = name
        self.type = type

        if attributes == None:
            attributes = []
        elif None in attributes:
            while None in attributes:
                attributes.remove(None)

        self.attributes = attributes
    def __repr__(self):
        txt = self.type + "_" + self.name
        #if self.type != "complex":
        #    txt = self.type+"_"+self.name
        #else:
        #    txt = self.name
        if len(self.attributes)>0 and self.attributes != []:
            for i in self.attributes:
                if i != None:
                    txt+="_"+str(i)
        txt.replace("'", "")
        return txt

    #Overrides the default implementation
    #Two species are equivalent if they have the same name, type, and attributes
    def __eq__(self, other):

        if isinstance(other, specie) and self.type == other.type and self.name == other.name and set(self.attributes) == set(other.attributes):
            return True
        else:
            return False

    def __hash__(self):
        return str.__hash__(repr(self))


#An abstract representation of a chemical reaction in a CRN
#A reaction has the form:
#   \sum_i n_i I_i --> \sum_i m_i O_i @ rate = k
#   where n_i is the count of the ith input, I_i, and m_i is the count of the ith output, O_i.
#If the reaction is reversible, the reverse reaction is also included:
#   \sum_i m_i O_i  --> \sum_i n_i I_i @ rate = k_rev
class reaction(object):
    def __init__(self, inputs, outputs, k, input_coefs = None, output_coefs = None, k_rev = 0, mass_action = True, rate_formula = None):

        self.mass_action = mass_action

        if not mass_action and rate_formula == None:
            raise ValueError("A chemical reaction which is not mass action must have a rate formula.")
        elif not mass_action:
            raise NotImplementedError("rate formula have not been implemented yet.")

        self.rate_formula = rate_formula

        #Check that inputs and outputs only contain species
        for s in inputs+outputs:
            if not isinstance(s, specie):
                raise ValueError("A non-species object was used as a specie")

        self.inputs = []
        self.outputs = []
        for s in inputs:
            if s not in self.inputs:
                self.inputs.append(s)
        for s in outputs:
            if s not in self.outputs:
                self.outputs.append(s)

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
        complexes_equal = (self.complex_set_equality(self.inputs, self.input_coefs, other.inputs, other.input_coefs)
            and self.complex_set_equality(self.outputs, self.output_coefs, other.outputs, other.output_coefs))
        rates_equal = (other.k == self.k and other.k_r == self.k_r)


        #must both be reactions with the same rates and numbers of inputs and outputs
        if not isinstance(other, reaction):
            return False

        if complexes_equal and rates_equal:
            return True
        elif complexes_equal:
            warn("Two reactions with the same inputs and outputs but different rates are formally different, but may be undesired:"+repr(self)+"and "+repr(other))
            return True

        # If the reactions are reversible inverses of eachother, one's forward reaction could be the other's reverse
        elif self.reversible and other.reversible:
            reverse_complex_equal = (self.complex_set_equality(self.inputs, self.input_coefs, other.outputs, other.output_coefs)
            and self.complex_set_equality(self.outputs, self.output_coefs, other.inputs, other.input_coefs))
            reverse_rates_equal = (other.k == self.k_r and other.k_r == self.k)
            if reverse_complex_equal and reverse_rates_equal:
                return True
            elif reverse_complex_equal:
                warn(
                    "Two reversible reactions with the same inputs and outputs (reversed) but different rates are formally equal, but may be undesired:" + repr(
                        self) + "and " + repr(other))
                return True
            else:
                return False
        else:
            return False

    #Checks to see if two formal complexes (reaction input or output sets) are equal
    def complex_set_equality(self, c1, c1_coefs, c2, c2_coefs):
        if len(c1) != len(c2):
            return False
        else:
            for i in range(len(c1)):
                s1 = c1[i]
                coef1 = c1_coefs[i]
                if s1 not in c2 or coef1 != c2_coefs[c2.index(s1)]:
                    return False
        return True

    def pyrepr(self):
        if self.reversible:
            return [
                ([repr(i) for i in self.inputs], self.input_coefs, [repr(i) for i in self.outputs], self.output_coefs, self.k),
                ([repr(i) for i in self.outputs], self.output_coefs, [repr(i) for i in self.inputs], self.input_coefs, self.k_r)
            ]
        else:
            return [([repr(i) for i in self.inputs], self.input_coefs, [repr(i) for i in self.outputs], self.output_coefs, self.k)]

#A chemical reaction network is a container of species and reactions
#chemical reaction networks can be compiled into SBML or represented conveniently as python tuple objects.
#reaction types:
#   mass action: standard mass action semantics where the propensity of a reaction is given by
#           deterministic propensity = k \Prod_{inputs i} [S_i]^a_i
#           stochastic propensity = k \Prod_{inputs i} (S_i)!/(S_i - a_i)!
#           where a_i is the stochiometric coefficient of species i
class chemical_reaction_network(object):
    def __init__(self, species, reactions):
        self.species, self.reactions = self.check_crn_validity(reactions, species)

        self.species2index = {}
        for i in range(len(self.species)):
            self.species2index[str(self.species[i])] = i



    def check_crn_validity(self, reactions, species):
        # Check to make sure species are valid and only have a count of 1
        checked_species = []
        for s in species:
            if not isinstance(s, specie):
                raise ValueError("A non-species object was used as a specie: recieved "+repr(s))
            if species.count(s) > 1:
                warn("Species "+str(s)+" duplicated in CRN definition. Duplicates have been removed.")
            if s not in checked_species:
                checked_species.append(s)
        species = checked_species

        # Check to make sure reactions are valid meaning:
        #   only have a count of 1
        #   all species in the inputs/outputs are also in the species list
        checked_reactions = []
        for r in reactions:
            if not isinstance(r, reaction):
                raise ValueError("A non-reaction object was used as a reaction")

            if reactions.count(r) > 1:
                warn("Reaction " + str(r) + " duplicated in CRN definitions. Duplicates have been removed.")

            if reaction not in checked_reactions:
                checked_reactions.append(reaction)

            for s in r.inputs:
                if s not in species:
                    warn("Reaction " + repr(r) + " contains a species " + repr(s) + " which is not in the CRN")

            for s in r.outputs:
                if s not in species:
                    warn("Reaction " + repr(r) + " contains a species " + repr(s) + " which is not in the CRN")

        return species, reactions


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

    def pyrepr(self):
        reactions = []
        for r in self.reactions:
            reactions += r.pyrepr()
        species = [str(s) for s in self.species]
        return species, reactions

    def species_index(self, species):
        if len(self.species2index) != len(self.species):
            self.species2index = {}
            for i in range(len(self.species)):
                self.species2index[str(self.species[i])] = i
        return self.species2index[str(species)]

    def initial_condition_vector(self, init_cond_dict):
        x0 = [0.0 for s in self.species]
        for i in range(len(self.species)):
            s = str(self.species[i])
            if s in init_cond_dict:
                x0[i] = init_cond_dict[s]
        return x0


    def generate_sbml_model(self, stochastic_model = False, **keywords):
        document, model = sbmlutil.create_sbml_model(**keywords)

        for s in self.species:
            sbmlutil.add_species(model, model.getCompartment(0), s)

        rxn_count = 0
        for r in self.reactions:
            rxn_id = "r"+str(rxn_count)
            sbmlutil.add_reaction(model, r.inputs, r.input_coefs, r.outputs, r.output_coefs, r.k, rxn_id,
                stochastic = stochastic_model, mass_action=r.mass_action)
            rxn_count += 1
            if r.reversible:
                sbmlutil.add_reaction(model, r.outputs, r.output_coefs, r.inputs, r.input_coefs, r.k_r, rxn_id,
                                      stochastic=stochastic_model, mass_action=r.mass_action)

        return document, model

    def write_sbml_file(self, file_name = None, **keywords):
        document, _ = self.generate_sbml_model(**keywords)
        sbml_string = sbmlutil.libsbml.writeSBMLToString(document)
        f = open(file_name, 'w')
        f.write(sbml_string)
        f.close()
        return f

    def simulate_with_bioscrape_deterministic(self, timepoints, file, initial_condition_dict):
        import bioscrape

        if isinstance(file, str):
            file_name = file
        else:
            file_name = file.name

        m = bioscrape.types.read_model_from_sbml(file_name)
        m.set_species(initial_condition_dict)

        s = bioscrape.simulator.ModelCSimInterface(m)
        s.py_prep_deterministic_simulation()
        s.py_set_initial_time(0)

        sim = bioscrape.simulator.DeterministicSimulator()
        result = sim.py_simulate(s, timepoints)

        return result, m
