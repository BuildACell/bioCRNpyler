from unittest import TestCase
from biocrnpyler import *
import sys

def Component_CRN_validation(CRN, component, mixture):
	#Helper function which ensures components are compiled into "reasonable" CRNs with
	#inputs and outputs (somwhere in the network) that match the functionality

	reaction_inputs = [w.species for r in CRN.reactions for w in r.inputs ]
	reaction_outputs = [w.species for r in CRN.reactions for w in r.outputs ]

	#The below cases test for different kinds of Components
	
	if isinstance(component, DNAassembly):
		G = component.dna
		T = component.transcript
		P = component.protein
		prom = component.promoter
		rbs = component.rbs

		#Test Expression (No Tx/Tl) Convention
		if isinstance(mixture, ExpressionExtract) or isinstance(mixture, ExpressionDilutionMixture):

			#transcript should be ignored in all reactions
			assert T not in reaction_outputs
			assert T not in reaction_inputs

			if  P is not None and prom is not None:
				assert G in reaction_inputs #dna must be an input
				assert P in reaction_outputs #protein must be an output
			else:
				assert P not in reaction_outputs #in this case, no expression
		else:
			if G is not None and T is not None and P is not None and rbs is not None and prom is not None: #Transcription and Translation
				assert G in reaction_inputs #dna must be an input
				assert T in reaction_inputs #transcript must be an input
				assert T in reaction_outputs #transcript must be an output
				assert P in reaction_outputs #protein must be an output
			elif G is not None and T is not None and prom is not None and (P is None or rbs is None): #Just Transcription
				assert G in reaction_inputs #dna must be an input
				assert T in reaction_outputs #protein must be an output
				assert P not in reaction_outputs #No protein output
			else:
				assert T not in reaction_outputs #Otherwise no transcription
				assert P not in reaction_outputs #No translation
	
	if isinstance(component, ChemicalComplex):
		assert component.get_species() in reaction_outputs #Complex can be formed
		assert all([s in reaction_inputs for s in component.internal_species]) #All species in the complex are inputs

	if isinstance(component, Enzyme):
		assert component.enzyme in reaction_inputs #enzyme should be an input
		assert component.enzyme in reaction_outputs #enzyme should be an output
		assert all([s in reaction_inputs for s in component.substrates]) #substrates should be inputs
		assert all([s in reaction_outputs for s in component.products]) #products should be outputs
	

class CombinatorialComponentMixtureTest(TestCase):

	def setUp(self) -> None:
		#Stores a tuple Class, argument keywords dictionary
		#which will be instantiated in each test case below
		self.component_classes = [
			(DNA, {"name":"dna"}),
			(RNA, {"name":"rna"}),
			(Protein, {"name":"protein"}),
			(ChemicalComplex, {"species":[Species("S1"), Species("S2")]}),
			(Enzyme, {"enzyme":Species("E"), "substrates": [Species("S")], "products":[Species("P")]}),
			(DNAassembly, {'name':"dna_assembly_v1", "promoter":"P", "rbs":"R", "transcript":"T", "protein":"X"}),
			(DNAassembly, {'name':"dna_assembly_v2", "promoter":"P", "rbs":"R", "transcript":None, "protein":None}),
			(DNAassembly, {'name':"dna_assembly_v3", "promoter":"P", "rbs":None}),
			(DNAassembly, {'name':"dna_assembly_v4", "promoter":None, "rbs":None}),
			(DNAassembly, {'name':"dna_assembly_v5", "promoter":"P", "rbs":"rbs", "transcript":None, "protein":"X"}),
			(DNAassembly, {'name':"dna_assembly_v6", "promoter":"P", "rbs":"rbs", "transcript":"T", "protein":None})
		]

		#Each of these promoters will be passed into 
		self.promoter_classes = [
			(Promoter, {"name":"promoter"}),
			(ActivatablePromoter, {"name":"activatable_promoter", "activator":Species("A")}),
			(RepressiblePromoter, {"name":"repressible_promoter", "repressor":Species("R")}),
			(RegulatedPromoter, {"name":"regulated_promoter", "regulators":[Species("R1"), Species("R2")]}),
			(CombinatorialPromoter, {"name":"combinatorial_promoter", "regulators":[Species("R1"), Species("R2")]})
		]

		self.rbs_classes = [
			(RBS, {"name":"rbs"})
			]

		self.mixture_classes = [
			(ExpressionExtract, {"name":"expression_extract"}),
			(SimpleTxTlExtract, {"name":"simple_tx_tl_extract"}),
			(TxTlExtract, {"name":"tx_tl_extract"}),
			(ExpressionDilutionMixture, {"name":"expression_dilution_mixture"}),
			(SimpleTxTlDilutionMixture, {'name':"simple_tx_tl_dilution_mixture"}),
			(TxTlDilutionMixture, {"name":"tx_tl_dilution_mixture"}),
			(EnergyTxTlExtract, {"name":"energy_tx_tl_extract"})
			]

		#Use default parameters
		self.parameters = {"kb":1.0, "ku":1.0, "ktx":1.0, "ktl":1.0, "kdeg":1.0, "kdil":1.0, "kexpress":1.0,"kcat":1.0, "K":10, 
		"cooperativity":2, "n":2, "k":1.0, "max_occ":10, "kleak":.01, "k1":1.0, "kbr":1.0, "k2":1.0, "kur":.1, "k_iso":.2, 
		"length":100, "vmax":5, "ktx_solo":1.0, "k_iso_r":10, "ktl_solo":1.0}

		self.transcription_mechs = [
			(SimpleTranscription, {}), 
			(Transcription_MM, {"rnap":Species("RNAP")}),
			(multi_tx, {"pol":Species("RNAP")}),
			]

		self.translation_mechs = [
			(SimpleTranslation, {}),
			(Translation_MM, {"ribosome":Species("Ribosome")}),
			(multi_tl, {"ribosome":Species("Ribosome")}),
			]

	#Instantiate every Component
	def test_component_instantiation(self):
		try:
			for (comp, args) in self.component_classes:
				C = comp(**args)
		except Exception as e:
			error_txt = f"Instantiating Component {comp} with args {args}. \n Unexpected Error: {str(e)}."
			raise Exception(error_txt)

	#Instantiate every promoter
	def test_promoter_instantiation(self):
		try:
			for (prom, args) in self.promoter_classes:
				P = prom(**args)
		except Exception as e:
			error_txt = f"Instantiating Promoter {prom} with args {args}. \n Unexpected Error: {str(e)}."
			raise Exception(error_txt)

	#Instantiate every RBS
	def test_rbs_instantiation(self):
		try:
			for (rbs, args) in self.rbs_classes:
				R = rbs(**args)
		except Exception as e:
			error_txt = f"Instantiating RBS {rbs} with args {args}. \n Unexpected Error: {str(e)}."
			raise Exception(error_txt)

	#Instantiate every Mixture
	def test_mixture_instantiation(self):
		try:
			for (mixture, args) in self.mixture_classes:
				M = mixture(**args)
		except Exception as e:
			error_txt = f"Instantiating Mixture {mixture} with args {args}. \n Unexpected Error: {str(e)}."
			raise Exception(error_txt)

	def test_mechanism_instantiation(self):
		try:
			for (mech, args) in self.transcription_mechs:
				M = mech(**args)
			for (mech, args) in self.translation_mechs:
				M = mech(**args)
		except Exception as e:
			error_txt = f"Instatiating Mechanism {mech}. \n Unexpected Error: {str(e)}."
			raise Exception(error_txt)

	#Instantiate Each Component in each Mixture
	def test_components_in_mixtures(self):
		
		try:
			for (comp, args) in self.component_classes:
				C = comp(**args)
				for (mixture, args) in self.mixture_classes:
					args["parameters"] = dict(self.parameters)
					M = mixture(**args, components = [C])

					if M.get_component(component = C) is None:
						raise AttributeError(f"{C} not found in M.components={M.components}")

					CRN = M.compile_crn()
					document, _ = CRN.generate_sbml_model()
					assert validate_sbml(document) == 0

					#Validate the CRN topology
					Component_CRN_validation(CRN = CRN, component = C, mixture = M)

		except Exception as e:
			error_txt = f"Instantiating Component {comp} in Mixture {mixture} with args {args}. \n Unexpected Error: {str(e)}."
			raise Exception(error_txt)

	#Instantiate Each set of Promoter-RBS in a DNAassembly in each Mixture
	def test_dna_assemblies_in_mixtures(self):
		try:
			for (prom, args) in self.promoter_classes:
				P = prom(**args)
				for (rbs, args) in self.rbs_classes:
					R = rbs(**args)
					for (mixture, args) in self.mixture_classes:
						A = DNAassembly(name = "G", promoter = P, rbs = R)
						args["parameters"] = dict(self.parameters)

						#test adding the component in the construtor
						M = mixture(**args, components = [A])
						if M.get_component(component = A) is None:
							raise AttributeError(f"{A} not found in M.components={M.components}")
						CRN = M.compile_crn()
						document, _ = CRN.generate_sbml_model()
						assert validate_sbml(document) == 0

						#Validate the CRN topology
						Component_CRN_validation(CRN = CRN, component = A, mixture = M)

						#Test adding the component after the constructor
						M2 = mixture(**args)
						M2.add_component(A)
						if M2.get_component(component = A) is None:
							raise AttributeError(f"{A} not found in M.components={M2.components}")

						CRN2 = M2.compile_crn()
						document, _ = CRN2.generate_sbml_model()
						assert validate_sbml(document) == 0

						#Validate the CRN topology
						Component_CRN_validation(CRN = CRN2, component = A, mixture = M)
						
		except Exception as e:
			error_txt = f"Instantiating Promoter {prom} & RBS {rbs} in a DNAassembly in Mixture {mixture} with args {args}. \n Unexpected Error: {str(e)}."
			raise Exception(error_txt)

	def test_dna_assemblies_with_combinations_of_mechanisms(self):
		try:
			for (prom, args) in self.promoter_classes:
				P = prom(**args)
				for (rbs, args) in self.rbs_classes:
					R = rbs(**args)
					for (mech_tx, args_tx) in self.transcription_mechs:
						for (mech_tl, args_tl) in self.translation_mechs:
							Mtx = mech_tx(**args_tx)
							Mtl = mech_tl(**args_tl)
							mechs = {Mtx.mechanism_type:Mtx, Mtl.mechanism_type:Mtl, "binding": One_Step_Cooperative_Binding()}
							A = DNAassembly(name = "G", promoter = P, rbs = R)
							M = Mixture(mechanisms = mechs, components = [A], parameters = self.parameters)

							CRN = M.compile_crn()
							document, _ = CRN.generate_sbml_model()
							assert validate_sbml(document) == 0

							#Validate the CRN topology
							Component_CRN_validation(CRN = CRN, component = A, mixture = M)

		except Exception as e:
			error_txt = f"Instantiating Promoter {prom} & RBS {rbs} in a DNAassembly in Mixure with mech_tx = {mech_tx} and mech_tl = {mech_tl}. \n Unexpected Error: {str(e)}."
			raise Exception(error_txt)

