from unittest import TestCase
from biocrnpyler import *
import sys


class CombinatorialComponentMixtureTest(TestCase):

	def setUp(self) -> None:
		#Stores a tuple Class, argument keywords dictionary
		#which will be instantiated in each test case below
		self.component_classes = [
			(DNA, {"name":"dna"}),
			(RNA, {"name":"rna"}),
			(Protein, {"name":"protein"}),
			(ChemicalComplex, {"species":[Species("S1"), Species("S2")]}),
			(Enzyme, {"enzyme":Species("E"), "substrate":Species("S"), "product":Species("P")}),
			(MultiEnzyme, {"enzyme":Species("E1"), "substrates":[Species("S1"), Species("S2")], "products":[Species("P1"), Species("P2")]}),
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
			(TxTlDilutionMixture, {"name":"tx_tl_dilution_mixture"})
			]

		#Use default parameters
		self.parameters = {"kb":1.0, "ku":1.0, "ktx":1.0, "ktl":1.0, "kdeg":1.0, "kdil":1.0, "kexpress":1.0,"kcat":1.0, "K":10, 
		"cooperativity":2, "n":2, "k":1.0, "max_occ":10, "kleak":.01, "k1":1.0, "kbr":1.0, "k2":1.0, "kur":.1, "k_iso":.2, "ktx_solo":1.0, "k_iso_r":10, "ktl_solo":1.0}

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
						A = DNAassembly("assembly", promoter = P, rbs = R)
						args["parameters"] = dict(self.parameters)

						#test adding the component in the construtor
						M = mixture(**args, components = [A])
						if M.get_component(component = A) is None:
							raise AttributeError(f"{A} not found in M.components={M.components}")
						CRN = M.compile_crn()

						#Test adding the component after the constructor
						M2 = mixture(**args)
						M2.add_component(A)
						if M2.get_component(component = A) is None:
							raise AttributeError(f"{A} not found in M.components={M2.components}")
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
							A = DNAassembly("assembly", promoter = P, rbs = R)
							M = Mixture(mechanisms = mechs, components = [A], parameters = self.parameters)

							CRN = M.compile_crn()
		except Exception as e:
			error_txt = f"Instantiating Promoter {prom} & RBS {rbs} in a DNAassembly in Mixure with mech_tx = {mech_tx} and mech_tl = {mech_tl}. \n Unexpected Error: {str(e)}."
			raise Exception(error_txt)

