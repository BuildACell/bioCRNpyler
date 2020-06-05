# mechanism.py - mechanism class for implementing TX-TL mechanisms
# RMM, 16 Aug 2018
#
# Mechanisms the means by which all reactions in a TX-TL reaction are
# established.  Mechanisms can be overridden to allow specialized
# processing of core reactions (eg, adding additional detail, using
# simplified models, etc.
#
# Mechanisms are established in the following order (lower levels
# override higher levels):
#
# Default extract mechanisms
#   Default mechanisms
#       Mechanisms passed to Component() [eg DNA Assembly]
#         Mechanisms based to Sub) [eg, DNA elements]
#
# This hierarchy allows reactions to be created without the user
# having to specify any alternative mechanisms (defaults will be
# used), but also allows the user to override all mechanisms used for
# every (e.g, by giving an alternative transcription
# mechanisms when setting up an extract) or individual mechanisms for
# a given (by passing an alternative mechanism just to that
# .
#
# Copyright (c) 2018, Build-A-Cell. All rights reserved.
# See LICENSE file in the project root directory for details.


from warnings import warn
from .chemical_reaction_network import Species, Reaction, ComplexSpecies, Multimer
from .component import Component
import itertools as it

class Mechanism(object):
    """Mechanism class for core mechanisms

    Core mechanisms within a mixture (transcription, translation, etc)

    The Mechanism class is used to implement different core
    mechanisms in TX-TL.  All specific core mechanisms should be
    derived from this class.

    """
    def __init__(self, name, mechanism_type=""):
        self.name = name
        self.mechanism_type = mechanism_type
        if mechanism_type == "" or mechanism_type is None:
            warn(f"Mechanism {name} instantiated without a type. This could "
                 "prevent the mechanism from being inherited properly.")

    def update_species(self):
        """
        the child class should implement this method
        :return: empty list
        """
        warn(f"Default Update Species Called for Mechanism = {self.name}.")
        return []

    def update_reactions(self, component = None, part_id = None):
        """
        the child class should implement this method
        :return: empty list
        """
        warn(f"Default Update Reactions Called for Mechanism = {self.name}.")
        return []

    def __repr__(self):
        return self.name


class EmptyMechanism(Mechanism):
    def __init__(self, name, mechanism_type):
        Mechanism.__init__(self, name=name, mechanism_type=mechanism_type)

    def update_species(self, **keywords):
        return []

    def update_reactions(self, **keywords):
        return []