
# Copyright (c) 2020, Build-A-Cell. All rights reserved.
# See LICENSE file in the project root directory for details.

from typing import List
from warnings import warn


class Mechanism(object):
    """Mechanism class for core mechanisms.

    Core mechanisms within a mixture (transcription, translation, etc)

    The Mechanism class is used to implement different core
    mechanisms in TX-TL.  All specific core mechanisms should be
    derived from this class.

    """
    def __init__(self, name: str, mechanism_type=""):
        """Initializes a Mechanism instance.

        :param name: name of the Mechanism
        :param mechanism_type: mechanism_type in string
        """
        self.name = name
        self.mechanism_type = mechanism_type
        if mechanism_type == "" or mechanism_type is None:
            warn(f"Mechanism {name} instantiated without a type. This could "
                 "prevent the mechanism from being inherited properly.")

    def update_species(self, component=None, part_id=None) -> List:
        """the child class should implement this method.

        :return: empty list
        """
        warn(f"Default Update Species Called for Mechanism = {self.name}.")
        return []

    def update_reactions(self, component=None, part_id=None) -> List:
        """ the child class should implement this method.

        :return: empty list
        """
        warn(f"Default Update Reactions Called for Mechanism = {self.name}.")
        return []

    def __repr__(self):
        return self.name


class EmptyMechanism(Mechanism):
    """For use when one needs a Mechanism to do nothing, such as translation in Expression Mixtures."""
    def __init__(self, name, mechanism_type):
        """Initializes an EmptyMechanism instance.

        :param name: name of the Mechanism
        :param mechanism_type: mechanism_type in string
        """
        Mechanism.__init__(self, name=name, mechanism_type=mechanism_type)

    def update_species(self, component=None, part_id=None, **keywords):
        return []

    def update_reactions(self, component=None, part_id=None, **keywords):
        return []
