#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

from typing import List
from warnings import warn


class ComponentEnumerator:
    def __init__(self, name: str):
        """Class to enumerate components.

        A component enumerator's job is to create new components in a
        process similar to mechanisms.

        """
        self.name = name

    def enumerate_components(self, component=None) -> List:
        """Method to enumerate components.

        This will create new components based on the input component
        somehow. The child class should implement this.

        :return: empty list

        """
        warn(
            f"Default update_components called for ComponentEnumerator = "
            f"{self.name}."
        )
        return []

    def __repr__(self):
        return self.name


class LocalComponentEnumerator(ComponentEnumerator):
    """Class to enumerate local components.

    A component enumerator's job is to create new components in a process
    similar to mechanisms.  A local component enumerator only cares about
    the single component that is passed in

    """

    def __init__(self, name: str):
        ComponentEnumerator.__init__(self, name=name)


class GlobalComponentEnumerator(ComponentEnumerator):
    def __init__(self, name: str):
        """Class to enumerate global components.

        A component enumerator's job is to create new components in a
        process similar to mechanisms.  A global component enumerator takes
        in every component that is in the mixture. This is for complex
        enumeration that cares about other components

        """
        ComponentEnumerator.__init__(self, name=name)
