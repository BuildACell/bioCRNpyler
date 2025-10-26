#  Copyright (c) 2020, Build-A-Cell. All rights reserved.
#  See LICENSE file in the project root directory for details.

from typing import List
from warnings import warn


class ComponentEnumerator:
    """Base class for enumerating new components from existing components.

    A `ComponentEnumerator` creates new components in a process similar to
    mechanisms. Component enumerators are used during CRN compilation to
    expand or transform components, generating derived components that are
    then compiled into species and reactions.

    Parameters
    ----------
    name : str
        Name identifier for the component enumerator.

    Attributes
    ----------
    name : str
        The name of the enumerator.

    See Also
    --------
    LocalComponentEnumerator : Enumerator for single-component processing.
    GlobalComponentEnumerator : Enumerator requiring all mixture components.
    Mechanism : Base class for reaction generation.

    Notes
    -----
    This is a base class that should be subclassed to implement specific
    enumeration behavior. The key method to override is
    `enumerate_components`.

    Component enumerators are used during the mixture compilation process
    to:

    - Generate derived components (e.g., transcripts from DNA)
    - Transform components based on context
    - Create component variants or states

    Examples
    --------
    Create a custom component enumerator:

    >>> class MyEnumerator(bcp.ComponentEnumerator):
    ...     def __init__(self):
    ...         super().__init__(name='MyEnumerator')
    ...
    ...     def enumerate_components(self, component=None):
    ...         # Custom enumeration logic
    ...         new_components = []
    ...         # ... generate new components ...
    ...         return new_components

    """

    def __init__(self, name: str):
        self.name = name

    def enumerate_components(self, component=None) -> List:
        """Enumerate new components from an input component.

        This method creates new components based on the input component.
        The base implementation returns an empty list and issues a warning.
        Subclasses should override this method to provide specific
        enumeration behavior.

        Parameters
        ----------
        component : Component, optional
            The input component to enumerate from. Can be None for
            enumerators that don't require an input component.

        Returns
        -------
        list of Component
            List of newly enumerated components. Base implementation
            returns an empty list.

        Warns
        -----
        UserWarning
            Issues a warning when the default implementation is called,
            indicating that a subclass should override this method.

        See Also
        --------
        Component.enumerate_components : Component-level enumeration method.

        Notes
        -----
        This method is called during CRN compilation as part of the
        component enumeration phase. Subclasses should implement specific
        logic to generate derived components.

        """
        warn(
            f"Default update_components called for ComponentEnumerator = "
            f"{self.name}."
        )
        return []

    def __repr__(self):
        """Return string representation of the enumerator.

        Returns
        -------
        str
            The name of the enumerator.

        """
        return self.name


class LocalComponentEnumerator(ComponentEnumerator):
    """Component enumerator that operates on individual components.

    A `LocalComponentEnumerator` processes components independently,
    creating new components based solely on the properties of a single
    input component. This is the most common type of enumerator, used when
    enumeration does not require knowledge of other components in the
    mixture.

    Parameters
    ----------
    name : str
        Name identifier for the local component enumerator.

    Attributes
    ----------
    name : str
        The name of the enumerator (inherited from ComponentEnumerator).

    See Also
    --------
    ComponentEnumerator : Base class for component enumerators.
    GlobalComponentEnumerator : Enumerator requiring all mixture components.

    Notes
    -----
    Local component enumerators are appropriate when:

    - Enumeration depends only on the component being processed
    - No cross-component interactions are needed
    - Components can be processed independently

    Common examples include:

    - Generating RNA transcripts from DNA components
    - Creating protein products from RNA components
    - Expanding DNA assemblies into constituent parts

    The `enumerate_components` method receives only the single component
    being processed.

    Examples
    --------
    Create a custom local enumerator:

    >>> class TranscriptEnumerator(bcp.LocalComponentEnumerator):
    ...     def __init__(self):
    ...         super().__init__(name='TranscriptEnumerator')
    ...
    ...     def enumerate_components(self, component=None):
    ...         if isinstance(component, bcp.DNA):
    ...             # Create RNA transcript from DNA
    ...             transcript = bcp.RNA(name=f'{component.name}_transcript')
    ...             return [transcript]
    ...         return []

    """

    def __init__(self, name: str):
        ComponentEnumerator.__init__(self, name=name)


class GlobalComponentEnumerator(ComponentEnumerator):
    """Component enumerator that operates on all mixture components.

    A `GlobalComponentEnumerator` has access to all components in the
    mixture, allowing for complex enumeration that depends on interactions
    or relationships between multiple components. This is used when
    enumeration decisions require global context.

    Parameters
    ----------
    name : str
        Name identifier for the global component enumerator.

    Attributes
    ----------
    name : str
        The name of the enumerator (inherited from ComponentEnumerator).

    See Also
    --------
    ComponentEnumerator : Base class for component enumerators.
    LocalComponentEnumerator : Enumerator for single-component processing.

    Notes
    -----
    Global component enumerators are appropriate when:

    - Enumeration depends on multiple components
    - Cross-component interactions must be considered
    - Global context or mixture-wide information is needed

    The `enumerate_components` method typically receives all components
    in the mixture, allowing the enumerator to make decisions based on the
    complete set of components.

    Common examples include:

    - Generating complexes between components
    - Creating interaction networks
    - Enumerating components based on global constraints

    **Performance Note:** Global enumerators may be more computationally
    expensive than local enumerators since they must consider all
    components.

    Examples
    --------
    Create a custom global enumerator:

    >>> class ComplexEnumerator(bcp.GlobalComponentEnumerator):
    ...     def __init__(self):
    ...         super().__init__(name='ComplexEnumerator')
    ...
    ...     def enumerate_components(self, component=None):
    ...         # Access all components (passed via mixture)
    ...         # Generate complexes between compatible components
    ...         new_complexes = []
    ...         # ... complex enumeration logic ...
    ...         return new_complexes

    Use in a mixture:

    >>> enumerator = ComplexEnumerator()
    >>> mixture = bcp.Mixture(
    ...     components=[comp1, comp2, comp3],
    ...     global_component_enumerators=[enumerator]
    ... )

    """

    def __init__(self, name: str):
        ComponentEnumerator.__init__(self, name=name)
