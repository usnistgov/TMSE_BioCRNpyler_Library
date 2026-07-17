#  See LICENSE file in the project root directory for details.

from typing import List, Union

from biocrnpyler.core.component import Component
from biocrnpyler.core.species import Species

from biocrnpyler.core.mechanism import Mechanism
from biocrnpyler.core.mixture import Mixture
from biocrnpyler.core.parameter import ParameterDatabase
from biocrnpyler.core.reaction import Reaction
from biocrnpyler.components.basic import RNA, Protein
from biocrnpyler.components.dna.rbs import RBS


class RNAassembly(RNA):
    """High-level representation of an RNA-based translation construct.

    An RNAassembly represents a complete RNA translation unit combining a 
    ribosome binding site (RBS), a coding sequence, and the resulting protein 
    product. This class provides a convenient interface for modeling translation 
    directly from an RNA template: RNA --> Protein, where the RBS controls the 
    translation process.

    Parameters
    ----------
    name : str
        Name of the RNA assembly.
    transcript : RNA, str, or None, optional
        The RNA species or name for the assembly. If None, a RNA species
        with `name` is created automatically.
    rbs : RBS, str, or None, optional
        The ribosome binding site component or name controlling translation.
        If None, no translation occurs. If a string, a default RBS is
        created.
    protein : Protein, str, or None, optional
        The protein product of translation. If None, a Protein species with
        `name` is created automatically.
    length : int, optional
        Length of the RNA sequence in base pairs.
    attributes : list of str, optional
        List of attribute tags for the assembly and its species.
    mechanisms : dict or list, optional
        Custom mechanisms for this assembly, overriding mixture defaults.
    compartment : Compartment, optional
        The compartment containing this assembly and its products.
    parameters : dict, optional
        Parameter values specific to this assembly.
    initial_concentration : float, optional
        Initial concentration of the RNA species.
    **kwargs
        Additional keyword arguments passed to the parent `RNA` class.

    Attributes
    ----------
    rna : Species
        The RNA species representing the genetic construct.
    rbs : RBS or None
        The ribosome binding site controlling translation.
    transcript : Species or None
        The RNA transcript produced by transcription.
    protein : Species or None
        The protein product of translation.

    See Also
    --------
    DNAassembly : High-level representation of a gene expression construct.
    RBS : Component representing ribosome binding sites.
    RNA : Base class for RNA components.
    Protein : Base class for protein components.

    Notes
    -----
    The RNAassembly models translation directly from an RNA species. Unlike 
    `DNAassembly`, it does not include a promoter or a transcription step. 
    It is ideal for modeling systems where RNA is introduced directly, such as 
    via transfection, viral delivery, or within an in vitro translation assay.

    Examples
    --------
    Create a basic RNA assembly with automatic species creation:

    >>> rna_construct = bcp.RNAassembly(
    ...     name='rna_gfp',
    ...     rbs='rbs1'
    ... )
    >>> rna_construct.transcript
    rna_rna_gfp
    >>> rna_construct.protein
    protein_rna_gfp

    Create an RNA assembly with custom species names:

    >>> rna_construct = bcp.RNAassembly(
    ...     name='rna_reporter',
    ...     rbs='rbs_strong',
    ...     transcript='mRNA_gfp',
    ...     protein='protein_gfp'
    ... )

    """

    def __init__(
        self,
        name: str,
        transcript=None,
        rbs=None,
        protein=None,
        length=None,
        attributes=None,
        mechanisms=None,
        compartment=None,
        parameters=None,
        initial_concentration=None,
        **kwargs,
    ):
        """Initialize a RNAassembly.

        See class docstring for parameter descriptions.

        """
        self.rbs = None
        self.transcript = None

        # This has to be called at the end so mechanisms are set for
        # the RBS, etc.
        RNA.__init__(
            self,
            name,
            length=length,
            mechanisms=mechanisms,
            parameters=parameters,
            initial_concentration=initial_concentration,
            attributes=attributes,
            compartment=compartment,
            **kwargs,
        )

        self.update_transcript(transcript)
        self.update_protein(protein)
        self.update_rbs(rbs, transcript=self.transcript, protein=self.protein)

    def get_species(self):
        """Get the primary RNA species of this assembly.

        Returns
        -------
        Species
            The RNA species representing this genetic construct.

        """
        return self.transcript

    def set_mixture(self, mixture: Mixture) -> None:
        """Set the mixture containing this component and its sub-components.

        Also propagates the mixture reference to the RBS
        components if they exist.

        Parameters
        ----------
        mixture : Mixture
            The mixture object that contains this assembly.

        """
        self.mixture = mixture
        if self.rbs is not None:
            self.rbs.set_mixture(mixture)

    def update_transcript(
        self, transcript: Union[None, RNA, str], attributes=None
    ) -> None:
        """Set or update the RNA transcript for this assembly.

        Creates an RNA species from the provided input and updates the
        transcript references in RBS components if they exist.

        Parameters
        ----------
        transcript : RNA, str, bool, or None
            The RNA component, species name, False, or None. If None, creates
            an RNA species using the assembly's name. If a string, creates a
            new RNA species with that name. If an RNA object, uses it
            directly. 
        attributes : list of str, optional
            Attribute tags to add to the RNA species.

        Notes
        -----
        This method automatically updates the `transcript` attribute of the
         RBS components to maintain consistency.

        """
        if transcript is None:
            self.transcript = self.set_species(
                self.name,
                material_type='rna',
                attributes=attributes,
                compartment=self.compartment,
            )
        else:
            self.transcript = self.set_species(
                transcript,
                material_type='rna',
                attributes=attributes,
                compartment=self.compartment,
            )

        if self.rbs is not None:
            self.rbs.transcript = self.transcript

    def update_protein(
        self, protein: Union[None, Protein, str], attributes=None
    ) -> None:
        """Set or update the protein product for this assembly.

        Creates a Protein species from the provided input and updates the
        protein references in the RBS components if they exist.

        Parameters
        ----------
        protein : Protein, str, or None
            The Protein component, species name, or None. If None, creates a
            Protein species using the assembly's name. If a string, creates a
            new Protein species with that name. If a Protein object, uses it
            directly.
        attributes : list of str, optional
            Attribute tags to add to the Protein species.

        Notes
        -----
        This method automatically updates the `protein` attribute of the
         RBS components to maintain consistency across the
        assembly.

        """
        if protein is None:
            self.protein = self.set_species(
                self.name,
                material_type='protein',
                attributes=attributes,
                compartment=self.compartment,
            )
        else:
            self.protein = self.set_species(
                protein,
                material_type='protein',
                attributes=attributes,
                compartment=self.compartment,
            )

        if self.rbs is not None:
            self.rbs.protein = self.protein

    def update_rbs(
        self,
        rbs: Union[RBS, str],
        transcript: RNA = None,
        protein: Protein = None,
    ) -> None:
        """Set or update the ribosome binding site component.

        Creates an RBS component from the provided input and propagates the
        assembly's parameters, mixture, and mechanisms to the RBS.

        Parameters
        ----------
        rbs : RBS, str, or None
            The RBS component, RBS name, or None. If None, no RBS is created.
            If a string, creates a default RBS with that name using
            `RBS.from_rbs`. If an RBS object, uses it directly.
        transcript : RNA, optional
            The RNA transcript containing the RBS. If provided, updates the
            assembly's transcript before creating the RBS.
        protein : Protein, optional
            The protein product of translation. If provided, updates the
            assembly's protein before creating the RBS.

        Notes
        -----
        This method automatically:

        - Propagates the assembly's parameter database to the RBS
        - Sets the RBS's mixture reference
        - Adds the assembly's mechanisms to the RBS (without overwriting
          existing RBS mechanisms)

        """
        if protein is not None:
            self.update_protein(protein)

        if transcript is not None:
            self.update_transcript(transcript)

        if rbs is not None:
            self.rbs = RBS.from_rbs(
                name=rbs,
                assembly=self,
                protein=self.protein,
                transcript=self.transcript,
            )
        else:
            self.rbs = None

        if self.rbs is not None:
            self.rbs.update_parameters(
                parameter_database=self.parameter_database,
                overwrite_parameters=False,
            )
            self.rbs.set_mixture(self.mixture)
            self.rbs.add_mechanisms(self.mechanisms, optional_mechanism=True)

    def update_species(self) -> List[Species]:
        """Generate all species associated with this assembly.

        Collects species from the RNA and RBS components during
        CRN compilation.

        Returns
        -------
        list of Species
            List containing the RNA species and all species generated by the
             RBS components.

        Notes
        -----
        This method is called during CRN compilation by
        `Mixture.compile_crn` to collect all chemical species generated by
        this assembly.

        """
        # :return: list of Species that a RNAassemlby instance holds
        species = []
        species.append(self.transcript)

        if self.rbs is not None:
            species += self.rbs.update_species()

        return species

    def update_reactions(self) -> List[Reaction]:
        """Generate all reactions associated with this assembly.

        Collects reactions from the RBS components during CRN
        compilation.

        Returns
        -------
        list of Reaction
            List of all reactions generated by the and RBS
            components, including translation, and regulatory
            reactions.

        Notes
        -----
        This method is called during CRN compilation by
        `Mixture.compile_crn` to collect all chemical reactions generated by
        this assembly.

        """
        # :return: list of Reactions that a RNAassemlby instance holds.
        reactions = []

        if self.rbs is not None:
            reactions += self.rbs.update_reactions()

        return reactions

    def update_parameters(
        self,
        parameter_file: str = None,
        parameters: ParameterDatabase = None,
        overwrite_parameters: bool = True,
    ) -> None:
        """Update parameters for the assembly and its sub-components.

        Propagates parameter updates to the RNA assembly itself and to the
         RBS components if they exist.

        Parameters
        ----------
        parameter_file : str, optional
            Path to a CSV or TSV parameter file to load.
        parameters : ParameterDatabase, optional
            ParameterDatabase object to merge with the assembly's parameters.
        overwrite_parameters : bool, default=True
            If True, new parameter values overwrite existing ones. If False,
            existing parameters are preserved.

        Notes
        -----
        This method calls `update_parameters` on:

        1. The parent RNA class (updating the RNA's parameters)
        2. The RBS component (if it exists)

        """
        RNA.update_parameters(
            self=self,
            parameter_file=parameter_file,
            parameters=parameters,
            overwrite_parameters=overwrite_parameters,
        )

        if self.rbs is not None:
            self.rbs.update_parameters(
                parameter_file=parameter_file,
                parameters=parameters,
                overwrite_parameters=overwrite_parameters,
            )

    def add_mechanism(
        self,
        mechanism: Mechanism,
        mech_type: str = None,
        overwrite: bool = False,
        optional_mechanism: bool = False,
    ) -> None:
        """Add a mechanism to the assembly and its sub-components.

        Adds the mechanism to the assembly's mechanism dictionary and
        propagates it to the RBS components without overwriting
        their existing mechanisms.

        Parameters
        ----------
        mechanism : Mechanism
            The mechanism object to add.
        mech_type : str, optional
            The mechanism type key. If None, uses the mechanism's
            `mechanism_type` attribute.
        overwrite : bool, default=False
            If True, overwrites existing mechanisms with the same type in the
            assembly. If False, raises ValueError for duplicate types.
        optional_mechanism : bool, default=False
            If True, suppresses ValueError when a mechanism key conflict
            occurs in the assembly and `overwrite` is False.

        Notes
        -----
        The mechanism is always added to the RBS with
        `optional_mechanism=True`, meaning it will never overwrite existing
        mechanisms in those components even if `overwrite=True`.

        """
        Component.add_mechanism(
            self,
            mechanism,
            mech_type=mech_type,
            overwrite=overwrite,
            optional_mechanism=optional_mechanism,
        )

        if self.rbs is not None:
            self.rbs.add_mechanism(
                mechanism, mech_type=mech_type, optional_mechanism=True
            )

    def __str__(self):
        return type(self).__name__ + ': ' + self.name

    def __repr__(self):
        txt = str(self)
        if self.rbs is not None:
            txt += '\n\t' + repr(self.rbs)
            txt += '\n\tprotein = ' + repr(self.protein)
        return txt