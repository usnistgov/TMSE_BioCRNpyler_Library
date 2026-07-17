#  See LICENSE file in the project root directory for details.

import warnings
from typing import Union

from biocrnpyler.core.component import Component
from biocrnpyler.core.species import Species
from biocrnpyler.components.dna.assembly import DNAassembly
from RNAassembly import RNAassembly

# Definitions used repeatedly throughout Classes
def identify_indices(dict, key):
    if len(dict[key]) == 2:
        index_i = dict[key][1]
        index_j = dict[key][1]
        suffix = '_'+str(index_i)
    elif len(dict[key])==3:
        index_i = dict[key][1]
        index_j = dict[key][2]

        if type(index_i) is list:
            suffix_i = '_'+ str(index_i[0])
            for indice in index_i[1:]:
                suffix_i += 'n'+str(indice)
            suffix= suffix_i + '_'+str(index_j)
        else:
            suffix = '_'+str(index_i)+'_'+str(index_j)
    else:
        warnings.warn(
            "Domain indices are out of bounds.",
            UserWarning,
            ) 
    return index_i, index_j, suffix

def parse_name(molecule):
    name_parts = molecule.name.split("_")
    if len(name_parts) > 2:
        part_idxi = name_parts[1]
        part_idxj = name_parts[2]
        suffix = '_'+part_idxi+'_'+part_idxj
        if 'n' in part_idxj:
            idxi_parts = part_idxj.split("n")
            part_idxj = [idxi_parts[0], idxi_parts[1]]
    else:
        part_idxi = name_parts[1]
        part_idxj = []
        suffix = '_'+part_idxi
    return part_idxi, part_idxj, suffix

def flatten(t):
    return [item for sublist in t for item in sublist]
    
# Classes
class LeakctRSD(Component):
    """Component representing an RNA molecule with transcriptional leak behavior.

    A `LeakctRSD` component represents a Gate RNA molecule ('Gate' in name) 
    allowing for a leak to an Output RNA to occur. It models the unintended 
    transcription (tx) leak of Gate RNA species to Output species sharing the 
    same indices. The component uses a 'leak' mechanism to generate reactions 
    based on this transcription leak.

    The leak reaction follows the mass-action schema:

        $$100 * Gate_RNA_dummy + n * Gate_RNA --> n * Output_RNA$$

    where 'Gate_RNA' represents the active RNA that produces the 'Output_RNA' 
    (leading to transcriptional leak), 'Gate_RNA_dummy' represents the baseline 
    or inactive RNA species in the mixture, and $n$ represents the stoichiometric 
    coefficient corresponding to the leak percentage.

    Parameters
    ----------
    input_mol : list
        A list containing the input molecules. The molecules inside must have 
        an `.rna` attribute.
    outputs : list, default = None
        List of output species corresponding with the leaky RNA molecules. 
        If None, it will generate outputs based on input molecules.
    compartment : Compartment, default = None
        The compartment in which the component and its leak reactions occur.
    attributes : dict or list, default = None
        Additional attributes to be associated with the component.
    **kwargs :
        Additional keyword arguments passed to the `Component` base class
        constructor.
    
    Attributes
    ----------
    leak_mol : dict
        Dictionary of molecules identified as having leakage behavior mapped 
        to their corresponding reactants and products.
    outputs : list of Species
        List of output RNA species generated.
    compartment : Compartment
        The compartment associated with the component.

    See Also
    --------
    StrandExchangeMol: Molecule capable of strand exchange.
    Component : Base class for biomolecular components.

    Notes
    -----
    Transcription leak primarily models the fractional unintended activation 
    of a Gate species that produces an Output with matching indices.

    Example
    --------
    Create a strand exchange molecule with transcription, RNA maturation
    (ribozyme cutting), and ctRSD leak:

    >>> input_dict = {'G{1,2}':['Gate',1,2],}
    >>> strand_mol = StrandExchangeMol(
    ...     mol_dict = input_dict, 
    ...     transcription = True, 
    ...     rna_maturation = True
    ...     leak = 'ctRSDleak',
    ... )
    >>> ribozyme =RNAmaturation(strand_mol)
    >>> leak_rna = LeakctRSD(strand_mol)

    Use with a mixture and mechanisms:

    >>> mixture = Mixture(
    ...     components = [strand_mol.mol, leak_rna, ribozyme],
    ...     mechanisms = {'transcription': SimpleTranscription(),
    ...                   'maturation': MaturationO1(), 'leak': GateTxLeak()},
    ...     parameters = {'ktx': 0.013, 'leak_percentage': 0.03,'kmat': 0.00417}
    ... )
    >>> crn = mixture.compile_crn()
    
    """

    def __init__(
        self, 
        input_mol: Union[list],
        outputs = None, #what would be the best way of implementing this? Need to define each, what if names are very different, how can we tell the computer? Dict might be better
        compartment = None,
        attributes = None,
        **kwargs,
        ):
                
        # Initialize 
        leak_dict = {}
        leak_rna = []

        if isinstance(input_mol,list):
            rna_mols=[]
            for item in input_mol:
                try:
                    rna_mols.extend(item.rna)
                except AttributeError:
                    raise TypeError("Input molecules do not have .rna attribute.")
        else:
            rna_mols=input_mol.rna
        
        rxn_count= 0
        for rna in rna_mols:
            if 'dummy' in rna.attributes: 
                index_i, index_j, suffix = parse_name(rna)

                if len(index_j)==1:
                    idx_rna =rna_mols.index( 
                                    Species('Gate'+suffix, 
                                            material_type= 'rna',
                                            compartment = compartment,
                                            attributes= ['pre']))
                    rna_notleak = rna_mols[idx_rna]

                    # if outputs is None:
                    output_leak=Species('Output'+suffix, 
                                    material_type = 'rna',
                                    compartment = compartment)
                    output_leak.i = index_i
                    output_leak.j = index_j
                    products= output_leak
                    
                else: 
                    idx_rna =rna_mols.index( 
                                Species('THErr'+suffix, 
                                        material_type= 'rna',
                                        compartment = compartment,
                                        attributes= ['pre']))     
                    rna_notleak = rna_mols[idx_rna]
                                                    
                    name_mol = index_j[1]
                    cds_output = RNAassembly(name_mol, rbs = "RBS",
                                                    compartment = compartment)

                    products=[cds_output.get_species()] 

                    if index_j[0] in (None, ''):
                        suffix = f'_{index_i}'
                    else:
                        suffix = f'_{index_i}_{index_j[0]}'

                    output_leak=Species('Output'+suffix, 
                                    material_type = 'rna',
                                    compartment = compartment)
                    output_leak.i = rna.i
                    output_leak.j = rna.j[0]
                    products.append(output_leak)

                #Add reaction to dictionary
                leak_dict.update({f'rxn_{str(rxn_count)}':
                                    [[rna_notleak, rna], # multiplier on the rna leak will be added in the mechanism part
                                        products]}) 
                rxn_count=+1 
                leak_rna.append(rna)
        #Mol
        self.leak_mol = leak_dict
        self.outputs = leak_rna
        self.compartment = compartment

        # Name the component
        if (name := kwargs.pop('name', None)) is None:
            name = 'mol_temp' 

        Component.__init__(
            self=self, name=name, compartment = compartment, attributes=attributes, **kwargs)

    def get_species(self):
        """Get the list species with leaky transcription.

        Returns
        -------
        list of Species
            List of species with leakage behavior.

        """
        return self
        
    def update_species(self):
        """Use 'leak' mechanism to generate output species.

        Uses the 'leak' mechanism to generate species having in the 
        transcriptional leak.

        Returns
        -------
        list of Species
            List of species generated by the leak mechanism.

        """

        mech_diff = self.get_mechanism('leak')
        return mech_diff.update_species(self.outputs)

    def update_reactions(self):
        """Use 'leak' mechanism to generate transcriptional leak reactions.

        Uses the 'leak' mechanism to generate reactions for leaky 
        transcription in ctRSD.

        Returns
        -------
        list of Reaction
            List of leak reactions (forward) for Gate DNA to Output RNA 
            species.

        """
        mech_diff = self.get_mechanism('leak')
        return mech_diff.update_reactions(
                self.leak_mol, 
                self.compartment,
                component=self, part_id=self.name)
    
 
class StrandExchangeMol(Component):
    """A component representing molecules that can undergo strand exchange.

    A `StrandExchangeMol` component models nucleic acid molecules (DNA and/or 
    RNA) involved in strand exchange networks. It allows for the transcription 
    of DNA species into RNA and subsequent RNA maturation activity (e.g., 
    ribozyme cleavage). The component uses a 'transcription' mechanism to 
    generate reactions based on the specified molecular interactions.

    Parameters
    ----------
    mol_dict : dict
        The dictionary describing the molecular inputs and indices. Must 
        be a dictionary with the structure `{'key1': ['species_type', index_i, index_j]}`, 
        where each molecular species has a unique key and at least one 
        index of type int. `species_type` is a string with options like 
        'Input', 'Gate', and 'TGate' for threshold gates.
    material_type : str, optional
        The material type for input species (e.g., 'dna', 'rna'). Defaults 
        to None (if transcription is enabled, this defaults to 'dna').
    transcription : bool, default=False
        Whether to enable DNA-to-RNA transcription.
    rna_maturation : bool, default=False
        Whether the transcribed RNA molecules undergo single-order maturation 
        (i.e., ribozyme cutting).
    leak : bool, default=False
        Whether to enable transcriptional leak behavior, which generates 
        additional dummy leak RNA species.
    compartment : Compartment, optional
        The compartment in which the component and its reactions occur.
    attributes : list of str, optional
        List of attribute tags to associate with the component.
    **kwargs
        Additional keyword arguments passed to the `Component` base class 
        constructor.

    Attributes
    ----------
    mol : list of Species
        List of input species generated from the input dictionary. If 
        transcription is enabled, these are generated as DNA assemblies.
    rna : list of Species or None
        List of RNA species if transcription is enabled; otherwise None.
    transcription : bool
        Flag indicating whether transcription is enabled.
    maturation : bool
        Flag indicating whether RNA molecules undergo single-order maturation.
    leak : bool
        Flag indicating whether transcriptional leakage behavior is enabled.
    compartment : Compartment
        The compartment associated with the component.

    See Also
    --------
    LeakctRSD : Component representing RNA molecules with leak behavior.
    Component : Base class for biomolecular components.

    Notes
    -----
    Strand exchange models can include transcription (DNA to RNA) and 
    potential ribozyme maturation activity.

    - Transcription without ribozyme: 
        $$ DNA --> RNA + DNA $$

    - Transcription with ribozyme: 
        $$ 
            DNA --> preRNA + DNA
            preRNA --> RNA
        $$
    
    If transcription is not enabled, the strand exchange component simply 
    creates species with the specified indices without generating 
    transcription-based reactions.
    
    Examples
    --------
    Create a strand exchange molecule with transcription and RNA maturation
    (ribozyme cutting):

    >>> input_dict = {'I{1,}': ['Input', 1], 'G{1,2}':['Gate',1,2],}
    >>> strand_mol = StrandExchangeMol(
    ...     mol_dict = input_dict, 
    ...     transcription = True, 
    ...     rna_maturation = True
    ... )

    Use with a mixture and transcription mechanism:

    >>> mixture = Mixture(
    ...     components = [strand_mol.mol],
    ...     mechanisms = {'transcription': SimpleTranscription()},
    ...     parameters = {'ktx': 0.013, 'kmat': 0.00417}
    ... )
    >>> crn = mixture.compile_crn()

    """

    def __init__(
        self, 
        mol_dict: Union[dict],
        material_type = None,
        transcription = False,
        rna_maturation = False,
        leak = False,
        compartment = None,
        attributes = None,
        **kwargs,
    ):

        # Checks mol_dict for compatibility with component
        if isinstance(mol_dict, (tuple, list)):
            if len(mol_dict) != 1 or not isinstance(mol_dict[0], dict):
                raise TypeError("Input must be a dictionary or a " \
                    "tuple/list containing exactly one dictionary.")
            self.mol_dict = mol_dict[0]
        elif isinstance(mol_dict, dict):
            self.mol_dict = mol_dict
        else:
            raise TypeError("Input must be a dictionary.")  
              
        if transcription is True:
            if material_type is None:
                material_type = 'dna'
            elif material_type != 'dna':
                warnings.warn(
                    "Conflicting parameters detected. The 'transcription' parameter takes precedence "
                    "over 'material_type'. Overriding the defined material type to match the transcription setting.",
                    UserWarning,
                    ) 
                material_type = 'dna'
        elif transcription is False and material_type is None:
            warnings.warn(
            "Componentent 'material_type' parameter is not defined. "
            "Execution will continue with default settings, but explicitly defining a 'material_type' is recommended.",
            UserWarning,
            ) 

        name_suffix = ''
        
        # Initialize mol list
        mol_species = []
        rna_species = []

        # Iterate over all molecules to create species
        for key in list(self.mol_dict.keys()):
            index_i, index_j, suffix = identify_indices(self.mol_dict, key)
            mol_name = self.mol_dict[key][0]
            
            if type(index_j) is list:
                if index_j[0] is None:
                    suffix = f'_{index_i}_n{index_j[1]}'
                else:
                    suffix = f'_{index_i}_{index_j[0]}n{index_j[1]}'

            if transcription is True:
                if rna_maturation is True: 
                    new_rna = Species(mol_name+suffix, material_type ='rna', attributes=['pre'],
                                      compartment = compartment) 
                    name_suffix = 'pre'
                else:
                    new_rna = Species(mol_name+suffix, material_type ='rna',
                                      compartment = compartment) 
                
                new_rna.i = index_i
                new_rna.j = index_j
                rna_species.append(new_rna)
                new_transcript = [new_rna]

                if leak == 'ctRSDleak':
                    leak_rna = Species('Output'+suffix, material_type ='rna',
                                      compartment = compartment, attributes = ['dummy'])
                    leak_rna.i = index_i
                    leak_rna.j = index_j
                    new_transcript.append(leak_rna)
                    rna_species.append(leak_rna)

                # make DNA
                new_dna = DNAassembly(new_rna.name, promoter = 'strong', transcript = new_transcript,
                                      compartment = compartment)
                self.rna = rna_species

            else:
                new_dna = Species(mol_name+suffix, material_type = material_type,
                                      compartment = compartment) 

            new_dna.i = index_i
            new_dna.j = index_j
            mol_species.append(new_dna)
        
        # Store the generated molecular species
        self.mol = mol_species
        self.transcription = transcription
        self.maturation = rna_maturation
        self.leak = leak
        self.compartment = compartment

        # Name the component
        if (name := kwargs.pop('name', None)) is None:
            name = name_suffix +'tmse_mol'

        Component.__init__(
            self=self, name=name, compartment = compartment, attributes=attributes, **kwargs)
        
    def get_species(self):
        """Get the input species for strand displacement.

        Returns
        -------
        Component
            The strand displacement component instance itself.

        """

        return self
        
    def update_species(self):
        """Use 'transcription' mechanism to generate strand displacement species.

        If transcription is enabled, uses the 'transcription' mechanism 
        to generate species including the DNA, RNA, and preRNA (if RNA maturation 
        is enabled). Otherwise, does not use any mechanism and simply returns 
        the strand displacement species directly.

        Returns
        -------
        list of Species
            List of species generated by the transcription mechanism, if 
            transcription is enabled. Otherwise, returns the list of species 
            from input molecules.

        """

        if self.transcription is True:
            mech_diff = self.get_mechanism('transcription')
            return mech_diff.update_species(self)
        else:
            return self

    def update_reactions(self):
        """Use 'transcription' mechanism to generate reactions.

        Uses the 'transcription' mechanism to generate reactions for the 
        transcription of strand displacement species, retaining indices 
        from DNA inputs. Additional RNA maturation reactions are generated if 
        RNA maturation is enabled.

        Returns
        -------
        list of Reaction
            List of transcription and RNA maturation reactions if transcription 
            is enabled, otherwise an empty list.

        """
        if self.transcription is True: 
            mech_diff = self.get_mechanism('transcription')
            return mech_diff.update_reactions(
                self.rna, self.compartment, #self.mol, 
                component=self, part_id=self.name
            )
        else:
            return []


class TMSE(Component):
    """Molecules that undergo two-domain mediated strand exchange.

    A `TMSE` (Two-domain Mediated Strand Exchange) component represents 
    an abstraction of a process involving molecular species interactions,
    specifically modeling strand exchanges between input strands, gate 
    strands, output strands, fuel strands, threshold gates, and reporters. 
    The component uses a 'tmse' mechanism to generate reactions based on 
    these strand interactions.

    Parameters
    ----------
    mol_species : list or dict
        A list of initial molecular species (or a dictionary defining them) 
        that will be processed and categorized into inputs, gates, outputs, 
        fuel strands, threshold gates, and reporters.
    compartment : Compartment, optional
        The compartment in which the component and its reactions occur.
    attributes : list of str, optional
        List of attribute tags to associate with the component.
    **kwargs
        Additional keyword arguments passed to the `Component` base class
        constructor.

    Attributes
    ----------
    transcripts : list of Species
        List of RNA transcripts associated with error-handling/leak species.
    inputs : list of Species
        List of input molecular species.
    gates : list of Species
        List of gate molecular species.
    outputs : list of Species
        List of output molecular species generated from interactions.
    reporters : list of Species
        List of reporter molecular species.
    gateouts : list of Species
        Intermediate molecular species resulting from gate interactions.
    signals : list of Species
        Signal products resulting from interactions with reporters.
    rxn_dict : dict
        Dictionary mapping reaction IDs to their corresponding reactants 
        and products.
    compartment : Compartment
        The compartment associated with the component.

    See Also
    --------
    StrandExchangeMol : Class to create RNA/DNA species for TMSE circuit.
    Component : Base class for biomolecular components.

    Notes
    -----
    The `TMSE` component models strand exchange networks by evaluating the 
    domain complementarities (using `i` and `j` indices) between all provided 
    strands. It categorizes them and maps out their valid interactions, storing 
    them in `rxn_dict`. It uses an assigned 'tmse' mechanism (e.g., `OneStepTMSE`,
    `TwoStepTMSE`, or `ThreeStepTMSE`) to generate the final CRN reactions.

    Examples
    --------
    Create RNA strands:

    >>> input_dict = {'I{1,}': ['Input', 1], 'G{1,2}':['Gate',1,2],}
    >>> strand_mol = StrandExchangeMol(
    ...     mol_dict = input_dict)

    Create a TMSE component:

    >>> tmse_component= TMSE(strand_mol.mol)

    Use with a mixture and TMSE mechanism:

    >>> mixture = Mixture(
    ...     components=[tmse_component],
    ...     mechanisms = {'tmse': OneStepTMSE()},
    ...     parameters = {'kx': 1e-06, 'krev':2.7e-07}
    ... )
    >>> crn = mixture.compile_crn()

    """

    def __init__(
            self,
            mol_species: Union[list, dict],
            compartment = None, 
            attributes=None,
            **kwargs,
    ):

        if isinstance(mol_species, dict):
            species_list= [] #Initialize reporter list
            for key in list(mol_species.keys()):
                strand_name = mol_species[key][0]
                index_i, index_j, suffix = identify_indices(mol_species, key)         
                spec= Species(strand_name+suffix, material_type = 'dna', compartment = compartment)
                spec.i = index_i
                spec.j = index_j
                species_list.append(spec)
        #else if list of species
        elif isinstance(mol_species, list):
            # Check if it is mix list of array and dicitonary
            species_list = []
            for item in mol_species:
                if isinstance(item, list):
                    species_list.extend(item)
                else:
                    species_list.append(item)

            if len(species_list) > 0 :
                # Check Species for indices
                for spec in species_list:
                    try: 
                        spec.i, spec.j
                    except AttributeError:
                        raise TypeError(f"Domain indices are not defined for species {spec}")
        else:
            raise TypeError("Input must be a dictionary or a list of species.")

        #Separate and sort different types of strands
        input_strands = []
        gate_strands = []
        fuel_strands = []
        tgate_strands =[]
        therr_strands =[]
        reporter_strands = []

        # for strand in species_list and [if 'pre' is not strand.attributes]:
        for strand in [s for s in species_list if 'pre' not in s.attributes 
                                                    and 'dummy' not in s.attributes]:
            # if 'Input' in strand.name:
            #     input_strands.append(strand)
            
            if 'Gate' in strand.name: # and 'u' not in strand.name:
                if 'T' in strand.name:
                    tgate_strands.append(strand)
                else:
                    gate_strands.append(strand)
            elif 'Fuel' in strand.name:
                fuel_strands.append(strand)
            elif 'THErr' in strand.name:
                therr_strands.append(strand)
            elif 'Reporter' in strand.name:
                reporter_strands.append(strand)
            elif 'Input' in strand.name or 'Output' in strand.name: # Explicit check instead of catch-all 'else'
                input_strands.append(strand)
            else:
                raise ValueError(f"Strand '{strand.name}' does not have a designated strand displacement type.")

        #Products and reaction dict {[subs],[prod],[rxn_type]}
        output_strands = []
        gateout_strands = []
        transcript_strands = []
        tsme_dict = {}
        rxn_num = 0 #Counter for reactions

        def check_strand_interactions(list_1, list_2, rxn_count, 
                                      output_strands, gateout_strands):
            attributes = None
            for strand1 in list_1:
                material_type1 = strand1.material_type
                for strand2 in list_2:
                    material_type2 = strand2.material_type

                    if material_type1 == material_type2:
                        Out_material_type = material_type1
                        GateOut_material_type = material_type1
                    else:
                        Out_material_type = material_type2
                        GateOut_material_type = material_type1+material_type2

                        if Out_material_type == 'dna': 
                            attributes= ['ss']

                    if strand1.j == strand2.i:
                        index_i, index_j = [strand2.i, strand2.j]
                        index_k = strand1.i
                        
                        suffix_O = '_'+str(index_i) +'_'+str(index_j)
                        suffix_GO = '_'+str(index_k) +'_'+str(index_i)

                        # Outputs
                        if any('Output'+suffix_O in x.name for x in output_strands):
                            idx_out =output_strands.index(
                                Species('Output'+suffix_O,
                                        material_type= Out_material_type,
                                        compartment = compartment)) 
                            output = output_strands[idx_out]
                        else:
                            output= Species('Output'+suffix_O, 
                                            material_type=Out_material_type,
                                            attributes=attributes,
                                            compartment = compartment)
                            output.i = index_i
                            output.j = index_j
                            output_strands.append(output)

                        # GateOut 
                        gateout = Species('GateOut'+suffix_GO, 
                                          material_type=GateOut_material_type,
                                          compartment = compartment) 
                        gateout.i= index_k
                        gateout.j = index_i

                        gateout_strands.append(gateout)
                        
                        #Add reaction to DF
                        tsme_dict.update({f'rxn_{str(rxn_count)}':
                                          [[strand2, strand1], 
                                           [output, gateout]]}) 
                        rxn_count+=1
                if rxn_num == 0:
                    TypeError('No interactions found.')
            return output_strands, gateout_strands, rxn_count
        
        def check_fuel_interactions(fuel_strands, gateout_strands, 
                                    rxn_count, rna0):
            fuelout_strands = [] 
            for strand1 in fuel_strands:
                material_type1=strand1.material_type
                for strand2 in gateout_strands:
                    material_type2=strand2.material_type

                    if material_type1==material_type2:
                        material_type =material_type1
                    else:
                        material_type= None

                    if strand2.j == strand1.i:
                        index_i, index_j = [strand1.i, strand2.j]
                        index_k = strand2.i
                        
                        suffix = '_'+str(index_k) +'_'+str(index_j)

                        if index_k == index_j: 
                            suffix = '_'+str(index_k) 

                        suffix_FO = '_'+str(index_i)

                        # Fuel species 
                        if any('Input'+suffix in x.name for x in rna0):
                            idx_recyl =rna0.index(
                                Species('Input'+suffix, 
                                        material_type= material_type,
                                        compartment = compartment))

                        elif any('Output'+suffix in x.name for x in rna0): 
                            idx_recyl =rna0.index(
                                Species('Output'+suffix,
                                        material_type= material_type, 
                                        compartment = compartment))
                        else:
                            warnings.warn(
                                "Fuel strands have no interactions.", 
                                UserWarning,
                                ) 

                        # FuelOut 
                        fuelout = Species('FuelOut'+suffix_FO, 
                                          material_type=material_type,
                                          compartment = compartment) 
                        fuelout.i= index_j
                        fuelout.j = index_j
                        fuelout_strands.append(fuelout)

                        rna_recyl = rna0[idx_recyl]
                        #Add reaction to DF
                        tsme_dict.update({f'rxn_{str(rxn_count)}':
                                          [[strand2, strand1],
                                           [rna_recyl, fuelout]]}) 
                        rxn_count+=1
            return rxn_count

        def check_tg_interactions(tgate_strands, output_strands, rxn_count):
            for strand1 in tgate_strands:
                material_type1=strand1.material_type
                for strand2 in output_strands:
                    material_type2=strand2.material_type
                    if material_type1==material_type2:
                        material_type =material_type1
                    else:
                        material_type= material_type1+material_type2 
                    if (strand1.i==strand2.j) or (strand1.i == strand2.i and strand1.j == strand2.j): #Check And statement here
                        index_i, index_j = [strand2.i, strand2.j]
                        suffix = '_'+str(index_i) +'_'+str(index_j)
                        suffix_FO = '_'+str(index_j)

                        # Waste
                        waste_0 = Species('Waste'+suffix_FO, 
                                          material_type=material_type,
                                          compartment = compartment)
                        waste_0.i= index_i

                        waste_1 = Species('Waste'+suffix, 
                                          material_type=material_type,
                                          compartment = compartment) 
                        waste_1.i= index_i
                        waste_1.j= index_j
                        
                        #Add reaction to DF
                        tsme_dict.update({f'rxn_{str(rxn_count)}':
                                          [[strand2, strand1], 
                                           [waste_0, waste_1]]}) 
                        rxn_count+=1
            return rxn_count
                
        #Initialize signal-related species
        signal_strands = []
        reporterout_strands = []

        def check_therr_interactions(input_list, therr_list, rxn_count):
            for therr_strand in therr_list:
                if therr_strand.material_type  == 'rna':
                    for strand in input_list:
                        if strand.j == therr_strand.i:
                            index_i, index_j = [therr_strand.i, therr_strand.j]                            
                            if type(index_j[1]) is str:
                                name_mol = index_j[1]
                            else:
                                name_mol = str(index_j[1])
                            
                            cds_output = RNAassembly(name_mol, rbs = "RBS",
                                                compartment = compartment)
                            transcript_strands.append(cds_output) #transcript_strands
                            product = cds_output.transcript #not giving the rnaassemble

                            if index_j[0] is None:
                                suffix = f'_{index_i}'
                            else:
                                suffix = f'_{index_i}_{index_j[0]}'

                            output = Species('Output'+suffix, 
                                            material_type='rna',
                                            compartment = compartment)
                            
                            output.i= index_i
                            output.j = index_j[0]
                            output_strands.append(output)
                            
                            #Add reaction to DF
                            tsme_dict.update({f'rxn_{str(rxn_count)}':
                                            [[strand, therr_strand], 
                                            [product, output]]}) 
                            rxn_count+=1
                else:
                    raise TypeError('THErr needs to be RNA')
                        
                if rxn_num == 0:
                    TypeError('No interactions found.')
            return output_strands, gateout_strands, rxn_count

        def reporter_strand_interactions(list_reporter, list_2, rxn_count, signal_strands): #check report not interacting
            for strand1 in list_reporter:
                material_type1 = strand1.material_type
                for strand2 in list_2:
                    material_type2=strand2.material_type
                    attributes = None

                    # determine the material type of signal and reporter out
                    if material_type1==material_type2: 
                        signal_material_type = material_type1
                        reporterout_material_type = material_type1
                    else: # Signal take material type of reporter
                        signal_material_type = material_type1
                        reporterout_material_type = material_type2 + material_type1
                        if material_type1 == 'dna':
                            attributes = ['ss']

                    if strand1.j == strand2.i or strand1.i==strand2.j:
                        if strand1.j == strand2.i:
                            index_i, index_j = [strand1.j, strand2.j]
                        else:
                            index_i, index_j = [strand2.j,strand1.j]

                        # Signals
                        signal= Species('Signal_'+str(index_i), 
                                        material_type = signal_material_type,
                                        attributes = attributes,
                                        compartment = compartment) #Using gate.j as naming convention
                        signal.i = index_i
                        signal.j = index_j
                        signal_strands.append(signal)

                        # ReporterOut
                        reporterout = Species('ReporterOut_'+str(index_i), material_type = reporterout_material_type,
                                              compartment = compartment) #Using gate.i as naming convention
                        reporterout.i = index_i
                        reporterout.j = index_j

                        reporterout_strands.append(reporterout)
                        
                        #Add reaction to dictionary
                        tsme_dict.update({f'rxn_{str(rxn_count)}':
                                            [[strand2, strand1], 
                                             [signal, reporterout]]}) 
                        rxn_count+=1
                if rxn_num == 0:
                    TypeError('No interactions found.')
        
            return signal_strands, rxn_count
        
        rxn_num = 0
        #Checking Input-Gates Strand interactions
        output_strands, gateout_strands, rxn_num= check_strand_interactions(
            input_strands, gate_strands, rxn_num, output_strands, gateout_strands)
            
        #Check for Gate-Output Strand interactions
        output_strands, gateout_strands, rxn_num= check_strand_interactions(
            output_strands, gate_strands,rxn_num, output_strands, gateout_strands)
               
        #Check for Fuel
        if len(fuel_strands)>0:
            rxn_num= check_fuel_interactions(fuel_strands, gateout_strands, 
                                             rxn_num, input_strands+output_strands)

        #Check for TG
        if len(tgate_strands)> 0:
            rxn_num= check_tg_interactions(tgate_strands, 
                                           input_strands+output_strands, rxn_num)
            
        #Check for THErr-Input Strand interactions ##Previously after reporter
        output_strands, gateout_strands, rxn_num= check_therr_interactions(
            output_strands+input_strands,therr_strands, rxn_num)
            
        #Check Reporter against Outputs
        if len(reporter_strands)> 0:
            signal_strands, rxn_num = reporter_strand_interactions(
                reporter_strands, output_strands+input_strands, rxn_num, signal_strands)
    
        #Does not check for Output-Output strand or THErr gate byproduct interactions.
        
        #Strand substrates and products
        self.transcripts = transcript_strands
        self.inputs = input_strands
        self.gates = gate_strands
        self.outputs = output_strands
        self.reporters = reporter_strands
        self.gateouts= gateout_strands
        self.signals = signal_strands
        self.rxn_dict = tsme_dict
        self.compartment =compartment

        # Name the component
        if (name := kwargs.pop('name', None)) is None:
            name = 'TMSE'

        Component.__init__(
            self=self, name=name, compartment = compartment, attributes=attributes, **kwargs)

    def get_species(self):
        """Get all species involved in strand displacement.

        Returns
        -------
        Species
            The species involved in strand displacement.

        """
        return self.transcripts
        
    def update_species(self):
        """Use the 'tmse' mechanism to generate species.

        Uses the assigned 'tmse' mechanism to generate the full list of 
        interacting strand species based on the inputs, gates, outputs, 
        and signals identified during initialization.

        Returns
        -------
        list of Species
            List of species returned by the TMSE mechanism.

        """
        mech_diff = self.get_mechanism('tmse')
        return mech_diff.update_species(self.inputs, self.gates, self.outputs, self.signals)

    def update_reactions(self):
        """Use the 'tmse' mechanism to generate strand exchange reactions.

        Uses the assigned 'tmse' mechanism (e.g., OneStepTMSE, TwoStepTMSE) 
        to compile the abstract interaction pairings mapped in `rxn_dict` 
        into concrete reaction steps.

        Returns
        -------
        list of Reaction
            List of reactions generated by the TMSE mechanism governing 
            the inputs, gates, outputs, signals, gateouts, and reporters.

        """
        mech_diff = self.get_mechanism('tmse')
        return mech_diff.update_reactions(
            self.inputs, self.gates,
            self.outputs, self.signals, self.gateouts, self.reporters, self.rxn_dict, self.compartment,
            component=self, part_id=self.name
        )

class RNAmaturation(Component): 
    """Component representing an RNA molecule that undergoes maturation processes.

    An `RNAmaturation` component represents an abstraction of a biological 
    process involving RNA maturation, specifically modeling the cleavage of 
    precursor RNA molecules into functional fragments. The component uses a 
    'maturation' mechanism to simulate the catalytic activity of RNA cleavage 
    (e.g., ribozyme activity).

    Parameters
    ----------
    mols : list or object
        List of objects (or a single object) containing RNA species to be 
        processed. The objects inside must have an `.rna` attribute containing 
        `Species` objects.
    attributes : list of str, optional
        List of attribute tags to associate with the component and its species.
    compartment : Compartment, optional
        The compartment in which the component and its reactions occur.
    **kwargs
        Additional keyword arguments passed to the `Component` base class
        constructor.

    Attributes
    ----------
    rna : list of Species
        Combined list of precursor (uncut) and processed (mature) RNA species 
        resulting from the maturation process.
    compartment : Compartment
        The compartment associated with the component.

    See Also
    --------
    Component : Base class for biomolecular components.

    Notes
    -----
    The component identifies precursor RNA by checking for the `"pre"` tag 
    in the species' attributes. It then creates a mature version without the 
    `"pre"` attribute. 
    
    The cleavage mechanism generates a unidirectional mass-action reaction:
    $$ pre-RNA --> RNA $$

    Examples
    --------
    Create a ribozyme cleavage of RNA:
    
    >>> input_dict = {'I{1,}': ['Input', 1], 'G{1,2}':['Gate',1,2],}
    >>> strand_mol = StrandExchangeMol(
    ...     mol_dict = input_dict, 
    ...     transcription = True, 
    ...     rna_maturation = True
    ... )
    >>> ribozyme = RNAmaturation(
    ...     mols = strand_mol
    ... )

    Use with a mixture and maturation mechanism:

    >>> mixture = Mixture(
    ...     components=[ribozyme],
    ...     mechanisms = {'maturation': MaturationO1()},
    ...     parameters = {'kmat': 0.00417}
    ... )
    >>> crn = mixture.compile_crn()

    """

    def __init__(
        self, 
        mols: Union[list],
        attributes = None,
        compartment = None,
        **kwargs,
    ):     
        # Initialize cut gates
        mod_rna = []

        if isinstance(mols,list):
            rna_mols=[]
            for item in mols:
                try:
                    rna_mols.extend(item.rna)
                except AttributeError:
                    raise TypeError("Input molecules do not have .rna attributes.")
        else:
            rna_mols=mols.rna

        for rna in rna_mols:
            if "pre" in rna.attributes:
                mod_rna.append(rna)
                index_i = rna.i
                index_j = rna.j
                post_rna = Species(rna.name, material_type='rna',
                                   compartment = compartment) 
                post_rna.i = index_i
                post_rna.j = index_j
                mod_rna.append(post_rna)
            if len(mod_rna) == 0:
                raise TypeError("Unable to identify pre-RNA.")

        #RNA
        self.rna = mod_rna
        self.compartment = compartment
        # Name the component
        if (name := kwargs.pop('name', None)) is None:
            name = 'ribozyme'

        Component.__init__(
            self=self, name=name, compartment = compartment, attributes=attributes, **kwargs)

    def get_species(self):
        """Get all RNA species, including those produced by maturation.

        Returns
        -------
        list of Species
            The RNA species including those involved in ribozyme cleavage.

        """
        return self.rna
        
    def update_species(self):
        """Use 'maturation' mechanism to generate processed species.

        Uses the 'maturation' mechanism to generate mature RNA from 
        precursor RNA (preRNA).

        Returns
        -------
        list of Species
            List of all RNA species, including the generated mature RNA 
            produced by the cleavage mechanism.

        """
        mech_diff = self.get_mechanism('maturation')
        return mech_diff.update_species(self.rna)

    def update_reactions(self):
        """Use 'maturation' mechanism to generate maturation reactions.

        Uses the 'maturation' mechanism to generate the specific reactions 
        converting precursor RNA into mature RNA.

        Returns
        -------
        list of Reaction
            List of RNA maturation (e.g., ribozyme cleavage) reactions.

        """
        mech_diff = self.get_mechanism('maturation') 
        return mech_diff.update_reactions(
            self.rna,  self.compartment, component=self, part_id=self.name
            )