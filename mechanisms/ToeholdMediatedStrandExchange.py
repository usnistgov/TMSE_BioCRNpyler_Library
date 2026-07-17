#  See LICENSE file in the project root directory for details.

from biocrnpyler.core.mechanism import Mechanism
from biocrnpyler.core.reaction import Reaction
from biocrnpyler.core.species import Species, Complex

# Classes
    
class MaturationO1(Mechanism):
    """
    Mechanism modeling ribozyme-mediated RNA cleavage.

    A 'maturation' mechanism that models cleavage of RNA molecules by ribozymes.
    Ribozymes are catalytic RNA molecules capable of cleaving other RNA molecules at specific sites.
    This mechanism can be used to simulate the effects of targeted RNA degradation or processing
    facilitated by ribozyme activity.

    The reaction follows the schema:
        $$ 'uncleaved_product' --> 'cleaved_product' $$

    where uncleaved_product represents the original RNA before cleavage occurs, and cleaved_product 
    is the resulting product after the cleavage process.

    Parameters
    ----------
    name : str, default='maturation'
        Name identifier for this mechanism instance.
    mechanism_type : str, default='maturation'
        Type classification of this mechanism.

    Attributes
    ----------
    name : str
        Name of the mechanism instance.
    mechanism_type : str
        Type classification ('maturation').

    See Also
    --------
    Mechanism : Base class for all mechanisms.

    Notes
    -----
    The ribozyme cleavage mechanism models targeted degradation or processing of RNA molecules by ribozymes. This can include:

    - Cleavage of specific RNA sequences by ribozymes with high specificity.
    - Production of shorter RNA fragments as a result of the cleavage process.
    - Regulation of gene expression through targeted RNA degradation.

    Required parameters for this mechanism:
    - 'kmat' : Reaction rate constant for ribozyme-mediated cleavage (forward rate).

    Examples
    --------
    Model ribozyme-mediated cleavage of mRNA:
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

    def __init__(self, name='maturation', mechanism_type="maturation", **kwargs):
        Mechanism.__init__(self, name, mechanism_type)

    def update_species(self, rna, **kwargs):
        """Generate species for the ribozyme cleavage mechanism.

        Returns the input RNA molecule involved in the reaction.

        Parameters
        ----------
        rna : Species or list of Species
            The input RNA molecule(s) undergoing cleavage.
        **kwargs
            Additional keyword arguments (unused).

        Returns
        -------
        list of Species
            List containing [rna].
        """
        return rna

    def update_reactions(
        self,
        rna,
        compartment,
        component=None,
        part_id=None,
        **kwargs,
    ):
        """Generate reactions for the ribozyme cleavage mechanism.

        Creates a reaction representing targeted RNA processing by ribozymes with specified rate constants.

        Parameters
        ----------
        rna : Species or list of Species
            The input RNA molecule(s) undergoing cleavage.
        component : Component, optional
            Component containing parameter values. Required if parameters are not provided directly.
        part_id : str, optional
            Identifier for parameter lookup. If None and component is provided, defaults to component.name.
        **kwargs
            Additional keyword arguments (unused).

        Returns
        -------
        list of Reaction
            List containing reactions representing the ribozyme cleavage process.

        Raises
        ------
        ValueError
            If component is None and kmat parameter is not provided.

        Notes
        -----
        The reaction represents targeted RNA degradation with a mass-action model:
          $ 'uncleaved_product' --> 'cleaved_product' $
        """
        # Get Parameters

        def parameter(rna):
            part_id = rna.name[1:]
            kmat = component.get_parameter(
                'kmat', part_id=part_id, mechanism=self
                )
            return kmat

        # Initialize reactions array
        rxns = []
        # Ribozyme cleavage
        for RNA in rna:
            if "pre" in RNA.attributes:
                kmat = parameter(RNA)
                idx_post = rna.index(Species(RNA.name, material_type='rna',
                                             compartment = compartment))
                post_rna = rna[idx_post]
                ribo_rxn = Reaction.from_massaction(
                    inputs=[RNA],
                    outputs=[post_rna],
                    k_forward=kmat  # 1/s
                    )
                rxns.append(ribo_rxn)

        return rxns
    
class GateTxLeak(Mechanism): 
    """Leakage mechanism for unintended transcription of Output strands from Gate strands.

    A 'leak' mechanism that models unintentional transcription of Gate DNA 
    to Output RNA. This mechanism can be used to simulate the effects of 
    unwanted leak reactions that lead to unintended interaction of other 
    stands in the mixture.

    The reaction follows the schema:

        $$ 100*'Gate_RNA_dummy' + n*'Gate_RNA' --> n*'Output_RNA' $$
    
    where 'Gate_RNA' represents the RNA that produces the 'Output_RNA' 
    (leading to the leak), 'Gate_RNA_dummy' represents the total RNA species 
    in the mixture (to maintain mass-balance), and $n$ represents the stoichiometric 
    coefficient corresponding to the leak percentage.
     
    Parameters
    ----------
    name : str, default='tx_leak'
        Name identifier for this mechanism instance.
    mechanism_type : str, default='leak'
        Type classification of this mechanism.

    Attributes
    ----------
    name : str
        Name of the mechanism instance.
    mechanism_type : str
        Type classification ('leak').

    See Also
    --------
    Mechanism : Base class for all mechanisms.

    Notes
    -----
    The leakage mechanism models unintended production of an output strand 
    species. 

    Required parameters for this mechanism:
    - 'ktx' : Reaction rate constant (forward rate).
    - 'leak_percentage' : Percentage of the reaction that leads to leakage.

    Examples
    --------
    Model transcriptional leak to output:

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
        self, name='tx_leak', mechanism_type="leak", **kwargs
    ):
        Mechanism.__init__(self, name, mechanism_type)

    def update_species( self, leak_mol, **kwargs,):
        """Generate species for the leakage mechanism.

        Returns the input molecule and output species involved in the 
        reaction.

        Parameters
        ----------
        leak_mol : Species or list of Species
            The input molecule(s) that can be leaky.
        outputs : list of Species
            List containing possible output molecules that could result 
            from the leak process.
        **kwargs
            Additional keyword arguments (unused).

        Returns
        -------
        list of Species
            List containing [leak_mol, outputs].

        """
        return leak_mol
    
    def update_reactions(
        self,
        leak_dict,
        compartment,
        component=None,
        part_id=None,
        **kwargs,
    ):
        """Generate reactions for the leak mechanism.

        Creates a reaction representing production of output molecules with 
        forward rate constants.

        Parameters
        ----------
        leak_mol : Species or list of Species
            The input molecule(s) undergoing leakage.
        outputs : list of Species
            List containing possible output molecules that could result 
            from the leak process.
        component : Component, optional
            Component containing parameter values.
        part_id : str, optional
            Identifier for parameter lookup. If None and component is 
            provided, defaults to component.name.
        **kwargs
            Additional keyword arguments (unused).

        Returns
        -------
        list of Reaction
            List containing reactions representing the leak process.

        Raises
        ------
        ValueError
            If component is None and ktx or leak_percentage parameters are 
            not provided.

        Notes
        -----
        The reaction represents unintended production with a mass-action 
        model:
        $$
            100*'Gate_RNA_dummy' + n*'Gate_RNA' --> n*'Output_RNA' 
            ('rates': 'ktx', 'leak_percentage') 
        $$

        """
        # Get Parameters
        def get_paramters(leak_mol, k_name, mechanism=self):
            part_id=leak_mol.name
            kvalue = component.get_parameter(
                k_name, part_id=part_id, mechanism=mechanism
            )
            return kvalue
        
        # Get name
        def parse_name(molecule):
            name_parts = molecule.name.split("_")
            if len(name_parts) > 2:
                part_idxi = name_parts[1]
                part_idxj = name_parts[2]
                suffix = '_'+part_idxi+'_'+part_idxj
            else: 
                part_idxi = name_parts[1]
                suffix = '_'+part_idxi
            return suffix
        
        #Initialize reactions array
        rxns = []
        # Leak reactions
        for key in list(leak_dict.keys()):
            substrates = leak_dict[key][0]
            product_leak = leak_dict[key][1]
            gate_species = substrates[1]

            ktx = get_paramters(gate_species, 'ktx', mechanism='transcription') 
            leak_percentage =  get_paramters(gate_species, 'leak_percentage') 
            kleak= int(leak_percentage.value*100) 

            leak_rx = Reaction.from_massaction(
                inputs = [kleak*[substrates[0]],100*[substrates[1]]],
                outputs = kleak*[product_leak],
                k_forward = ktx.value*2) # 1/s
            
            rxns.append(leak_rx)
        return rxns    

class   OneStepTMSE(Mechanism):
    """Mechanism modeling one-step TMSE (Toehold Mediated Strand Exchange) 
    reactions.

    A 'one_step_tmse' mechanism that models TMSE reactions in a single step. 
    This mechanism is used to simulate the interaction between RNA molecules 
    with toeholds where one strand displaces another through base-pairing 
    interactions. This can be applied in various contexts like gene regulation, 
    RNA sensing, or dynamic nucleic acid networks.

    The reaction follows the schema:

        $$ 'strand1' + 'strand2' <--> 'product1' + 'product2'  $$

    where 'strand1' and 'strand2' represent the reacting nucleic acid molecules 
    (such as an input strand and a gate strand), and 'product1' and 'product2' 
    represent the resulting species (such as Output and Gateout) 
    following the toehold-mediated exchange.

    Parameters
    ----------
    name : str, default='one_step_tmse'
        Name identifier for this mechanism instance.
    mechanism_type : str, default="tmse"
        Type classification of this mechanism.

    Attributes
    ----------
    name : str
        Name of the mechanism instance.
    mechanism_type : str
        Type classification ('tmse').

    See Also
    --------
    Mechanism : Base class for all mechanisms.
    TwoStepTMSE: Two-step TMSE mechanism.
    ThreeStepTMSE: Three-step TMSE mechanism.

    Notes
    -----
    The one-step TMSE mechanism models the strand displacement reactions 
    that occur through toeholds. This can include:

    - Interaction between RNA molecules with complementary toehold regions.
    - Production of product strands as a result of strand exchange.

    Required parameters for this mechanism:
    - 'kx' : Reaction rate constant for the forward reaction (forward rate).
    - 'krev' : Reaction rate constant for the reverse reaction (reverse rate, can be 0 if irreversible).
    - 'kth' : Threshold parameter used when waste products are involved.
    - 'krep' : Reporter parameter used when reporters are involved.

    Examples
    --------
    Model TMSE reaction:

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
            self, name='one_step_tmse', mechanism_type="tmse", **kwargs
    ):
        Mechanism.__init__(self, name, mechanism_type)

    def update_species( self,  inputs, gates, outputs, signals, **kwargs,):
        """Generate species for the one-step TMSE mechanism.

        Returns the all RNA molecules involved in the reaction.

        Parameters
        ----------
        inputs : Species or list of Species
            The input RNA molecule(s) undergoing strand displacement.
        gates : Species or list of Species
            Gate RNA molecule(s) with toeholds facilitating reactions.
        outputs : Species or list of Species
            The output RNA molecules resulting from the reaction.
        gateouts : Species or list of Species
            Waste products resulting from the reaction (optional).
        **kwargs
            Additional keyword arguments (unused).

        Returns
        -------
        list of Species
            List containing [inputs, gates, outputs, gateouts].

        """
        return [inputs, gates, outputs, signals]
    
    def update_reactions(
        self,
        inputs,
        gates,
        outputs,
        signals,
        gateouts,
        reporters, 
        rxn_dict,
        compartment,
        component=None,
        part_id=None,
        **kwargs,
    ):
        """Generate reactions for the one-step TMSE mechanism.

        Creates reactions representing strand displacement with specified 
        rate constants. Handles both irreversible and reversible reactions 
        depending on parameter values.

        Parameters
        ----------
        inputs : Species or list of Species
            The input RNA molecule(s) undergoing strand displacement.
        gates : Species or list of Species
            Gate RNA molecule(s) with toeholds facilitating reactions.
        outputs : Species or list of Species
            The output RNA molecules resulting from the reaction.
        signals : Species or list of Species
            Signal products resulting from the reaction (optional).
        gateouts : Species or list of Species
            Waste products resulting from the reaction (optional).
        reporters : Species or list of Species
            The reporter RNA molecules with toeholds facilitating reactions.
        rxn_dict : dict
            Dictionary containing tuples of reactants and products for each reaction.
        component : Component, optional
            Component containing parameter values. Required if parameters are not provided directly.
        part_id : str, optional
            Identifier for parameter lookup. If None and component is provided, defaults to component.name.
        **kwargs
            Additional keyword arguments (unused).

        Returns
        -------
        list of Reaction
            List containing reactions representing the strand displacement process.

        Raises
        ------
        ValueError
            If component is None and kx parameter is not provided.

        Notes
        -----
        The reaction represents strand displacement with a mass-action model. The equation below provides 
        a generic reaction schema; the specific molecular species and applicable rate constants depend on 
        the actual strand types involved in the interaction:
        $$ 
            'strand1' + 'strand2' <--> 'product1' + 'product2' 
            ('rates': 'kx', 'krev', 'kth', 'krep')
        $$

        If krev is zero, the TMSE is irreversible. If a threshold gate is 
        present waste products are produced instead of output and gateout. 
        Furthermore if fuel is present a previous strand is produced and 
        fuelout. Reactions involving fuel and threshold gate are irreversible.

        """
        #Get Paramters
        def get_paramters(subs):
            for sub in subs:
                if 'Gate' in sub.name:
                    part_id = sub.name
                elif 'THErr' in sub.name:
                    part_id =sub.name
                elif 'Fuel' in sub.name:
                    part_id =sub.name
            kx= component.get_parameter(
                'kx', part_id=part_id, mechanism=self
                ) 
            krev= component.get_parameter(
                'krev', part_id=part_id, mechanism=self
                )
            return kx,krev
        
        def get_paramters_thresh(subs):
            for sub in subs:
                if 'TGate' in sub.name:
                    part_id = sub.name
            kth= component.get_parameter(
                'kth', part_id=part_id, mechanism=self
                ) 
            return kth
        
        def get_paramters_krep(subs):
            for sub in subs:
                if 'Reporter' in sub.name:
                    part_id = sub.name
            krep= component.get_parameter(
                'krep', part_id=part_id, mechanism=self
                ) 
            return krep

        # Initialize reactions array
        rxns = []

        for key in list(rxn_dict.keys()):
            if any('Waste' in x.name for x in rxn_dict[key][1]):
                # Get Parameters
                kth = get_paramters_thresh(rxn_dict[key][0])
                rxn= Reaction.from_massaction(
                                inputs = rxn_dict[key][0],
                                outputs = rxn_dict[key][1],
                                k_forward= kth)
                
            elif any('Reporter' in x.name for x in rxn_dict[key][0]):
                # Get Parameters
                krep = get_paramters_krep(rxn_dict[key][0])
                rxn= Reaction.from_massaction(
                                inputs = rxn_dict[key][0],
                                outputs = rxn_dict[key][1],
                                k_forward= krep)
            else:
                # Get Parameters
                kx, kx_rev = get_paramters(rxn_dict[key][0])
                if kx_rev.value == 0:
                    rxn= Reaction.from_massaction(
                                    inputs = rxn_dict[key][0],
                                    outputs = rxn_dict[key][1],
                                    k_forward= kx)
                else:
                    rxn= Reaction.from_massaction(
                                    inputs = rxn_dict[key][0],
                                    outputs = rxn_dict[key][1],
                                    k_forward= kx, k_reverse=kx_rev)
            rxns.append(rxn)
        return rxns

class TwoStepTMSE(Mechanism):
    """Mechanism modeling two-step TMSE (Toehold Mediated Strand Exchange) 
    reactions.

    A 'two_step_tmse' mechanism that models TMSE reactions in two steps. 
    This mechanism is used to simulate the interaction between RNA molecules 
    with toeholds where one strand displaces another through an intermediate 
    complex state before releasing the products.

    The reaction follows the schema:

        $$ 'strand1' + 'strand2' <--> 'complex' <--> 'product1' + 'product2' $$

    where 'strand1' and 'strand2' represent the reacting nucleic acid molecules 
    (such as an input strand and a gate strand), 'complex' represents the bound 
    intermediate state, and 'product1' and 'product2' represent the resulting 
    species (such as Output and Gateout) following the toehold-mediated exchange.

    Parameters
    ----------
    name : str, default='two_step_tmse'
        Name identifier for this mechanism instance.
    mechanism_type : str, default="tmse"
        Type classification of this mechanism.

    Attributes
    ----------
    name : str
        Name of the mechanism instance.
    mechanism_type : str
        Type classification ('tmse').

    See Also
    --------
    Mechanism : Base class for all mechanisms.
    OneStepTMSE: One-step TMSE mechanism.
    ThreeStepTMSE: Three-step TMSE mechanism.

    Notes
    -----
    The two-step TMSE mechanism models the strand displacement reactions that 
    occur through toeholds, explicitly including the intermediate bound complex. 
    This can include:

    - Interaction between RNA molecules with complementary toehold regions to form 
        a complex.
    - Production of output products and gateout byproducts as a result of strand 
        exchange dissociating from the complex.

    Required parameters for this mechanism:
    - 'kb' : Binding forward rate constant.
    - 'kb_rev' : Binding reverse rate constant (can be 0 if irreversible).
    - 'ku' : Unbinding forward rate constant.
    - 'ku_rev' : Unbinding reverse rate constant (can be 0 if irreversible).
    - 'kth' : Threshold parameter used when waste products are involved.
    - 'krep' : Reporter parameter used when reporters are involved.

    Examples
    --------
    Model two-step TMSE reaction:

    >>> input_dict = {'I{1,}': ['Input', 1], 'G{1,2}':['Gate',1,2],}
    >>> strand_mol = StrandExchangeMol(
    ...     mol_dict = input_dict)

    Create a TMSE component:

    >>> tmse_component= TMSE(strand_mol.mol)

    Use with a mixture and TMSE mechanism:

    >>> mixture = Mixture(
    ...     components=[tmse_component],
    ...     mechanisms = {'tmse': OneStepTMSE()},
    ...     parameters = {'kx': 1e-06, 'krev':2.7e-07, 'kb':1e-6, 
     ...                  'kb_rev':2.7e-07,'ku':1e-6, 'ku_rev':2.7e-07}
    ... )
    >>> crn = mixture.compile_crn()

    """
    
    def __init__(
            self, name='two_step_tmse', mechanism_type="tmse", **kwargs
    ):
        Mechanism.__init__(self, name, mechanism_type)

    def update_species( self,  inputs, gates, outputs, gateouts, **kwargs,):
        """Generate species for the two-step TMSE mechanism.

        Returns the all RNA molecules involved in the reaction.

        Parameters
        ----------
        inputs : Species or list of Species
            The input RNA molecule(s) undergoing strand displacement.
        gates : Species or list of Species
            Gate RNA molecule(s) with toeholds facilitating reactions.
        outputs : Species or list of Species
            The output RNA molecules resulting from the reaction.
        gateouts : Species or list of Species
            Waste products resulting from the reaction (optional).
        **kwargs
            Additional keyword arguments (unused).

        Returns
        -------
        list of Species
            List containing [inputs, gates, outputs, gateouts].

        """
        return [inputs, gates, outputs, gateouts]
    
    def update_reactions(
        self,
        inputs,
        gates,
        outputs,
        signals,
        gateouts,
        reporters, 
        rxn_dict,
        compartment,
        component=None,
        part_id=None,
        complex=None,
        **kwargs,
    ):
        """Generate reactions for the two-step TMSE mechanism.

        Creates reactions representing strand displacement with specified 
        rate constants. Handles both irreversible and reversible reactions 
        depending on parameter values.

        Parameters
        ----------
        inputs : Species or list of Species
            The input RNA molecule(s) undergoing strand displacement.
        gates : Species or list of Species
            Gate RNA molecule(s) with toeholds facilitating reactions.
        outputs : Species or list of Species
            The output RNA molecules resulting from the reaction.
        signals : Species or list of Species
            Signal products resulting from the reaction (optional).
        gateouts : Species or list of Species
            Waste products resulting from the reaction (optional).
        reporters : Species or list of Species
            The reporter RNA molecules with toeholds facilitating reactions.
        rxn_dict : dict
            Dictionary containing tuples of reactants and products for each reaction.
        component : Component, optional
            Component containing parameter values. Required if parameters are not 
            provided directly.
        part_id : str, optional
            Identifier for parameter lookup. If None and component is provided, defaults to component.name.
        **kwargs
            Additional keyword arguments (unused).

        Returns
        -------
        list of Reaction
            List containing reactions representing the strand displacement process.

        Raises
        ------
        ValueError
            If component is None and ksd parameter is not provided.

        Notes
        -----
        The reaction represents strand displacement with a mass-action model. The equation below provides 
        a generic reaction schema; the specific molecular species and applicable rate constants depend on 
        the actual strand types involved in the interaction:
        $$ 
            'input' + 'gate' <--> 'complex' <--> 'output' + 'gateout'  
            (rates: 'kb', 'kb_rev', 'ku', 'ku_rev', 'kth', 'krep')
        $$

        If reverse rates ('kb_rev' or 'ku_rev') are zero, the corresponding step is 
        irreversible. If a threshold gate is present, waste products are produced 
        instead of standard output and gateout.If fuel is present, a previous strand 
        is produced along with fuelout. Reactions involving fuel and threshold gates 
        are irreversible.

        """
        #Get Paramters
        def get_paramters(subs):
            for sub in subs:
                if 'Gate' in sub.name:
                    part_id = sub.name
                elif 'Fuel' in sub.name:
                    part_id =sub.name
                elif 'THErr' in sub.name:
                    part_id =sub.name
            # binding
            kb= component.get_parameter(
                'kb', part_id=part_id, mechanism=self
            ) 
            kb_rev= component.get_parameter(
                'kb_rev', part_id=part_id, mechanism=self
            ) 

            # unbinding
            ku= component.get_parameter(
                'ku', part_id=part_id, mechanism=self
            ) 
            ku_rev= component.get_parameter(
                'ku_rev', part_id=part_id, mechanism=self
            ) 
            return kb, kb_rev, ku, ku_rev
        
        def get_paramters_thres(subs):
            for sub in subs:
                if 'TGate' in sub.name:
                    part_id = sub.name
            kth= component.get_parameter(
                'kth', part_id=part_id, mechanism=self
            ) #1e-06,
            return kth
        
        def get_paramters_krep(subs):
            for sub in subs:
                if 'Reporter' in sub.name:
                    part_id = sub.name
            krep= component.get_parameter(
                'krep', part_id=part_id, mechanism=self
            ) 
            return krep

        # Initialize reactions array
        rxns = []

        for key in list(rxn_dict.keys()):
            if any('Waste' in x.name for x in rxn_dict[key][1]):
                # Get Parameters
                kth = get_paramters_thres(rxn_dict[key][0])
                rxn= Reaction.from_massaction(
                                inputs = rxn_dict[key][0],
                                outputs = rxn_dict[key][1],
                                k_forward= kth)
                rxns.append(rxn)
            elif any('Reporter' in x.name for x in rxn_dict[key][0]):
                # Get Parameters
                krep = get_paramters_krep(rxn_dict[key][0])
                rxn= Reaction.from_massaction(
                                inputs = rxn_dict[key][0],
                                outputs = rxn_dict[key][1],
                                k_forward= krep)
                rxns.append(rxn)
            else:
                #Create complex
                if complex is None:
                    complexI = Complex([rxn_dict[key][0]], compartment= compartment)
                else:
                    complexI = complex

                # Get Parameters
                kb_x, kb_rev, ku_x, ku_rev = get_paramters(rxn_dict[key][0])

                # binding reaction
                if kb_rev.value == 0:
                    rxn_b = Reaction.from_massaction(inputs = rxn_dict[key][0],
                                                outputs = [complexI],
                                                k_forward= kb_x)
                else:
                    rxn_b = Reaction.from_massaction(inputs = rxn_dict[key][0],
                                                outputs = [complexI],
                                                k_forward= kb_x, k_reverse=kb_rev)
                rxns.append(rxn_b)

                # unbinding reaction
                if ku_rev.value == 0:
                    rxn_u = Reaction.from_massaction(inputs = [complexI],
                                                    outputs = rxn_dict[key][1],
                                                    k_forward= ku_x)
                else:
                    rxn_u = Reaction.from_massaction(inputs = [complexI],
                                                     outputs = rxn_dict[key][1],
                                                     k_forward= ku_x, k_reverse=ku_rev)
                rxns.append(rxn_u)
        return rxns
    
class ThreeStepTMSE(Mechanism):
    """Mechanism modeling three-step TMSE (Toehold Mediated Strand Exchange) 
    reactions.

    A 'three_step_tmse' mechanism that models TMSE reactions in three steps: 
    binding, branch migration, and unbinding. This mechanism is used to 
    simulate the interaction between RNA molecules with toeholds where one 
    strand displaces another through multiple intermediate complex states.

    The reaction follows the schema:

        $$ 'strand1' + 'strand2' <--> 'complex1' <--> 'complex2' <--> 'output' + 'product1' + 'product2' $$

    where 'strand1' and 'strand2' represent the reacting nucleic acid molecules 
    (such as an input strand and a gate strand), 'complex1' and 'complex2' 
    represent the intermediate binding and branch-migrating states, and 'product1' 
    and 'product2' represent the resulting species (such as Output and Gateout) 
    following the toehold-mediated exchange.

    Parameters
    ----------
    name : str, default='three_step_tmse'
        Name identifier for this mechanism instance.
    mechanism_type : str, default="tmse"
        Type classification of this mechanism.

    Attributes
    ----------
    name : str
        Name of the mechanism instance.
    mechanism_type : str
        Type classification ('tmse').

    See Also
    --------
    Mechanism : Base class for all mechanisms.
    OneStepTMSE: One-step TMSE mechanism.
    TwoStepTMSE: Two-step TMSE mechanism.

    Notes
    -----
    The three-step TMSE mechanism models strand displacement reactions 
    that occur through toeholds, explicitly including intermediate states 
    to more accurately represent binding, branch migration, and dissociation. 
    This can include:

    - Interaction between RNA molecules with complementary toehold regions to form an initial bound complex.
    - Branch migration of the strands transitioning between intermediate complex states.
    - Production of output products and gateout byproducts as a result of 
    strand exchange dissociating from the final complex.

    Required parameters for this mechanism:
    - 'kb' : Binding forward rate constant.
    - 'kb_rev' : Binding reverse rate constant (can be 0 if irreversible).
    - 'kx' : Branch migration forward rate constant.
    - 'kx_rev' : Branch migration reverse rate constant (can be 0 if irreversible).
    - 'ku' : Unbinding forward rate constant.
    - 'ku_rev' : Unbinding reverse rate constant (can be 0 if irreversible).
    - 'kth' : Threshold parameter used when waste products are involved.
    - 'krep' : Reporter parameter used when reporters are involved.

    Examples
    --------
    Model TMSE reaction:

    >>> input_dict = {'I{1,}': ['Input', 1], 'G{1,2}':['Gate',1,2],}
    >>> strand_mol = StrandExchangeMol(
    ...     mol_dict = input_dict)

    Create a TMSE component:

    >>> tmse_component= TMSE(strand_mol.mol)

    Use with a mixture and TMSE mechanism:

    >>> mixture = Mixture(
    ...     components=[tmse_component],
    ...     mechanisms = {'tmse': OneStepTMSE()},
    ...     parameters = {'kx': 1e-06, 'krev':2.7e-07, 'kb':1e-6, 
     ...                  'kb_rev':2.7e-07,'ku':1e-6, 'ku_rev':2.7e-07}
    ... )
    >>> crn = mixture.compile_crn()

    """
    
    def __init__(
            self, name='three_step_tmse', mechanism_type="tmse", **kwargs
    ):
        Mechanism.__init__(self, name, mechanism_type)

    def update_species( self,  inputs, gates, outputs, gateouts, **kwargs,):
        """Generate species for the three-step TMSE mechanism.

        Returns the all RNA molecules involved in the reaction.

        Parameters
        ----------
        inputs : Species or list of Species
            The input RNA molecule(s) undergoing strand displacement.
        gates : Species or list of Species
            Gate RNA molecule(s) with toeholds facilitating reactions.
        outputs : Species or list of Species
            The output RNA molecules resulting from the reaction.
        gateouts : Species or list of Species
            Waste products resulting from the reaction (optional).
        **kwargs
            Additional keyword arguments (unused).

        Returns
        -------
        list of Species
            List containing [inputs, gates, outputs, gateouts].

        """
        return [inputs, gates, outputs, gateouts]
    
    def update_reactions(
        self,
        inputs,
        gates,
        outputs,
        signals,
        gateouts,
        reporters, 
        rxn_dict,
        compartment,
        component=None,
        part_id=None,
        complex=None,
        complex2=None,
        **kwargs,
    ):
        """Generate reactions for the one-step TMSE mechanism.

        Creates reactions representing strand displacement with specified 
        rate constants. Handles both irreversible and reversible reactions 
        depending on parameter values.

        Parameters
        ----------
        inputs : Species or list of Species
            The input RNA molecule(s) undergoing strand displacement.
        gates : Species or list of Species
            Gate RNA molecule(s) with toeholds facilitating reactions.
        outputs : Species or list of Species
            The output RNA molecules resulting from the reaction.
        signals : Species or list of Species
            Signal products resulting from the reaction (optional).
        gateouts : Species or list of Species
            Waste products resulting from the reaction (optional).
        reporters : Species or list of Species
            The reporter RNA molecules with toeholds facilitating reactions.
        rxn_dict : dict
            Dictionary containing tuples of reactants and products for each reaction.
        component : Component, optional
            Component containing parameter values. Required if parameters are not provided directly.
        part_id : str, optional
            Identifier for parameter lookup. If None and component is provided, defaults to component.name.
        **kwargs
            Additional keyword arguments (unused).

        Returns
        -------
        list of Reaction
            List containing reactions representing the strand displacement process.

        Raises
        ------
        ValueError
            If component is None and ksd parameter is not provided.

        Notes
        -----
        The reaction represents strand displacement with a mass-action model. The 
        equation below provides a generic reaction schema; the specific molecular 
        species and applicable rate constants depend on the actual strand types 
        involved in the interaction:

        $$ 
            'strand1' + 'strand2'  <-->  'complex1'  <-->  'complex2'  <-->  'product1' + 'product2'$$
            (rates: 'kb', 'kb_rev', 'kx', 'kx_rev', 'ku', 'ku_rev', 'kth', 'krep')
        $$

        If reverse rates ('kb_rev', 'kx_rev', or 'ku_rev') are zero, the 
        corresponding step is irreversible. If a threshold gate is present, waste 
        products are produced instead of standard output and gateout. If fuel is 
        present, a previous strand is produced along with fuelout. Reactions 
        involving fuel and threshold gates are irreversible.

        """
        #Get Paramters
        def get_paramters(subs):
            for sub in subs:
                if 'Gate' in sub.name:
                    part_id = sub.name
                elif 'Fuel' in sub.name:
                    part_id =sub.name
                elif 'THErr' in sub.name:
                    part_id =sub.name
            # binding
            kb= component.get_parameter(
                'kb', part_id=part_id, mechanism=self
            ) 
            kb_rev= component.get_parameter(
                'kb_rev', part_id=part_id, mechanism=self
            ) 

            # strand-exchange
            kx= component.get_parameter(
                'kx', part_id=part_id, mechanism=self
            ) 
            kx_rev= component.get_parameter(
                'krev', part_id=part_id, mechanism=self
            ) 

            # unbinding
            ku= component.get_parameter(
                'ku', part_id=part_id, mechanism=self
            ) 
            ku_rev= component.get_parameter(
                'ku_rev', part_id=part_id, mechanism=self
            ) 
            return kb, kb_rev, kx, kx_rev, ku, ku_rev
        
        def get_paramters_thres(subs):
            for sub in subs:
                if 'TGate' in sub.name:
                    part_id = sub.name
            kth= component.get_parameter(
                'kth', part_id=part_id, mechanism=self
            ) #1e-06,
            return kth
        
        def get_paramters_krep(subs):
            for sub in subs:
                if 'Reporter' in sub.name:
                    part_id = sub.name
            krep= component.get_parameter(
                'krep', part_id=part_id, mechanism=self
            ) 
            return krep

        # Initialize reactions array
        rxns = []

        for key in list(rxn_dict.keys()):
            if any('Waste' in x.name for x in rxn_dict[key][1]):
                # Get Parameters
                kth = get_paramters_thres(rxn_dict[key][0])
                rxn= Reaction.from_massaction(
                                inputs = rxn_dict[key][0],
                                outputs = rxn_dict[key][1],
                                k_forward= kth)
                rxns.append(rxn)
            elif any('Reporter' in x.name for x in rxn_dict[key][0]):
                # Get Parameters
                krep = get_paramters_krep(rxn_dict[key][0])
                rxn= Reaction.from_massaction(
                                inputs = rxn_dict[key][0],
                                outputs = rxn_dict[key][1],
                                k_forward= krep)
                rxns.append(rxn)
            else:
                #Create complex
                if complex is None:
                    complexI = Complex([rxn_dict[key][0]], compartment= compartment)
                else:
                    complexI = complex
                
                if complex2 is None:
                    complexII = Complex([rxn_dict[key][1]], compartment= compartment)
                else:
                    complexII = complex2

                # Get Parameters
                kb_x, kb_rev, kx, kx_rev, ku_x, ku_rev = get_paramters(rxn_dict[key][0])

                # binding reaction
                if kb_rev.value == 0:
                    rxn_b = Reaction.from_massaction(inputs = rxn_dict[key][0],
                                                   outputs = [complexI],
                                                   k_forward= kb_x)
                else:
                    rxn_b = Reaction.from_massaction(inputs = rxn_dict[key][0],
                                                   outputs = [complexI],
                                                   k_forward= kb_x, k_reverse=kb_rev)
                rxns.append(rxn_b)

                # branch mitgration
                if kx_rev.value == 0:
                    rxn_x = Reaction.from_massaction(inputs = [complexI],
                                                     outputs = [complexII],
                                                     k_forward= kx)
                else:
                    rxn_x = Reaction.from_massaction(inputs = [complexI],
                                                     outputs = [complexII],
                                                     k_forward= kx, k_reverse=kx_rev)
                rxns.append(rxn_x)

                # unbinding reaction
                if ku_rev.value == 0:
                    rxn_u = Reaction.from_massaction(inputs = [complexI],
                                                     outputs = rxn_dict[key][1],
                                                     k_forward= ku_x)
                else:
                    rxn_u = Reaction.from_massaction(inputs = [complexII],
                                                     outputs = rxn_dict[key][1],
                                                     k_forward= ku_x, k_reverse=ku_rev)
                rxns.append(rxn_u)
        return rxns