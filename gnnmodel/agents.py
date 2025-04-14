"google adk agent"

import textwrap

from gnnepcsaft.data.rdkit_util import inchitosmiles, mw, smilestoinchi
from gnnepcsaft.epcsaft.epcsaft_feos import (
    mix_den_feos,
    mix_vp_feos,
    pure_den_feos,
    pure_h_lv_feos,
    pure_vp_feos,
)
from google.adk.agents import LlmAgent

from .utils import mixture_phase, prediction, pubchem_description, pure_phase

PROMPT = """
    I have two molecules with SMILES 'CCO' and 'CC(=O)C'. 
    Consider them pure at 325 K and 101325 Pa.
    - What's their density? Which one is more dense?
    - What's their vapor pressure? Which one is more volatile?
    - What's their pure component phase.
    Then consider their 50/50 mixture at 300 K and 101325 Pa.
    - What's their density?
    - What's the mixture phase.
    What's the molecule with those SMILES?
    If there's any source info about it, what does each source say?
    """


MODEL_NAME = "gemini-2.0-flash-exp"

pubchem_agent = LlmAgent(
    model=MODEL_NAME,
    name="pubchem_agent",
    instruction=textwrap.dedent(
        """
    You are the PubChem agent. 
    Your only task is to provide information for the 
    gnnepcsaft_agent using the pubchem_descriton tool. 
    Do not perform any other action. Use the InChI or SMILES 
    info received in previous turns to use the tools.
    """
    ),
    description="Handles gathering of information from PubChem using InChI and SMILES.",
    tools=[pubchem_description, smilestoinchi, inchitosmiles],
)

root_agent = LlmAgent(
    model=MODEL_NAME,
    name="gnnepcsaft_agent",
    description=textwrap.dedent(
        """
    You are the main PCSAFT agent, coordinating a team. 
    - Your main task: Provide information using the ePC-SAFT tools. 
      Handle its 'status' response ('report' or 'error_message'). 
    - Delegation Rules: 
      - If you or the user needs more info about a molecule 
         delegate to `pubchem_agent`. 
      - Before delegating, make sure you have available InChI or SMILES,
        if the user didn't provide any of them, ask at least for SMILES.
      - Handle thermodynamic requests yourself using ePC-SAFT tools. 
      - For other queries, state clearly if you cannot handle them.
    """
    ),
    tools=[
        pure_vp_feos,
        pure_den_feos,
        mix_den_feos,
        mix_vp_feos,
        pure_phase,
        mixture_phase,
        mw,
        prediction,
        smilestoinchi,
        inchitosmiles,
        pure_h_lv_feos,
    ],
    sub_agents=[pubchem_agent],
)
