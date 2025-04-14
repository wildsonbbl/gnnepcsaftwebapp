"google adk agent"

import os
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

os.environ["GOOGLE_GENAI_USE_VERTEXAI"] = "False"

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

chemistry_agent = LlmAgent(
    model=MODEL_NAME,
    name="chemistry_agent",
    instruction=textwrap.dedent(
        """
    You are a chemistry expert that uses InChI and SMILES to analyse molecules. 
    Your only task is to find out the molecule structure, name, and general properties. 
    To be able to do that you are gonna need to take into account all the 
    organic groups known in chemistry, the difference between 
    a Lewis acid and base, the concept of a hydrogen bond, a hydrogen bond donor 
    and hydrogen bond acceptor. Once you have all this information 
    gathered, you will be able to analyse the molecule. 
    You can use theses questions to help yourself describe the molecule in detail.
      - What organic groups does it have?
      - Is it a Lewis acid or base or both? 
      - Can it do hydrogen bonds? 
      - Is it a hydrogen bond donor or acceptor?

    Do not perform any other action. Use the InChI or SMILES 
    info received in previous turns to make the analysis.
    """
    ),
    description="A chemistry expert that handles analyses of molecules from InChI and SMILES.",
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
         first delegate to `pubchem_agent`. 
      - Before delegating, make sure you have available InChI or SMILES,
        if the user didn't provide any of them, ask at least for SMILES.
      - If `pubchem_agent` doesn't find much info on PubChem, then delegate to `chemistry_agent`.
      - If the user wants analysis of a InChI or SMILES, delegate to `chemistry_agent`.
      - Always make sure you have InChI or SMILES before delegating.
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
    sub_agents=[pubchem_agent, chemistry_agent],
)
