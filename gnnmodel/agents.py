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
