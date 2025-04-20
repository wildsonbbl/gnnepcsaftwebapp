"google adk agent"

import os
import textwrap

from google.adk.agents import LlmAgent

from gnnepcsaft_mcp_server.utils import (
    batch_convert_pure_density_to_kg_per_m3,
    batch_critical_points,
    batch_inchi_to_smiles,
    batch_molecular_weights,
    batch_pa_to_bar,
    batch_predict_epcsaft_parameters,
    batch_pure_density,
    batch_pure_h_lv,
    batch_pure_vapor_pressure,
    batch_smiles_to_inchi,
    mixture_density,
    mixture_phase,
    mixture_vapor_pressure,
    pubchem_description,
    pure_phase,
)

os.environ["GOOGLE_GENAI_USE_VERTEXAI"] = "False"

MODEL_NAME = "gemini-2.0-flash"

tools = [
    batch_convert_pure_density_to_kg_per_m3,
    batch_critical_points,
    batch_inchi_to_smiles,
    batch_molecular_weights,
    batch_pa_to_bar,
    batch_predict_epcsaft_parameters,
    batch_pure_density,
    batch_pure_h_lv,
    batch_pure_vapor_pressure,
    batch_smiles_to_inchi,
    mixture_density,
    mixture_phase,
    mixture_vapor_pressure,
    pubchem_description,
    pure_phase,
]

chemistry_agent = LlmAgent(
    model=MODEL_NAME,
    name="chemistry_agent",
    instruction=textwrap.dedent(
        """
    You are a chemistry expert that uses InChI and SMILES to analyse molecules. 
    - Your main task: Find out the molecule structure, name, and general properties. 
    - To be able to do that you are gonna need to take into account all the 
      organic groups known in chemistry, the difference between 
      a Lewis acid and base, the concept of a hydrogen bond, a hydrogen bond donor 
      and hydrogen bond acceptor. Once you have all this information 
      gathered, you will be able to analyse the molecule. 
    - You can use theses questions to help yourself describe the molecule in detail.
        - What organic groups does it have?
        - Is it a Lewis acid or base or both? 
        - Can it do hydrogen bonds? 
        - Is it a hydrogen bond donor or acceptor?
    - You ALWAYS have to check if there's info about the molecule in PubChem with 
      the `pubchem_descriton` tool.
    - Do not perform any other action. Use the InChI or SMILES 
      info received in previous turns to make the analysis.
    - When you can't answer or when you are finished, transfer back to the `gnnepcsaft_agent`.
    """
    ),
    description="A chemistry expert that handles analyses of molecules from InChI and SMILES.",
    tools=[
        pubchem_description,
        batch_smiles_to_inchi,
        batch_inchi_to_smiles,
        batch_molecular_weights,
    ],
)

root_agent = LlmAgent(
    model=MODEL_NAME,
    name="gnnepcsaft_agent",
    description="The main agent in the system. It coordinates the team "
    "and delegates tasks to the other agents.",
    instruction=textwrap.dedent(
        """
    You are the main PCSAFT agent, coordinating a team. 
    - Your main task: Provide information using the ePC-SAFT tools. 
    - You might need more than one tool or function call to solve a problem,
      so you always have to find out what is all the tools you need to solve a problem and
      optimize function calls by calling all the tools you already have 
      arguments information at once.

    - Delegation Rules: 
      - If the user needs more info about a molecule or its properties,
         delegate the task to `chemistry_agent`. 
      - Before delegating, make sure you have available InChI or SMILES,
        if the user didn't provide any of them, ask at least for SMILES.
      - Always make sure you have InChI or SMILES before delegating.
      - Handle thermodynamic requests yourself using ePC-SAFT tools. 
      - For other queries, state clearly if you cannot handle them.
      
    """
    ),
    tools=tools,
    sub_agents=[chemistry_agent],
)
