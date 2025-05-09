"google adk agent"

import os
import textwrap
from typing import List, Optional

from google.adk.agents import LlmAgent
from google.adk.models.lite_llm import LiteLlm

from gnnepcsaft_mcp_server.plot_utils import v3000_mol_block
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

from .agents_utils import get_gemini_models, get_ollama_models

os.environ["GOOGLE_GENAI_USE_VERTEXAI"] = "False"


# Default model
DEFAULT_MODEL = "gemini-2.0-flash"

# Available models list - can be easily extended
# Initialize with a default list in case fetching fails or returns empty
AVAILABLE_MODELS = [
    "gemini-2.5-flash-preview-04-17",
    "gemini-2.5-pro-preview-05-06",
    "gemini-2.5-pro-exp-03-25",
    "gemini-2.0-flash",
    "gemini-2.0-flash-lite",
    "gemini-1.5-flash",
    "gemini-1.5-pro",
    "gemini-1.5-flash-8b",
]

gemini_models_data = get_gemini_models()
if (
    gemini_models_data
    and "models" in gemini_models_data
    and gemini_models_data["models"]
):
    AVAILABLE_MODELS = gemini_models_data["models"]

ollama_models_data = get_ollama_models()
if ollama_models_data and "models" in ollama_models_data:
    ollama_model_names = [
        f"ollama_chat/{model['name']}" for model in ollama_models_data["models"]
    ]
    AVAILABLE_MODELS.extend(ollama_model_names)

# Ensure DEFAULT_MODEL is in the list, if not, set to the first available or a fallback
if AVAILABLE_MODELS:
    if DEFAULT_MODEL not in AVAILABLE_MODELS:
        DEFAULT_MODEL = AVAILABLE_MODELS[0]
elif (
    not DEFAULT_MODEL
):  # If AVAILABLE_MODELS is empty and DEFAULT_MODEL was also dynamic
    DEFAULT_MODEL = "gemini-2.0-flash"  # Fallback if all fetching fails


all_tools = [
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
    v3000_mol_block,
]


def create_chemistry_agent(model_name=DEFAULT_MODEL):
    """Create a chemistry agent with the specified model"""
    return LlmAgent(
        model=(
            model_name if model_name.startswith("gemini") else LiteLlm(model=model_name)
        ),
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


async def create_root_agent(
    model_name: str = DEFAULT_MODEL, tools: Optional[List] = None
):
    """Create a root agent with the specified model"""
    chemistry_agent_ = create_chemistry_agent(model_name)

    if tools is None:
        tools_ = all_tools
    else:
        tools_ = tools

    return LlmAgent(
        model=(
            model_name if model_name.startswith("gemini") else LiteLlm(model=model_name)
        ),
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
        tools=tools_,
        sub_agents=[chemistry_agent_],
    )
