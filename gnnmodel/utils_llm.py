"Module for utils to work with LLMs"

import os
import textwrap
from json import loads
from urllib.parse import quote
from urllib.request import HTTPError, urlopen

from langchain_core.prompts import ChatPromptTemplate
from langchain_google_genai import ChatGoogleGenerativeAI
from markdown import markdown
from pydantic import SecretStr

from . import logger


def resume_mol(inchi: str, smiles: str, api_key: SecretStr):
    "Describe the molecule with google's gemini."

    llm = ChatGoogleGenerativeAI(
        model="gemma-3-27b-it",
        api_key=api_key if api_key else SecretStr(os.environ.get("GEMINI_API_KEY", "")),
    )

    ans = pubchem_description(inchi)

    query = """
            You are a chemistry expert who is given this InChI '{inchi}' and this SMILES '{smiles}' to analyse.
            Make sure to answer each one of the bellow questions. 
            To be able to do that you are gonna need to take into account all the 
            organic groups known in chemistry, 
            the difference between a Lewis acid and base, the concept of a 
            hydrogen bond, a hydrogen bond donor
            and hydrogen bond acceptor. Once you have all this 
            information gathered, you will be able to answer. You can use
            the data collected from PubChem bellow as reference too.

            PubChem Description: 
            {ans}

            QUESTIONS:
            First, describe the molecule with this InChI in detail.
            
            Then, answer the following questions about this molecule with 
            a detailed explanation of the answer:
            
               - Is it a Lewis acid or base or both? 
               - Can it do hydrogen bonds? 
               - Is it a hydrogen bond donor or acceptor?
            """
    query = textwrap.dedent(query)

    prompt = ChatPromptTemplate.from_template(query)

    # pylint: disable=E1131
    chain = prompt | llm
    # pylint: enable=E1131
    response = chain.invoke({"smiles": smiles, "inchi": inchi, "ans": ans})

    if isinstance(response.content, str):
        return markdown(response.content)
    return response.content


def pubchem_description(inchi: str) -> dict[str, str]:
    """
    Check if the molecule is in PubChem and return its description.

    Args:
        inchi (str): The InChI of the molecule.
    """
    url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchi/description/json?inchi="
        + quote(inchi, safe="")
    )
    try:
        with urlopen(url) as ans:
            ans = loads(ans.read().decode("utf8").strip())
    except (TypeError, HTTPError, ValueError):
        ans = {"result": "no data available on this molecule in PubChem."}
    return ans


def is_api_key_valid(api_key: str) -> bool:
    "Check if the API key is valid."
    url = f"https://generativelanguage.googleapis.com/v1/models?key={quote(api_key)}"
    try:
        with urlopen(url) as ans:
            ans = ans.read().decode("utf8").rstrip()
            return True
    except HTTPError as e:
        logger.error(e)
        return False
