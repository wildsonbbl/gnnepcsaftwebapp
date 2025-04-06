"Module for utils to work with LLMs"

import textwrap
from urllib.parse import quote
from urllib.request import HTTPError, urlopen

from langchain_core.prompts import ChatPromptTemplate
from langchain_google_genai import ChatGoogleGenerativeAI


def resume_mol(inchi: str):
    "Describe the molecule with google's gemini."

    llm = ChatGoogleGenerativeAI(
        model="gemini-pro"
    )  # needs GOOGLE_API_KEY env variable

    url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchi/description/json?inchi="
        + quote(inchi, safe="")
    )
    try:
        with urlopen(url) as ans:
            ans = ans.read().decode("utf8").rstrip()
    except (TypeError, HTTPError, ValueError):
        ans = "no data available."

    query = """
            You are a chemistry expert who is given this InChI {inchi} to analyse.
            Make sure to answer each one of the bellow questions. 
            To be able to do that you are gonna need to take into account all the 
            organic groups known in chemistry, 
            the difference between a Lewis acid and base, the concept of a 
            hydrogen bond, a hydrogen bond donor
            and hydrogen bond acceptor. Once you have all this 
            information gathered, you will be able to answer. You can use
            the passage bellow as reference.

            QUESTIONS: '
            First, describe the molecule with this InChI in detail.
            
            Then, answer the following questions about this molecule with 
            a detailed explanation of the answer:
            
               - Is it a Lewis acid or base or both? 
               - Can it do hydrogen bonds? 
               - Is it a hydrogen bond donor or acceptor?'
               
            PASSAGE: '{ans}'

            """
    query = textwrap.dedent(query)

    prompt = ChatPromptTemplate.from_template(query)

    # pylint: disable=E1131
    chain = prompt | llm
    # pylint: enable=E1131
    response = chain.invoke({"inchi": inchi, "ans": ans})

    return response.content
