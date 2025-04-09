"""Module for LLM agents"""

import inspect
import io
import os
import re
import textwrap
import traceback
from contextlib import redirect_stdout
from typing import Any, Callable, List, Literal, Optional, Union

from gnnepcsaft.epcsaft.epcsaft_feos import (
    mix_den_feos,
    mix_vp_feos,
    pure_den_feos,
    pure_vp_feos,
)
from langchain_community.tools.tavily_search import TavilySearchResults
from langchain_core.messages import HumanMessage
from langchain_core.rate_limiters import InMemoryRateLimiter
from langchain_core.runnables import RunnableConfig
from langchain_core.tools import tool
from langchain_google_genai import ChatGoogleGenerativeAI
from langchain_groq import ChatGroq
from langgraph.checkpoint.memory import MemorySaver
from langgraph.prebuilt import create_react_agent
from pydantic import SecretStr


def resume_mol_agent(inchi, smiles, llm):
    "Agent to describe the molecule."

    memory = MemorySaver()
    search = TavilySearchResults(
        max_results=2,
        include_domains=["pubchem.ncbi.nlm.nih.gov", "wikipedia.org", "chemrxiv.org"],
    )
    agent = create_react_agent(
        llm,
        tools=[search],
        checkpointer=memory,
    )
    query = f"""
            You are a chemistry expert who is given this 
            InChI '{inchi}' and 
            this SMILES '{smiles}' to analyse. 
            Make sure to answer each one of the bellow questions. 
            To be able to do that you are gonna need to take into account all the 
            organic groups known in chemistry, 
            the difference between a Lewis acid and base, the concept of a 
            hydrogen bond, a hydrogen bond donor 
            and hydrogen bond acceptor. Once you have all this 
            information gathered, you will be able to answer. 

            QUESTIONS: 
            First, from the InChI and SMILES try to come out with the name
             of the molecule, and with this name make a search to find out more about it. 
            Make a summary of what you found out.
            
            Then, answer the following questions about this molecule with 
            a detailed explanation of the answer: 
               - Is it a Lewis acid or base or both? 
               - Can it do hydrogen bonds? 
               - Is it a hydrogen bond donor or acceptor? 
               
            """
    query = textwrap.dedent(query)
    config = RunnableConfig(configurable={"thread_id": "abc123"})
    step = {"messages": [HumanMessage(content=query)]}
    for step in agent.stream(
        step,
        config,
        stream_mode="values",
    ):
        step["messages"][-1].pretty_print()


def pure_phase(
    vapor_pressure: float, system_pressure: float
) -> Literal["liquid", "vapor"]:
    """
    Given the vapor pressure and system pressure, return the phase of the molecule.
    Both pressures must be in the same unit.

    Args:
        vapor_pressure (float): The calculated vapor pressure of the pure component.
        system_pressure (float): The actual system pressure.

    """
    assert isinstance(vapor_pressure, (int, float)), "vapor_pressure must be a number"
    assert isinstance(system_pressure, (int, float)), "system_pressure must be a number"
    assert vapor_pressure > 0, "vapor_pressure must be positive"
    assert system_pressure > 0, "system_pressure must be positive"

    return "liquid" if vapor_pressure < system_pressure else "vapor"


def mixture_phase(
    bubble_point: float,
    dew_point: float,
    system_pressure: float,
) -> Literal["liquid", "vapor", "two-phase"]:
    """
    Given the bubble/dew point of the mixture and the system pressure,
    return the phase of the mixture.
    All pressures must be in the same unit.

    Args:
        bubble_point (float): The calculated bubble point of the mixture.
        dew_point (float): The calculated dew point of the mixture.
        system_pressure (float): The actual system pressure.
    """
    assert isinstance(bubble_point, (int, float)), "bubble_point must be a number"
    assert isinstance(dew_point, (int, float)), "dew_point must be a number"
    assert isinstance(system_pressure, (int, float)), "system_pressure must be a number"
    assert bubble_point > 0, "bubble_point must be positive"
    assert dew_point > 0, "dew_point must be positive"
    assert system_pressure > 0, "system_pressure must be positive"
    return (
        "liquid"
        if bubble_point < system_pressure
        else ("two-phase" if dew_point <= system_pressure else "vapor")
    )


def custom_agent(
    llm: Union[ChatGoogleGenerativeAI, ChatGroq],
    message: str,
    fn_list: Optional[List[Callable[..., Any]]] = None,
):
    "Agent to use functions"

    if fn_list is None:
        fn_list = [
            pure_vp_feos,
            pure_den_feos,
            mix_den_feos,
            mix_vp_feos,
            pure_phase,
            mixture_phase,
        ]

    python_docs = get_python_docs(fn_list)

    instruction_prompt_with_function_calling = """
    At each turn, if you decide to invoke any of the function(s), 
    you MUST put it in the format: 
    
    ```tool_code
    [func_name1(params_name1=params_value1, params_name2=params_value2...), func_name2(params), ...]
    ```
    
    You SHOULD NOT include any other text in the response if you call a function.
    
    The user will give the response to a method with the format: 
    
    ```tool_output 
    [func_name1 response, func_name2 response, ...]
    ``` 
    
    Use it to call more tools or generate a helpful, friendly response. 

    The python methods described below are imported and available, 
    you can only use defined methods. The generated code should be 
    readable and efficient.

    The following Python methods are available:

    ```python
    {python_docs}
    ```

    USER: 
    {message}
    """
    human_message = HumanMessage(
        content=textwrap.dedent(instruction_prompt_with_function_calling).format(
            python_docs=python_docs, message=textwrap.dedent(message)
        )
    )
    human_message.pretty_print()
    messages = [
        human_message,
    ]

    llm_response = llm.invoke(messages)
    messages += [llm_response]
    llm_response.pretty_print()
    call_response = extract_tool_call(llm_response.content)
    while call_response is not None:
        HumanMessage(content=call_response).pretty_print()
        messages += [HumanMessage(content=call_response)]
        llm_response = llm.invoke(messages)
        messages += [llm_response]
        llm_response.pretty_print()
        call_response = extract_tool_call(llm_response.content)


def get_python_docs(fn_list):
    "Get the python docs from the functions."
    python_docs = ""
    for fn in fn_list:
        python_docs += (
            fn.__name__
            + str(inspect.signature(fn))
            + ":\n"
            + str(inspect.getdoc(fn))
            + "\n\n"
        )

    return python_docs


def extract_tool_call(text) -> Optional[str]:
    "tool call from the response for gemma."

    pattern = r"```tool_code\s*(.*?)\s*```"
    match = re.search(pattern, text, re.DOTALL)
    if match:
        code = match.group(1).strip()
        # Capture stdout in a string buffer
        f = io.StringIO()
        with redirect_stdout(f):
            try:
                result = eval(code)  # pylint: disable=W0123
            except RuntimeError:
                result = f"Error on function call:\n\n{traceback.format_exc()}"
            except TypeError:
                result = (
                    traceback.format_exc()
                    + "\n\nCheck what's the output type of each function and try again."
                )
        output = f.getvalue()
        r = result if output == "" else output
        return f"```tool_output\n{str(r)}\n```"
    return None


def langgraph_agent(_prompt, llm, fn_list=None):
    """Agent to check PCSAFT functions"""

    if fn_list is None:
        fn_list = [
            pure_vp_feos,
            pure_den_feos,
            mix_den_feos,
            mix_vp_feos,
            pure_phase,
            mixture_phase,
        ]

    memory = MemorySaver()
    agent = create_react_agent(
        llm,
        tools=[tool(fn, parse_docstring=True) for fn in fn_list],
        checkpointer=memory,
    )
    query = textwrap.dedent(_prompt)
    config = RunnableConfig(configurable={"thread_id": "abc123"})
    for step in agent.stream(
        {"messages": [HumanMessage(content=query)]},
        config,
        stream_mode="values",
    ):
        step["messages"][-1].pretty_print()
    return agent


if __name__ == "__main__":

    PROMPT = """
       I have the PCSAFT parameters for molecule A as [2.87, 2.96, 187.37, 0.0559, 2460.62, 0.0, 1, 1]
       and for molecule B as [2.7758, 3.2416, 230.6054, 0.0002, 1448.9628, 0.0, 1, 0]. 
       Consider them pure at 325 K and 101325 Pa.
       - What's their density? Which one is more dense?
       - What's their vapor pressure? Which one is more volatile?
       After this check their pure component phase.

       Then consider their 50/50 mixture at 300 K and 101325 Pa.
       - What's their density?
       After this, check the mixture phase.

       Do just what's asked.
       """

    llama3_70b = ChatGroq(
        model="llama3-70b-8192",
        rate_limiter=InMemoryRateLimiter(requests_per_second=25 / 60),
    )

    llama4_maverick = ChatGroq(
        model="meta-llama/llama-4-maverick-17b-128e-instruct",
        rate_limiter=InMemoryRateLimiter(requests_per_second=25 / 60),
    )

    llama_33_70b = ChatGroq(
        model="llama-3.3-70b-versatile",
        rate_limiter=InMemoryRateLimiter(requests_per_second=25 / 60),
    )

    gemini_20_flash = ChatGoogleGenerativeAI(
        model="gemini-2.0-flash",
        api_key=SecretStr(os.environ.get("GEMINI_API_KEY", "")),
        rate_limiter=InMemoryRateLimiter(requests_per_second=5 / 60),
    )

    gemma_3_27b_it = ChatGoogleGenerativeAI(
        model="gemma-3-27b-it",
        api_key=SecretStr(os.environ.get("GEMINI_API_KEY", "")),
        rate_limiter=InMemoryRateLimiter(requests_per_second=25 / 60),
    )

    # langgraph_agent(PROMPT, gemini_20_flash)

    custom_agent(llama3_70b, PROMPT)
