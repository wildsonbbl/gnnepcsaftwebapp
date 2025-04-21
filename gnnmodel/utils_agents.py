"""Module for LLM agents"""

import inspect
import io
import os
import re
import textwrap
from ast import literal_eval
from contextlib import redirect_stdout
from json import dumps
from typing import Any, Callable, Dict, List, Literal, Optional, Union

from gnnepcsaft.data.rdkit_util import mw, smilestoinchi
from gnnepcsaft.epcsaft.epcsaft_feos import (
    mix_den_feos,
    mix_vp_feos,
    pure_den_feos,
    pure_vp_feos,
)
from langchain_community.tools.tavily_search import TavilySearchResults
from langchain_core.language_models import BaseChatModel
from langchain_core.messages import BaseMessage, HumanMessage
from langchain_core.rate_limiters import InMemoryRateLimiter
from langchain_core.runnables import RunnableConfig
from langchain_core.tools import BaseTool, tool
from langchain_google_genai import ChatGoogleGenerativeAI
from langgraph.checkpoint.memory import MemorySaver
from langgraph.prebuilt import create_react_agent
from pydantic import SecretStr, ValidationError

from .utils import predict_epcsaft_parameters
from .utils_llm import pubchem_description


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
            You are given this InChI '{inchi}' and this SMILES '{smiles}' to analyse. 
            Find out what it is and search more info on it.
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


def custom_agent(llm: BaseChatModel, message: str, tools: Dict[str, BaseTool]):
    "Agent to use functions"

    python_docs = dumps(
        [tools[key].args_schema.model_json_schema() for key in tools],  # type: ignore
    )

    instruction_prompt_with_function_calling = (
        textwrap.dedent(
            """
    At each turn, if you decide to call any of the function(s), 
    you MUST put it in the format: 
    
    ```tool_code
    [{'tool_name': 'func_name1', 'arguments': {'arg1': value1, 'arg2': value2, ...}},
    {'tool_name': 'func_name2', 'arguments': {'arg1': value1, 'arg2': value2, ...}}, ...]
    ```
    
    You SHOULD NOT include any other text in the response if you call a function.
    
    After calling a functions(s), just wait the next turn to receive the response with the format: 
    
    ```tool_output 
    [func_name1 response, func_name2 response, ...]
    ``` 
    
    Use it to call more tools or generate a helpful, friendly response. 

    The python methods described below are imported and available, 
    you can only use defined methods. The generated code should be 
    readable and efficient.

    The following Python methods are available:
    """
        )
        + f"""

    ```python\n{python_docs}
    ```

    USER:\n{message}
    """
    )
    human_message = HumanMessage(content=instruction_prompt_with_function_calling)
    human_message.pretty_print()
    messages = []
    messages += [human_message]

    messages = agent_loop(llm, messages, tools)
    return messages


def agent_loop(
    llm: BaseChatModel,
    messages: List[Union[BaseMessage, HumanMessage]],
    tools: Dict[str, BaseTool],
):
    "Loop of the agent to call the functions"
    llm_response = llm.invoke(messages)
    messages += [llm_response]
    llm_response.pretty_print()
    assert isinstance(llm_response.content, str)
    call_response = extract_tool_call(llm_response.content, tools)
    while call_response is not None:
        tool_message = HumanMessage(content=call_response)
        tool_message.pretty_print()
        messages += [tool_message]
        llm_response = llm.invoke(messages)
        messages += [llm_response]
        llm_response.pretty_print()
        assert isinstance(llm_response.content, str)
        call_response = extract_tool_call(llm_response.content, tools)
    return messages


def reviewer_agent(
    llm: BaseChatModel, messages: List[BaseMessage], tools: Dict[str, BaseTool]
):
    """
    Review the response of the previous agent
    """
    human_message = HumanMessage(
        content="This is all the work of another agent. \
          Make a review of it and correct it if necessary."
    )
    messages += [human_message]
    messages = agent_loop(llm, messages, tools)
    return messages


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


def extract_tool_call(text: str, tools: Dict[str, BaseTool]) -> Optional[str]:
    "extract tools call from the agent response."

    pattern = r"```tool_code\s*(.*?)\s*```"
    match = re.search(pattern, text, re.DOTALL)
    if match:
        code = match.group(1).strip()
        # Capture stdout in a string buffer
        f = io.StringIO()
        with redirect_stdout(f):
            try:
                tool_calls = literal_eval(code)
                result = []
                for tool_call in tool_calls:
                    try:
                        tool_name = tool_call["tool_name"]
                        args = tool_call["arguments"]
                        fn_result = tools[tool_name].invoke(input=args)
                        result.append(fn_result)
                    except (
                        RuntimeError,
                        SyntaxError,
                        AssertionError,
                        TypeError,
                        IndexError,
                        ValidationError,
                        ValueError,
                    ) as e:
                        result.append(
                            f"ERROR on this function call: {e}."
                            f" You MUST fix the issues for next turn."
                        )
            except (
                RuntimeError,
                SyntaxError,
                AssertionError,
                TypeError,
                ValueError,
            ) as e:
                result = (
                    f"ERROR on tool_code: {e}."
                    f" You MUST fix the issues for next turn."
                    f" Make sure to write only python code in ```toll_code ```."
                )
        output = f.getvalue()
        r = result if output == "" else output
        return f"```tool_output\n{str(r)}\n```"
    return None


def langgraph_agent(
    _prompt: str, llm: BaseChatModel, fn_list: List[Callable[..., Any]]
):
    """Agent to check PCSAFT functions"""

    memory = MemorySaver()
    agent = create_react_agent(
        llm,
        tools=[tool(fn, parse_docstring=True) for fn in fn_list],
        checkpointer=memory,
    )
    query = textwrap.dedent(_prompt)
    config = RunnableConfig(configurable={"thread_id": "abc123"})
    step = {"messages": [HumanMessage(content=query)]}
    for step in agent.stream(
        {"messages": [HumanMessage(content=query)]},
        config,
        stream_mode="values",
    ):
        step["messages"][-1].pretty_print()
    return step


if __name__ == "__main__":

    PROMPT = """
    I have two molecules with SMILES 'CCO' and 'CC(=0)C'. 
    Consider them pure at 325 K and 101325 Pa.
    - What's their density? Which one is more dense?
    - What's their vapor pressure? Which one is more volatile?
    - What's their pure component phase.

    Do just what's asked.
"""

    PROMPT2 = """
    Then consider their 50/50 mixture at 300 K and 101325 Pa.
    - What's their density?
    - What's the mixture phase.
"""

    PROMPT3 = """
    What's the molecule with those SMILES?
    If there's any source info about it, what does each source say?
    """

    from langchain_groq import ChatGroq  # pylint: disable=import-error
    from langchain_mistralai.chat_models import (  # pylint: disable=import-error
        ChatMistralAI,
    )

    _fn_list = [
        pure_vp_feos,
        pure_den_feos,
        mix_den_feos,
        mix_vp_feos,
        pure_phase,
        mixture_phase,
        pubchem_description,
        mw,
        smilestoinchi,
        predict_epcsaft_parameters,
    ]

    _tools = {fn.__name__: tool(fn, parse_docstring=True) for fn in _fn_list}

    llama3_70b = ChatGroq(
        model="llama3-70b-8192",
        rate_limiter=InMemoryRateLimiter(requests_per_second=25 / 60),
    )

    llama4_maverick = ChatGroq(
        model="meta-llama/llama-4-maverick-17b-128e-instruct",
        rate_limiter=InMemoryRateLimiter(requests_per_second=25 / 60),
    )
    llama4_scout = ChatGroq(
        model="meta-llama/llama-4-scout-17b-16e-instruct",
        rate_limiter=InMemoryRateLimiter(requests_per_second=25 / 60),
    )

    llama33_70b = ChatGroq(
        model="llama-3.3-70b-versatile",
        rate_limiter=InMemoryRateLimiter(requests_per_second=25 / 60),
    )

    gemini20_flash = ChatGoogleGenerativeAI(
        model="gemini-2.0-flash",
        api_key=SecretStr(os.environ.get("GEMINI_API_KEY", "")),
        rate_limiter=InMemoryRateLimiter(requests_per_second=5 / 60),
    )

    gemma3_27b_it = ChatGoogleGenerativeAI(
        model="gemma-3-27b-it",
        api_key=SecretStr(os.environ.get("GEMINI_API_KEY", "")),
        rate_limiter=InMemoryRateLimiter(requests_per_second=25 / 60),
    )

    deepseek_r1 = ChatGroq(
        model="deepseek-r1-distill-llama-70b",
        rate_limiter=InMemoryRateLimiter(requests_per_second=25 / 60),
    )

    mistral = ChatMistralAI(
        model_name="mistral-large-2411",
        rate_limiter=InMemoryRateLimiter(requests_per_second=1),
    )

    # langgraph_agent(PROMPT, gemini20_flash, _fn_list)

    # agent_messages = custom_agent(mistral, PROMPT, _tools)
    # prompt_message = HumanMessage(content=PROMPT2)
    # prompt_message.pretty_print()
    # agent_messages += [prompt_message]
    # agent_messages = agent_loop(llama4_maverick, agent_messages, _tools)
    # prompt_message = HumanMessage(content=PROMPT3)
    # prompt_message.pretty_print()
    # agent_messages += [prompt_message]
    # agent_messages = agent_loop(gemma3_27b_it, agent_messages, _tools)
    # reviewer_messages = reviewer_agent(gemini20_flash, agent_messages, _tools)
