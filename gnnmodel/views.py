"request handler."

import json
import os
import os.path as osp

from django.conf import settings
from django.http import JsonResponse
from django.shortcuts import render
from django.views.decorators.http import require_http_methods

from .forms import GoogleAPIKeyForm, InChIorSMILESareaInput, InChIorSMILESinput
from .models import ChatSession
from .utils import (
    available_params,
    build_mixture_context,
    build_pure_context,
    get_pred,
    init_mixture_forms,
    init_pure_forms,
    process_mixture_post,
    process_pure_post,
)
from .utils_llm import is_api_key_valid, resume_mol

file_dir = osp.dirname(__file__)


# Create your views here.
def pure(request):
    "handle request for pure substance"

    if request.method == "POST":
        forms = init_pure_forms(request.POST)
        post_data = process_pure_post(forms)
        context = build_pure_context(forms, post_data)
    else:
        forms = init_pure_forms()
        context = build_pure_context(forms)
    return render(request, "pure.html", context)


def batch(request):
    "handle request for batch of substances"
    if settings.PLATFORM == "webapp":
        return render(request, "batch-webapp.html")
    pred_list = []
    output = False
    if request.method == "POST":
        form = InChIorSMILESareaInput(request.POST)
        if form.is_valid():
            inchi_list, smiles_list = form.cleaned_data["text_area"]

            for smiles, inchi in zip(smiles_list, inchi_list):
                pred_list.append([round(para, 5) for para in get_pred(smiles, inchi)])
            output = True
    else:
        form = InChIorSMILESareaInput()

    context = {
        "form": form,
        "available_params": available_params,
        "parameters_list": pred_list,
        "output": output,
    }

    return render(request, "batch.html", context)


def mixture(request):
    "handle request for mixture"
    if settings.PLATFORM == "webapp":
        return render(request, "mixture-webapp.html")
    if request.method == "POST":
        forms = init_mixture_forms(request.POST)
        post_data = process_mixture_post(forms)
        context = build_mixture_context(post_data)
    else:
        context = build_mixture_context()
    return render(request, "mixture.html", context)


def homepage(request):
    "handle request"
    return render(request, "homepage.html")


def authorpage(request):
    "handle request"
    return render(request, "author.html")


def about(request):
    "handle request for about page"
    return render(request, "about.html")


def description(request):
    "handle request for molecule description"

    html_output = ""

    if request.method == "POST":
        form = InChIorSMILESinput(request.POST)
        google_api_key_form = GoogleAPIKeyForm(request.POST)
        if google_api_key_form.is_valid() and form.is_valid():
            smiles, inchi = form.cleaned_data["query"]
            html_output = resume_mol(
                inchi, smiles, google_api_key_form.cleaned_data["google_api_key"]
            )
    else:
        form = InChIorSMILESinput()
        google_api_key_form = GoogleAPIKeyForm()

    return render(
        request,
        "description.html",
        {
            "output": html_output,
            "google_api_key_form": google_api_key_form,
        },
    )


def chat(request):
    "handle request for chat"

    if settings.PLATFORM == "webapp":
        return render(request, "chat-webapp.html")

    show_form = True
    if request.method == "POST":
        google_api_key_form = GoogleAPIKeyForm(request.POST)
        if google_api_key_form.is_valid():
            show_form = False
    else:
        google_api_key_form = GoogleAPIKeyForm()
        google_api_key = os.getenv("GOOGLE_API_KEY")
        if google_api_key is not None and is_api_key_valid(google_api_key):
            show_form = False

    return render(
        request,
        "chat.html",
        {"show_form": show_form, "google_api_key_form": google_api_key_form},
    )


@require_http_methods(["GET"])
def get_sessions(request):  # pylint: disable=unused-argument
    """Get all sessions"""
    sessions = list(
        ChatSession.objects.values("session_id", "name", "created_at", "updated_at")
    )
    return JsonResponse({"sessions": sessions})


@require_http_methods(["POST"])
def create_session(request):
    """Create a new session"""
    data = json.loads(request.body)
    name = data.get("name", "New Session")
    session = ChatSession.objects.create(name=name)
    return JsonResponse({"session_id": str(session.session_id), "name": session.name})


@require_http_methods(["DELETE"])
def delete_session(request, session_id):  # pylint: disable=unused-argument
    """Delete a session"""
    try:
        session = ChatSession.objects.get(session_id=session_id)
        session.delete()
        return JsonResponse({"success": True})
    except ChatSession.DoesNotExist:
        return JsonResponse(
            {"success": False, "error": "Session not found"}, status=404
        )
