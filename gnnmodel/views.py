"request handler."

import os.path as osp

from django.conf import settings
from django.http import HttpResponse
from django.shortcuts import render

from .forms import InChIorSMILESareaInput
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
    pred_list = []
    output = False
    if request.method == "POST":
        form = InChIorSMILESareaInput(request.POST)
        if form.is_valid():
            _, smiles_list = form.cleaned_data["text_area"]

            for smiles in smiles_list:
                pred_list.append([round(para, 5) for para in get_pred(smiles)])
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
    if request.method == "POST":
        forms = init_mixture_forms(request.POST)
        post_data = process_mixture_post(forms)
        context = build_mixture_context(post_data)
    else:
        context = build_mixture_context()
    return render(request, "mixture.html", context)


def about(request):
    "handle request for about page"
    return render(request, "about.html")


def service_worker(request):  # pylint: disable=unused-argument
    "serve root service worker"
    with open(
        settings.STATIC_ROOT / "js/serviceworker.js", encoding="utf-8"
    ) as serviceworker_file:
        return HttpResponse(
            serviceworker_file.read(),
            content_type="application/javascript",
        )
