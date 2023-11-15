"request handler."
import re

import torch
from django.shortcuts import render

from data.graph import from_InChI, from_smiles
from train.models import PNAPCSAFTDUMMY, PnaconvsParams, ReadoutMLPParams
from train.utils import calc_deg

from .forms import InChIorSMILESinput

deg = calc_deg("ramirez", "/workspaces/ePC-SAFT")
device = torch.device("cpu")
model = PNAPCSAFTDUMMY(
    hidden_dim=128,
    pna_params=PnaconvsParams(
        propagation_depth=2,
        pre_layers=1,
        post_layers=3,
        deg=deg,
    ),
    mlp_params=ReadoutMLPParams(num_mlp_layers=1, num_para=3),
)
model.to("cpu")

checkpoint = torch.load(
    "/workspaces/ePC-SAFT/train/checkpoints/model5-0_50e6.pth", map_location="cpu"
)

model.load_state_dict(checkpoint["model_state_dict"])

model.eval()


def prediction(query: str) -> tuple[torch.Tensor, bool]:
    "Predict ePC-SAFT parameters."

    inchi_check = re.search("^InChI=", query)

    if inchi_check:
        gh_fn = from_InChI
    else:
        gh_fn = from_smiles

    try:
        graph = gh_fn(query).to(device)
        with torch.no_grad():
            pred = model.forward(graph)[0]
        output = True
    except (ValueError, TypeError, AttributeError, IndexError) as e:
        print(e)
        print("error for query:", query)
        pred = [None]
        output = False

    return pred, output


def index(request):
    "handle request"
    available_params = [
        "Segment number",
        "Segment diameter (Å)",
        "Dispersion energy (K)",
    ]
    pred = None
    query = ""
    output = False
    if request.method == "POST":
        form = InChIorSMILESinput(request.POST)

        if form.is_valid():
            query = form.cleaned_data["query"]
            pred, output = prediction(query)

    else:
        form = InChIorSMILESinput()

    context = {
        "form": form,
        "predicted_para": [
            (paraname, round(para.item(), 2))
            for para, paraname in zip(pred, available_params)
        ]
        if output
        else [(None, None)],
        "query": query,
        "output": output,
    }

    return render(request, "pred.html", context)


# Create your views here.
def homepage(request):
    "handle request"
    return render(request, "homepage.html")
