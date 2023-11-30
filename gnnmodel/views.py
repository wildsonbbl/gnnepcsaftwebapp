"request handler."
import os.path as osp
import re

import gnnepcsaft
import torch
from django.shortcuts import render
from gnnepcsaft.data.graph import from_InChI, smilestoinchi
from gnnepcsaft.train.models import PNAPCSAFT, PnaconvsParams, ReadoutMLPParams
from gnnepcsaft.train.utils import calc_deg

from .forms import InChIorSMILESinput
from .utils import plotdata

file_dir = osp.dirname(__file__)
workdir = osp.dirname(gnnepcsaft.__file__)


deg = calc_deg("ramirez", workdir)
device = torch.device("cpu")
model = PNAPCSAFT(
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

checkpoint = torch.load(file_dir + "/static/model5-11_00e6.pth", map_location="cpu")

model.load_state_dict(checkpoint["model_state_dict"])

model.eval()


def prediction(query: str) -> tuple[torch.Tensor, bool, str]:
    "Predict ePC-SAFT parameters."

    inchi_check = re.search("^InChI=", query)
    inchi = query
    if not inchi_check:
        try:
            inchi = smilestoinchi(query)
        except (ValueError, TypeError, AttributeError, IndexError) as e:
            print(e)
            print("error for query:", query)

    try:
        graph = from_InChI(inchi).to(device)
        with torch.no_grad():
            pred = model.forward(graph)[0]
        output = True
    except (ValueError, TypeError, AttributeError, IndexError) as e:
        print(e)
        print("error for query:", query)
        pred = [None]
        output = False

    return pred, output, inchi


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
    plotden, plotvp = "", ""
    if request.method == "POST":
        form = InChIorSMILESinput(request.POST)

        if form.is_valid():
            query = form.cleaned_data["query"]
            pred, output, inchi = prediction(query)
            plotden, plotvp = plotdata(pred.numpy(), inchi)

    else:
        form = InChIorSMILESinput()
    imgtype = "image/png"
    den_uri = f"data:{imgtype};base64,{plotden}"
    vp_uri = f"data:{imgtype};base64,{plotvp}"

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
        "plotden": plotden,
        "plotvp": plotvp,
        "den_uri": den_uri,
        "vp_uri": vp_uri,
    }

    return render(request, "pred.html", context)


# Create your views here.
def homepage(request):
    "handle request"
    return render(request, "homepage.html")


def authorpage(request):
    "handle request"
    return render(request, "author.html")
