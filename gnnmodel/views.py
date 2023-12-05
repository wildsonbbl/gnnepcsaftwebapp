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
from .models import GnnepcsaftPara
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

checkpoint = torch.load(file_dir + "/static/model.pth", map_location="cpu")

model.load_state_dict(checkpoint["model_state_dict"])

model.eval()


def prediction(query: str) -> tuple[torch.Tensor, bool, str]:
    "Predict ePC-SAFT parameters."

    inchi = checking_inchi(query)

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


def checking_inchi(query: str) -> str:
    "Check if query is inchi and return an inchi."
    inchi_check = re.search("^InChI=", query)
    inchi = query
    if not inchi_check:
        try:
            inchi = smilestoinchi(query)
        except (ValueError, TypeError, AttributeError, IndexError) as e:
            print(e)
            print("error for query:", query)
    return inchi


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
            # pylint: disable=no-member
            comp = GnnepcsaftPara.objects.filter(inchi=inchi).all()
            # pylint: enable=no-member
            db_update(pred, inchi, comp)

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


def db_update(pred, inchi, comp):
    "Updates the gnnepcsaft db."
    if len(comp) == 0:
        new_comp = GnnepcsaftPara(
            inchi=inchi, m=pred[0], sigma=pred[1], e=pred[2], counting=1
        )
        new_comp.save()
    else:
        stored_comp = comp[0]
        stored_comp.counting += 1
        stored_comp.save()


# Create your views here.
def homepage(request):
    "handle request"
    return render(request, "homepage.html")


def authorpage(request):
    "handle request"
    return render(request, "author.html")
