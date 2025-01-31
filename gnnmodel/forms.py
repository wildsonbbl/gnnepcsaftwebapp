"Django forms."
import re

from django import forms
from django.core.exceptions import ValidationError
from django.utils.translation import gettext_lazy as _
from gnnepcsaft.data.ogb_utils import smiles2graph
from gnnepcsaft.data.rdkit_util import inchitosmiles


class BootstrapForm(forms.Form):
    "To add bootstrap class to django form"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        for visible in self.visible_fields():
            visible.field.widget.attrs["class"] = "form-control"


class InChIorSMILESinput(BootstrapForm):
    "Form to receive InChI/SMILES from user."
    query = forms.CharField(
        label="Type/Paste InChI or SMILES",
        strip=True,
        empty_value="InChI or SMILES",
        required=True,
        initial="CCO",
    )

    def clean_query(self):
        "check valid input."
        data = self.cleaned_data["query"]

        inchi_check = re.search("^InChI=", data)

        if inchi_check:
            data = inchitosmiles(data, False, False)
        try:
            smiles2graph(data)
        except (ValueError, TypeError, AttributeError, IndexError) as e:
            raise ValidationError(_("Invalid InChI/SMILES.")) from e
        return data
