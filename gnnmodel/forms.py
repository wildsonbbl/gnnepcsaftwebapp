import re

from django import forms
from django.core.exceptions import ValidationError
from django.utils.translation import gettext_lazy as _

from model.data.graph import from_InChI, from_smiles


class BootstrapForm(forms.Form):
    "To add bootstrap class to django form"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        for visible in self.visible_fields():
            visible.field.widget.attrs["class"] = "form-control"


class InChIorSMILESinput(BootstrapForm):
    "Form to receive InChI/SMILES from user."
    query = forms.CharField(
        strip=True,
        empty_value="InChI or SMILES",
        required=True,
        help_text="InChI or SMILES",
        initial="CCO",
    )

    def clean_query(self):
        "check valid input."
        data = self.cleaned_data["query"]

        inchi_check = re.search("^InChI=", data)

        if inchi_check:
            gh_fn = from_InChI
        else:
            gh_fn = from_smiles

        try:
            gh_fn(data)
        except (ValueError, TypeError, AttributeError, IndexError) as e:
            raise ValidationError(_("Invalid InChI/SMILES.")) from e
        return data
