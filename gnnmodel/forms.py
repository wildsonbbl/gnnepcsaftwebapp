"Django forms."

import re

from django import forms
from django.core.exceptions import ValidationError
from django.utils.translation import gettext_lazy as _
from gnnepcsaft.data.ogb_utils import smiles2graph
from gnnepcsaft.data.rdkit_util import inchitosmiles, smilestoinchi


class InChIorSMILESinput(forms.Form):
    "Form to receive InChI/SMILES from user."

    query = forms.CharField(
        label="Type/Paste InChI or SMILES",
        strip=True,
        empty_value="InChI or SMILES",
        required=True,
        initial="CCO",
        widget=forms.TextInput(
            attrs={"class": "form-control", "aria-label": "Type/Paste InChI or SMILES"}
        ),
    )

    def clean_query(self):
        "check valid input."
        data = self.cleaned_data["query"]

        inchi_check = re.search("^InChI=", data)

        if inchi_check:
            try:
                data = inchitosmiles(data, False, False)
            except ValueError as e:
                raise ValidationError(_("Invalid InChI/SMILES.")) from e
        else:
            try:
                __ = smilestoinchi(data, False, False)
            except ValueError as e:
                raise ValidationError(_("Invalid InChI/SMILES.")) from e
        try:
            smiles2graph(data)
        except ValueError as e:
            raise ValidationError(_("Invalid InChI/SMILES.")) from e
        return data


class CustomPlotConfigForm(forms.Form):
    "Form to receive custom plot config."

    temp_min = forms.FloatField(
        label="Minimum Temperature (K)",
        required=False,
        initial=300.0,
        widget=forms.NumberInput(attrs={"class": "form-control"}),
    )
    temp_max = forms.FloatField(
        label="Maximum Temperature (K)",
        required=False,
        initial=400.0,
        widget=forms.NumberInput(attrs={"class": "form-control"}),
    )
    pressure = forms.FloatField(
        label="Pressure (Pa)",
        required=False,
        initial=101325.0,
        widget=forms.NumberInput(attrs={"class": "form-control"}),
    )


class CustomPlotCheckForm(forms.Form):
    "Form to check custom plot."

    custom_plot_checkbox = forms.BooleanField(
        label="Custom Plot",
        label_suffix="",
        required=False,
        initial=False,
        widget=forms.CheckboxInput(
            attrs={
                "class": "form-check-input",
                "aria-label": "Custom Plot",
            }
        ),
    )


class RhoCheckForm(forms.Form):
    "Form to check density."

    rho_checkbox = forms.BooleanField(
        label="Density (mol / m³)",
        label_suffix="",
        required=False,
        initial=True,
        widget=forms.CheckboxInput(
            attrs={
                "class": "form-check-input",
                "aria-label": "Density (mol / m³)",
            }
        ),
    )


class VPCheckForm(forms.Form):
    "Form to check vapor pressure."

    vp_checkbox = forms.BooleanField(
        label="Vapor pressure (Pa)",
        label_suffix="",
        required=False,
        initial=True,
        widget=forms.CheckboxInput(
            attrs={
                "class": "form-check-input",
                "aria-label": "Vapor pressure (Pa)",
            }
        ),
    )
