"Django forms."

import re

from django import forms
from django.core.exceptions import ValidationError
from django.utils.translation import gettext_lazy as _
from gnnepcsaft.data.ogb_utils import smiles2graph
from gnnepcsaft.data.rdkit_util import assoc_number, inchitosmiles, smilestoinchi


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
        "check valid input and output SMILES."
        data = self.cleaned_data["query"]

        inchi_check = re.search("^InChI=", data)

        if inchi_check:
            try:
                smiles = inchitosmiles(data, False, False)
                inchi = smilestoinchi(smiles, False, False)
                smiles2graph(smiles)
                assoc_number(inchi)
            except ValueError as e:
                raise ValidationError(_("Invalid InChI/SMILES.")) from e
        else:
            try:
                inchi = smilestoinchi(data, False, False)
                smiles = inchitosmiles(inchi, False, False)
                smiles2graph(smiles)
                assoc_number(inchi)
            except ValueError as e:
                raise ValidationError(_("Invalid InChI/SMILES.")) from e

        return smiles, inchi


class CustomPlotConfigForm(forms.Form):
    "Form to receive custom plot config."

    temp_min = forms.FloatField(
        label="Minimum Temperature (K)",
        min_value=1.0,
        required=False,
        initial=300.0,
        widget=forms.NumberInput(attrs={"class": "form-control"}),
    )
    temp_max = forms.FloatField(
        label="Maximum Temperature (K)",
        min_value=1.0,
        required=False,
        initial=400.0,
        widget=forms.NumberInput(attrs={"class": "form-control"}),
    )
    pressure = forms.FloatField(
        label="Pressure (Pa)",
        min_value=1.0,
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


class HlvCheckForm(forms.Form):
    "Form to check enthalpy of vaporization."

    h_lv_checkbox = forms.BooleanField(
        label="Enthalpy of vaporization (kJ/mol)",
        label_suffix="",
        required=False,
        initial=False,
        widget=forms.CheckboxInput(
            attrs={
                "class": "form-check-input",
                "aria-label": "Enthalpy of vaporization (kJ/mol)",
            }
        ),
    )


class SlvCheckForm(forms.Form):
    "Form to check entropy of vaporization."

    s_lv_checkbox = forms.BooleanField(
        label="Entropy of vaporization (J/mol/K)",
        label_suffix="",
        required=False,
        initial=False,
        widget=forms.CheckboxInput(
            attrs={
                "class": "form-check-input",
                "aria-label": "Entropy of vaporization (J/mol/K)",
            }
        ),
    )


class PhaseDiagramCheckForm(forms.Form):
    "Form to check phase diagram."

    phase_diagram_checkbox = forms.BooleanField(
        label="Phase diagrams",
        label_suffix="",
        required=False,
        initial=False,
        widget=forms.CheckboxInput(
            attrs={
                "class": "form-check-input",
                "aria-label": "Phase diagrams",
            }
        ),
    )
