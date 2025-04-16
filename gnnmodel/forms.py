"Django forms."

import os
import re

from django import forms
from django.core.exceptions import ValidationError
from django.utils.translation import gettext_lazy as _
from gnnepcsaft.data.ogb_utils import smiles2graph
from gnnepcsaft.data.rdkit_util import assoc_number, inchitosmiles, smilestoinchi
from pydantic import SecretStr

from .utils_llm import is_api_key_valid


class InChIorSMILESinput(forms.Form):
    "Form to receive InChI/SMILES from user."

    query = forms.CharField(
        strip=True,
        required=True,
        widget=forms.TextInput(
            attrs={
                "class": "form-control",
                "aria-label": "Type/Paste InChI or SMILES",
                "placeholder": "Type/Paste InChI or SMILES",
            }
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


class InChIorSMILESareaInput(forms.Form):
    "Form to receive InChI/SMILES from user."

    text_area = forms.CharField(
        label="Type/Paste a list of InChI or SMILES",
        strip=True,
        required=True,
        widget=forms.Textarea(
            attrs={
                "class": "form-control my-2",
                "aria-label": "Type/Paste InChI or SMILES",
                "placeholder": "One InChI or SMILES per line, example:"
                "\n\nCCO\nCC\nO\nCC(O)C",
            }
        ),
    )

    def clean_text_area(self):
        "check valid input and output SMILES."
        data: str = self.cleaned_data["text_area"]

        lines = data.split("\n")
        inchi_list, smiles_list = [], []
        for line in lines:
            line = line.strip()
            inchi_check = re.search("^InChI=", line)
            if inchi_check:
                try:
                    smiles = inchitosmiles(line, False, False)
                    inchi = smilestoinchi(smiles, False, False)
                    smiles2graph(smiles)
                    assoc_number(inchi)
                    smiles_list.append(smiles)
                    inchi_list.append(inchi)
                except ValueError as e:
                    raise ValidationError(_(f"Invalid InChI/SMILES: {line}")) from e
            else:
                try:
                    inchi = smilestoinchi(line, False, False)
                    smiles = inchitosmiles(inchi, False, False)
                    smiles2graph(smiles)
                    assoc_number(inchi)
                    smiles_list.append(smiles)
                    inchi_list.append(inchi)
                except ValueError as e:
                    raise ValidationError(_(f"Invalid InChI/SMILES: {line}")) from e
        return inchi_list, smiles_list


class InChIorSMILESareaInputforMixture(forms.Form):
    "Form to receive InChI/SMILES + mole fractions from user."

    text_area = forms.CharField(
        label="Type/Paste a list of InChI/SMILES | mole fractions",
        strip=True,
        required=True,
        widget=forms.Textarea(
            attrs={
                "class": "form-control my-2",
                "aria-label": "Type/Paste InChI or SMILES",
                "placeholder": "One 'InChI/SMILES Mole Fraction' per line,"
                " example:\n\nCCO 0.33\nCC 0.33\nInChI=1S/C3H8/c1-3-2/h3H2,1-2H3 0.33",
            }
        ),
    )

    def clean_text_area(self):
        "check valid input and output SMILES."
        data: str = self.cleaned_data["text_area"]

        lines = data.split("\n")
        inchi_list, smiles_list = [], []
        mole_fraction_list = []
        for full_line in lines:
            try:
                query, mole_fraction = full_line.strip().split(" ", maxsplit=1)
                query = query.strip()
                mole_fraction = mole_fraction.strip()
            except ValueError as e:
                raise ValidationError(_(f'Missing " " for line: {full_line}')) from e

            inchi_check = re.search("^InChI=", query)
            if inchi_check:
                try:
                    smiles = inchitosmiles(query, False, False)
                    inchi = smilestoinchi(smiles, False, False)
                    smiles2graph(smiles)
                    assoc_number(inchi)
                    smiles_list.append(smiles)
                    inchi_list.append(inchi)
                except ValueError as e:
                    raise ValidationError(
                        _(f'Invalid InChI/SMILES: "{query}" in line "{full_line}"')
                    ) from e
            else:
                try:
                    inchi = smilestoinchi(query, False, False)
                    smiles = inchitosmiles(inchi, False, False)
                    smiles2graph(smiles)
                    assoc_number(inchi)
                    smiles_list.append(smiles)
                    inchi_list.append(inchi)
                except ValueError as e:
                    raise ValidationError(
                        _(f'Invalid InChI/SMILES: "{query}" in line "{full_line}"')
                    ) from e
            try:
                mole_fraction_list.append(float(mole_fraction))
            except ValueError as e:
                raise ValidationError(
                    _(f'Invalid Mole Fraction: "{mole_fraction}" in line "{full_line}"')
                ) from e
        return inchi_list, smiles_list, mole_fraction_list


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


class STCheckForm(forms.Form):
    "Form to check surface tension."

    st_checkbox = forms.BooleanField(
        label="Surface tension (mN/m)",
        label_suffix="",
        required=False,
        initial=False,
        widget=forms.CheckboxInput(
            attrs={
                "class": "form-check-input",
                "aria-label": "Surface tension (mN/m)",
            }
        ),
    )


class GoogleAPIKeyForm(forms.Form):
    "Form to receive Google API Key."

    google_api_key = forms.CharField(
        strip=True,
        required=True,
        empty_value="",
        widget=forms.PasswordInput(
            attrs={
                "class": "form-control",
                "aria-label": "Gemini API Key",
                "placeholder": "Paste your Gemini API key here or set env variable GOOGLE_API_KEY",
            }
        ),
    )

    def clean_google_api_key(self):
        "check valid Google API Key."
        google_api_key = self.cleaned_data["google_api_key"]
        if google_api_key:
            if not is_api_key_valid(google_api_key):
                raise ValidationError(_("Invalid Gemini API Key"))
        elif os.environ.get("GOOGLE_API_KEY") is None:
            raise ValidationError(
                _("Gemini API Key is required for AI generated content")
            )
        elif not is_api_key_valid(os.environ.get("GOOGLE_API_KEY", "")):
            raise ValidationError(_("Invalid Gemini API Key"))
        os.environ["GOOGLE_API_KEY"] = google_api_key
        return SecretStr(google_api_key)
