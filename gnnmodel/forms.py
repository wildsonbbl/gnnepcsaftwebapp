"Django forms."

import re

from django import forms
from django.conf import settings
from django.core.exceptions import ValidationError
from django.utils.translation import gettext_lazy as _
from gnnepcsaft.data.ogb_utils import smiles2graph
from gnnepcsaft.data.rdkit_util import assoc_number, inchitosmiles, smilestoinchi


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
        if settings.PLATFORM == "webapp" and len(lines) > 10:
            raise ValidationError(_("Maximum 10 substances allowed for web app."))
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
                " example:\n\nCCO 0.33\nCC 0.33\nInChI=1S/C3H8/c1-3-2/h3H2,1-2H3 0.33\n"
                "k12 k13 k23\n\n(Note: last line is kij values)",
            }
        ),
    )

    def clean_text_area(self):
        "check valid input and output SMILES."
        data: str = self.cleaned_data["text_area"]

        lines = data.split("\n")
        kij = lines.pop()
        kij = kij.strip().split(" ")
        if len(kij) != (len(lines) ** 2 - len(lines)) / 2:
            raise ValidationError(
                _(
                    f"Number of Kij values ({len(kij)}) must "
                    f"be equal to {int((len(lines) ** 2 - len(lines)) / 2)}."
                )
            )
        # check if all kij are float
        try:
            kij = [float(i) for i in kij]
        except ValueError as e:
            raise ValidationError(
                _("All Kij values in last line must be valid numbers.")
            ) from e
        for k in kij:
            if k < -1.0 or k > 1.0:
                raise ValidationError(_("Kij values must be between -1 and 1."))

        if settings.PLATFORM == "webapp" and len(lines) > 10:
            raise ValidationError(_("Maximum 10 components allowed for web app."))
        inchi_list, smiles_list = [], []
        mole_fraction_list = []
        for full_line in lines:
            try:
                query, mole_fraction = full_line.strip().split(" ", maxsplit=1)
                query = query.strip()
                mole_fraction = mole_fraction.strip()
            except ValueError as e:
                raise ValidationError(
                    _(f"Missing mole fraction for line: {full_line}")
                ) from e

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
        return inchi_list, smiles_list, mole_fraction_list, kij


class CustomPlotConfigForm(forms.Form):
    "Form to receive custom plot config."

    temp_min = forms.FloatField(
        label="",
        min_value=0.0,
        required=False,
        widget=forms.NumberInput(
            attrs={
                "class": "form-control m-1",
                "aria-label": "Minimum Temperature (K)",
                "placeholder": "Minimum Temperature (K)",
            }
        ),
    )

    temp_max = forms.FloatField(
        label="",
        min_value=0.0,
        required=False,
        widget=forms.NumberInput(
            attrs={
                "class": "form-control m-1",
                "aria-label": "Maximum Temperature (K)",
                "placeholder": "Maximum Temperature (K)",
            }
        ),
    )

    pressure = forms.FloatField(
        label="",
        min_value=0.0,
        required=False,
        widget=forms.NumberInput(
            attrs={
                "class": "form-control m-1",
                "aria-label": "Pressure (Pa)",
                "placeholder": "Pressure (Pa)",
            }
        ),
    )

    npoints = forms.IntegerField(
        label="",
        min_value=2,
        required=False,
        widget=forms.NumberInput(
            attrs={
                "class": "form-control m-1",
                "aria-label": "Number of data points",
                "placeholder": "Number of data points",
            }
        ),
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
        label="Density (rho-T)",
        label_suffix="",
        required=False,
        initial=True,
        widget=forms.CheckboxInput(
            attrs={
                "class": "form-check-input",
                "aria-label": "Density (rho-T)",
            }
        ),
    )


class VPCheckForm(forms.Form):
    "Form to check vapor pressure."

    vp_checkbox = forms.BooleanField(
        label="Vap. Pres. (P-T)",
        label_suffix="",
        required=False,
        initial=True,
        widget=forms.CheckboxInput(
            attrs={
                "class": "form-check-input",
                "aria-label": "Vap. Pres. (P-T)",
            }
        ),
    )


class HlvCheckForm(forms.Form):
    "Form to check enthalpy of vaporization."

    h_lv_checkbox = forms.BooleanField(
        label="Enthalpy of Vap. (H-T)",
        label_suffix="",
        required=False,
        initial=False,
        widget=forms.CheckboxInput(
            attrs={
                "class": "form-check-input",
                "aria-label": "Enthalpy of Vap H-T)",
            }
        ),
    )


class SlvCheckForm(forms.Form):
    "Form to check entropy of vaporization."

    s_lv_checkbox = forms.BooleanField(
        label="Entropy of Vap. (S-T)",
        label_suffix="",
        required=False,
        initial=False,
        widget=forms.CheckboxInput(
            attrs={
                "class": "form-check-input",
                "aria-label": "Entropy of Vap. (S-T)",
            }
        ),
    )


class TernaryLLECheckForm(forms.Form):
    "Form to check liquid-liquid equilibrium."

    ternary_lle_checkbox = forms.BooleanField(
        label="Ternary VLE/LLE (x1-x2)",
        label_suffix="",
        required=False,
        initial=False,
        widget=forms.CheckboxInput(
            attrs={
                "class": "form-check-input",
                "aria-label": "Ternary VLE/LLE (x1-x2)",
            }
        ),
    )


class BinaryLLECheckForm(forms.Form):
    "Form to check liquid-liquid equilibrium."

    binary_lle_checkbox = forms.BooleanField(
        label="Binary LLE/VLE (T-x-x or T-x-y)",
        label_suffix="",
        required=False,
        initial=False,
        widget=forms.CheckboxInput(
            attrs={
                "class": "form-check-input",
                "aria-label": "Binary LLE/VLE (T-x-x or T-x-y)",
            }
        ),
    )


class BinaryVLECheckForm(forms.Form):
    "Form to check vapor-liquid equilibrium."

    binary_vle_checkbox = forms.BooleanField(
        label="Binary VLE (T-x-y + x-y)",
        label_suffix="",
        required=False,
        initial=False,
        widget=forms.CheckboxInput(
            attrs={
                "class": "form-check-input",
                "aria-label": "Binary VLE (T-x-y + x-y)",
            }
        ),
    )


class PhaseDiagramCheckForm(forms.Form):
    "Form to check phase diagram."

    phase_diagram_checkbox = forms.BooleanField(
        label="Phase Diag. (T-rho + P-rho)",
        label_suffix="",
        required=False,
        initial=False,
        widget=forms.CheckboxInput(
            attrs={
                "class": "form-check-input",
                "aria-label": "Phase Diag. (T-rho + P-rho)",
            }
        ),
    )


class STCheckForm(forms.Form):
    "Form to check surface tension."

    st_checkbox = forms.BooleanField(
        label="Surface tension (sigma-T)",
        label_suffix="",
        required=False,
        initial=False,
        widget=forms.CheckboxInput(
            attrs={
                "class": "form-check-input",
                "aria-label": "Surface tension (sigma-T)",
            }
        ),
    )
