"Django forms."

import re
from typing import NamedTuple

from django import forms
from django.conf import settings
from django.core.exceptions import ValidationError
from django.utils.translation import gettext_lazy as _
from gnnepcsaft.data.ogb_utils import smiles2graph
from gnnepcsaft.data.rdkit_util import assoc_number, inchitosmiles, smilestoinchi


def _convert_query_to_identifiers(query: str, error_message: str) -> tuple[str, str]:
    """Convert a SMILES/InChI query into canonical SMILES and InChI."""

    try:
        if re.search("^InChI=", query):
            smiles = inchitosmiles(query, False, False)
            inchi = smilestoinchi(smiles, False, False)
        else:
            inchi = smilestoinchi(query, False, False)
            smiles = inchitosmiles(inchi, False, False)
        smiles2graph(smiles)
        assoc_number(inchi)
    except ValueError as e:
        raise ValidationError(_(error_message)) from e
    return smiles, inchi


def _parse_kij_values(kij_line: str, component_count: int) -> list[float]:
    """Parse the last line of the mixture textarea as kij values."""

    kij_values = kij_line.strip().split()
    expected_kij_count = (component_count**2 - component_count) / 2
    if len(kij_values) != expected_kij_count:
        raise ValidationError(
            _(
                f"Number of Kij values ({len(kij_values)}) must "
                f"be equal to {int(expected_kij_count)}."
            )
        )

    try:
        kij_values = [float(value) for value in kij_values]
    except ValueError as e:
        raise ValidationError(_("All Kij values must be valid numbers.")) from e

    for kij_value in kij_values:
        if kij_value < -1.0 or kij_value > 1.0:
            raise ValidationError(_("Kij values must be between -1 and 1."))

    return kij_values


def _parse_mole_fractions(raw_input: str, component_count: int) -> list[float]:

    mole_fractions = raw_input.strip().split()
    mole_fractions_list = []

    if component_count != len(mole_fractions):
        raise ValidationError(
            _(
                f"Number of mole fractions ({len(mole_fractions)}) must"
                f" match number of components ({component_count})"
            )
        )

    for mole_fraction in mole_fractions:
        try:
            mole_fractions_list.append(float(mole_fraction))
        except ValueError as e:
            raise ValidationError(_(f'Invalid Mole Fraction: "{mole_fraction}"')) from e

    total = sum(mole_fractions_list)
    if abs(total - 1.0) > 0.01:
        raise ValidationError(_(f"Mole fractions must sum to ~1.0, got {total:.3f}"))

    return mole_fractions_list


class InChIorSMILESinput(forms.Form):
    "Form to receive InChI/SMILES from user."

    query = forms.CharField(
        label="",
        strip=True,
        required=True,
        widget=forms.TextInput(
            attrs={
                "class": "form-control m-1",
                "aria-label": "Type/Paste InChI or SMILES",
                "placeholder": "Type/Paste InChI or SMILES",
            }
        ),
    )

    def clean_query(self):
        "check valid input and output SMILES."
        data = self.cleaned_data["query"]
        return _convert_query_to_identifiers(data, "Invalid InChI/SMILES.")


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
            smiles, inchi = _convert_query_to_identifiers(
                line.strip(), f"Invalid InChI/SMILES: {line.strip()}"
            )
            smiles_list.append(smiles)
            inchi_list.append(inchi)
        return inchi_list, smiles_list


class InChIorSMILESareaInputforMixture(forms.Form):
    "Form to receive InChI/SMILES + mole fractions from user."

    smiles_inchi_list = forms.CharField(
        label="",
        strip=True,
        required=True,
        widget=forms.TextInput(
            attrs={
                "class": "form-control my-2",
                "aria-label": "Type/Paste a list of InChI/SMILES",
                "placeholder": "SMILES/InChI separated by empty space (e.g., 'O=C=O O ...')",
            }
        ),
    )

    mole_fractions = forms.CharField(
        label="",
        strip=True,
        required=False,
        widget=forms.TextInput(
            attrs={
                "class": "form-control my-2",
                "aria-label": "Type/Paste a list of mole fractions",
                "placeholder": 'Mole Fractions separated by empty space (e.g., "0.3 0.2 ...")',
            }
        ),
    )

    kijs = forms.CharField(
        label="",
        strip=True,
        required=False,
        widget=forms.TextInput(
            attrs={
                "class": "form-control my-2",
                "aria-label": "Type/Paste a list of kij values",
                "placeholder": "kij values separated by empty space (e.g., k12 k13 k23 ...)",
            }
        ),
    )

    def clean_smiles_inchi_list(self):
        "check valid input and output SMILES."
        raw_input: str = self.data["smiles_inchi_list"]

        smiles_inchi_list = raw_input.split(" ")
        if len(smiles_inchi_list) < 2:
            raise ValidationError(_("Provide at least two 'SMILES/InChI'"))

        if settings.PLATFORM == "webapp" and len(smiles_inchi_list) > 10:
            raise ValidationError(_("Maximum 10 components allowed for web app."))
        inchi_list, smiles_list = [], []

        for smiles_or_inchi in smiles_inchi_list:
            smiles, inchi = _convert_query_to_identifiers(
                smiles_or_inchi,
                f'Invalid InChI/SMILES: "{smiles_or_inchi}"',
            )
            smiles_list.append(smiles)
            inchi_list.append(inchi)
        return inchi_list, smiles_list

    def clean_mole_fractions(self):
        "clean mole fractions"
        _inchi_list, smiles_list = self.clean_smiles_inchi_list()

        raw_input: str = self.cleaned_data["mole_fractions"]

        if raw_input:
            return _parse_mole_fractions(raw_input, len(smiles_list))
        return [round(1.0 / len(smiles_list), 4) for _ in smiles_list]

    def clean_kijs(self):
        "clean kij values"
        _inchi_list, smiles_list = self.clean_smiles_inchi_list()
        raw_input: str = self.cleaned_data["kijs"]
        if raw_input:
            return _parse_kij_values(raw_input, len(smiles_list))
        return []


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


class TernaryVLEpxCheckForm(forms.Form):
    "Form to check VLE equilibrium."

    ternary_vlepx_checkbox = forms.BooleanField(
        label="Ternary VLE P-x1 (T and x2/x3 fixed)",
        label_suffix="",
        required=False,
        initial=False,
        help_text="<p class='form-text'>- More volatile component must "
        "be listed first in P-x1 calculation<br>"
        "- Solvent ratio must be between 0 and 1</p>",
        widget=forms.CheckboxInput(
            attrs={
                "class": "form-check-input",
                "aria-label": "Ternary VLE P-x1 (T and x2/x3 fixed)",
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
        help_text="<p class='form-text'>- Mole fractions and minimum temperature are used "
        "as starting value for calculations</p>",
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


class BinaryVLEpxyCheckForm(forms.Form):
    "Form to check vapor-liquid equilibrium."

    binary_vlepxy_checkbox = forms.BooleanField(
        label="Binary VLE (P-x-y)",
        label_suffix="",
        required=False,
        initial=False,
        help_text="<p class='form-text'>- More volatile component must "
        "be listed first in P-x-y calculation </p>",
        widget=forms.CheckboxInput(
            attrs={
                "class": "form-check-input",
                "aria-label": "Binary VLE (P-x-y)",
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


class KijCheckForm(forms.Form):
    "Form to check estimate kij"

    kij_checkbox = forms.BooleanField(
        label="Estimate kij value for binary mixture",
        label_suffix="",
        required=False,
        initial=False,
        help_text="<p class='form-text'>- All available VLE exp. data used in optimization</p>",
        widget=forms.CheckboxInput(
            attrs={
                "class": "form-check-input",
                "aria-label": "Estimate kij value for binary mixture",
            }
        ),
    )


class PureForms(NamedTuple):
    """Container for pure-page forms."""

    form: InChIorSMILESinput
    plot_config: CustomPlotConfigForm
    rho_checkbox: RhoCheckForm
    vp_checkbox: VPCheckForm
    h_lv_checkbox: HlvCheckForm
    s_lv_checkbox: SlvCheckForm
    phase_diagram_checkbox: PhaseDiagramCheckForm
    st_checkbox: STCheckForm


class MixtureForms(NamedTuple):
    """Container for mixture-page forms."""

    form: InChIorSMILESareaInputforMixture
    plot_config: CustomPlotConfigForm
    rho_checkform: RhoCheckForm
    vp_checkform: VPCheckForm
    binary_vle_checkform: BinaryVLECheckForm
    binary_lle_checkform: BinaryLLECheckForm
    ternary_lle_checkform: TernaryLLECheckForm
    binary_vlepxy_checkform: BinaryVLEpxyCheckForm
    ternary_vlepx_checkform: TernaryVLEpxCheckForm
    kij_checkform: KijCheckForm
