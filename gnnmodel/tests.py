"""Module for django tests."""

import os
from unittest.mock import patch

from django.http import QueryDict
from django.test import Client, TestCase
from django.urls import reverse

from .forms import (
    CustomPlotConfigForm,
    InChIorSMILESareaInput,
    InChIorSMILESareaInputforMixture,
    InChIorSMILESinput,
)


class ViewsTestCase(TestCase):
    """Test case for views."""

    def setUp(self):
        """Set up test environment."""
        self.client = Client()
        # Sample data for testing
        self.valid_inchi = "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
        self.valid_smiles = "CCO"  # Ethanol
        self.valid_parameters = [
            2.8733,
            2.9623,
            187.377,
            0.0559,
            2460.6204,
            0.0,
            1,
            1,
        ]
        # Store the original API key if it exists
        self.original_api_key = os.environ.get("GOOGLE_API_KEY")

    def tearDown(self):
        """Clean up after tests."""
        # Restore the original API key if it existed
        if self.original_api_key:
            os.environ["GOOGLE_API_KEY"] = self.original_api_key
        elif "GOOGLE_API_KEY" in os.environ:
            del os.environ["GOOGLE_API_KEY"]

    @patch("gnnmodel.views.build_pure_context")
    @patch("gnnmodel.views.process_pure_post")
    @patch("gnnmodel.views.init_pure_forms")
    def test_pure_get(
        self, mock_init_pure_forms, mock_process_pure_post, mock_build_pure_context
    ):
        """Test GET request to pure view."""
        mock_init_pure_forms.return_value = "forms"
        mock_build_pure_context.return_value = {
            "output": False,
            "form": InChIorSMILESinput(),
            "plot_config": CustomPlotConfigForm(),
        }
        response = self.client.get(reverse("pure"))
        self.assertEqual(response.status_code, 200)
        self.assertTemplateUsed(response, "pure.html")
        self.assertFalse(response.context["output"])
        self.assertIsInstance(response.context["form"], InChIorSMILESinput)
        self.assertIsInstance(response.context["plot_config"], CustomPlotConfigForm)
        mock_init_pure_forms.assert_called_once()
        mock_process_pure_post.assert_not_called()
        mock_build_pure_context.assert_called_once_with("forms")

    @patch("gnnmodel.views.build_pure_context")
    @patch("gnnmodel.views.process_pure_post")
    @patch("gnnmodel.views.init_pure_forms")
    def test_pure_post_valid(
        self, mock_init_pure_forms, mock_process_pure_post, mock_build_pure_context
    ):
        """Test POST request to pure view with valid data."""
        mock_init_pure_forms.return_value = "forms"
        mock_process_pure_post.return_value = {"output": True}
        mock_build_pure_context.return_value = {
            "output": True,
            "form": InChIorSMILESinput(),
            "plot_config": CustomPlotConfigForm(),
        }

        form_data = {
            "query": self.valid_inchi,
            "custom_plot_checkbox": "False",
            "rho_checkbox": "True",
            "vp_checkbox": "True",
            "h_lv_checkbox": "False",
            "s_lv_checkbox": "False",
            "phase_diagram_checkbox": "False",
            "st_checkbox": "False",
            "temp_min": "300.0",
            "temp_max": "400.0",
            "pressure": "101325.0",
        }
        form_querydict = QueryDict("", mutable=True)
        form_querydict.update(form_data)

        response = self.client.post(reverse("pure"), form_data)
        self.assertEqual(response.status_code, 200)
        self.assertTemplateUsed(response, "pure.html")
        self.assertTrue(response.context["output"])
        mock_init_pure_forms.assert_called_once_with(form_querydict)
        mock_process_pure_post.assert_called_once_with("forms")
        mock_build_pure_context.assert_called_once_with("forms", {"output": True})

    @patch("gnnmodel.views.build_mixture_context")
    @patch("gnnmodel.views.process_mixture_post")
    @patch("gnnmodel.views.init_mixture_forms")
    def test_mixture_get(
        self,
        mock_init_mixture_forms,
        mock_process_mixture_post,
        mock_build_mixture_context,
    ):
        """Test GET request to mixture view."""
        mock_init_mixture_forms.return_value = "forms"
        mock_build_mixture_context.return_value = {
            "output": False,
            "form": InChIorSMILESareaInputforMixture(),
        }
        response = self.client.get(reverse("mixture"))
        self.assertEqual(response.status_code, 200)
        self.assertTemplateUsed(response, "mixture.html")
        self.assertFalse(response.context["output"])
        self.assertIsInstance(
            response.context["form"], InChIorSMILESareaInputforMixture
        )
        mock_init_mixture_forms.assert_not_called()
        mock_process_mixture_post.assert_not_called()
        mock_build_mixture_context.assert_called_once_with()

    @patch("gnnmodel.views.build_mixture_context")
    @patch("gnnmodel.views.process_mixture_post")
    @patch("gnnmodel.views.init_mixture_forms")
    def test_mixture_post_valid(
        self,
        mock_init_mixture_forms,
        mock_process_mixture_post,
        mock_build_mixture_context,
    ):
        """Test POST request to mixture view with valid data."""
        mock_init_mixture_forms.return_value = "forms"
        mock_process_mixture_post.return_value = {"output": True}
        mock_build_mixture_context.return_value = {
            "output": True,
            "form": InChIorSMILESareaInputforMixture(),
        }

        form_data = {
            "text_area": f"{self.valid_inchi} 0.5\n{self.valid_inchi} 0.5",
            "temp_min": "298.15",
            "temp_max": "350.0",
            "pressure": "101325.0",
        }
        form_querydict = QueryDict("", mutable=True)
        form_querydict.update(form_data)

        response = self.client.post(reverse("mixture"), form_data)
        self.assertEqual(response.status_code, 200)
        self.assertTemplateUsed(response, "mixture.html")
        self.assertTrue(response.context["output"])
        mock_init_mixture_forms.assert_called_once_with(form_querydict)
        mock_process_mixture_post.assert_called_once_with("forms")
        mock_build_mixture_context.assert_called_once_with({"output": True})

    @patch("gnnmodel.views.get_pred")
    def test_batch_get(self, mock_get_pred):
        """Test GET request to batch view."""
        response = self.client.get(reverse("batch"))
        self.assertEqual(response.status_code, 200)
        self.assertTemplateUsed(response, "batch.html")
        self.assertFalse(response.context["output"])
        self.assertIsInstance(response.context["form"], InChIorSMILESareaInput)
        self.assertEqual(mock_get_pred.call_count, 0)

    @patch("gnnmodel.views.get_pred")
    def test_batch_post_valid(self, mock_get_pred):
        """Test POST request to batch view with valid data."""
        # Mock return values
        mock_get_pred.return_value = self.valid_parameters

        # Create form data with multiple compounds
        form_data = {
            "text_area": f"{self.valid_inchi}\n{self.valid_inchi}",
        }

        response = self.client.post(reverse("batch"), form_data)
        self.assertEqual(response.status_code, 200)
        self.assertTemplateUsed(response, "batch.html")
        self.assertEqual(mock_get_pred.call_count, 2)  # Called twice for two compounds
        mock_get_pred.assert_called_with(self.valid_smiles)

    def test_about(self):
        """Test about page view."""
        response = self.client.get(reverse("about"))
        self.assertEqual(response.status_code, 200)
        self.assertTemplateUsed(response, "about.html")
