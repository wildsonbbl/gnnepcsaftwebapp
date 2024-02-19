"Module for url config."
from django.urls import path

from . import views

urlpatterns = [
    path("estimator/", views.estimator, name="estimator"),
    path("", views.homepage, name="homepage"),
    path("author/", views.authorpage, name="author"),
    path("estimator/description", views.description, name="description"),
]
