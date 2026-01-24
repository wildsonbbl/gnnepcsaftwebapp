"Module for url config."

from django.urls import path, re_path

from . import views

urlpatterns = [
    path("", views.pure, name="pure"),
    path("about/", views.about, name="about"),
    path("batch/", views.batch, name="batch"),
    path("mixture/", views.mixture, name="mixture"),
    re_path(r"^serviceworker\.js$", views.service_worker, name="serviceworker"),
]
