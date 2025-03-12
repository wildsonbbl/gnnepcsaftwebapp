"Module for url config."

from django.urls import path

from . import views

urlpatterns = [
    path("", views.estimator, name="estimator"),
    path("about/", views.about, name="about"),
    # path("", views.homepage, name="homepage"),
    # path("author/", views.authorpage, name="author"),
]
