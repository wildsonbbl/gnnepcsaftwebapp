from django.urls import path

from . import views

urlpatterns = [
    path("estimator/", views.index, name="estimator"),
    path("", views.homepage, name="homepage"),
    path("author/", views.authorpage, name="author"),
]
