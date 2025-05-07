"Module for url config."

from django.urls import path, re_path

from . import views

urlpatterns = [
    path("", views.pure, name="pure"),
    path("about/", views.about, name="about"),
    path("batch/", views.batch, name="batch"),
    path("mixture/", views.mixture, name="mixture"),
    # path("description/", views.description, name="description"),
    path("chat/", views.chat, name="chat"),
    # path("", views.homepage, name="homepage"),
    # path("author/", views.authorpage, name="author"),
    path("api/sessions/", views.get_sessions, name="get_sessions"),
    path("api/sessions/create/", views.create_session, name="create_session"),
    path(
        "api/sessions/<uuid:session_id>/delete/",
        views.delete_session,
        name="delete_session",
    ),
    re_path(r"^serviceworker\.js$", views.service_worker, name="serviceworker"),
    path("offline/", views.offline, name="offline"),
]
