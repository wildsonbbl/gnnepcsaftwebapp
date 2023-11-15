from django.db import models


class PropImage(models.Model):
    "Images for property plots."
    den = models.FileField(default="images/den.png")
    vp = models.FileField(default="images/vp.png")
