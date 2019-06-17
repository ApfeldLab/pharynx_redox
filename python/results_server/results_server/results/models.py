from django.db import models

# Create your models here.

class Experiment(models.Model):
    date = models.DateField('experiment date')
    data = models.BinaryField()
    experimenter = models.CharField(max_length=100)
