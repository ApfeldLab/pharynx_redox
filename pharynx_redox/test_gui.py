import sys
import time

from PyQt5.QtWidgets import QApplication
from PyQt5.QtCore import QTime

app = QApplication(sys.argv)


try:
    due = QTime.currentTime()
    message = "Alert!"
    if len(sys.argv) < 2:
        raise ValueError
    hours, mins = sys.argv[1].split(":")
    due = QTime(int(hours), int(mins))
    if not due.isValid():
        raise ValueError
    if len(sys.argv) > 2:
        message = " ".join(sys.argv[2:])
except ValueError:
    message = "Usage: alert.pyw HH:MM [optional message]"  # 24H clock

while QTime.currentTime() < due:
    time.sleep(20)

label = QLabel("<font color=red size=72><b>" + message + "</b><font>")
label.setWindowFlags(Qt.)