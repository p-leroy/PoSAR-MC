import sys
import random
from PySide2.QtWidgets import (QApplication, QLabel, QPushButton, QVBoxLayout, QWidget)
from PySide2.QtCore import Slot, Qt

class MyWidget(QWidget):
    def __init__(self):
        QWidget.__init__(self)

        self.hello = ["Hallo Welt", "你好，世界", "Hei maailma", "Hola Mundo", "Привет мир"]

        self.button = QPushbutton("Click me!")
        self.text = QLabel("Hello World!")
        self.text.setAlignment(Qt.AlignCenter)

if __name__ == "__main__":
    app = QApplication(sys.argv)

    widget = MyWidget()
    widget.resize(800, 600)
    widget.show()

    sys.exit(app.exec())
