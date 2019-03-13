import numpy as np
import pyqtgraph
from PyQt5 import QtCore, QtGui

dialogList = []
class MessageDialog(QtGui.QDialog):
    def __init__(self, parent, message="", duration=3000):
        if isinstance(parent, str):
            message = parent
            parent = None
        super(MessageDialog, self).__init__(parent=parent)
        layout  = QtGui.QVBoxLayout(self)
        # text = QtGui.QLabel("<font size='6'>{}</font>".format(message), self)
        text = QtGui.QLabel(self)
        text.setTextFormat(QtCore.Qt.RichText)
        text.setText("<font size='6'>{}</font>".format(message))
        text.setWordWrap(True)
        layout.addWidget(text)
        self.setLayout(layout)
        self.setModal(False)

        dialogList.append(self)

        if duration:
            self.timer = QtCore.QTimer.singleShot(duration, self.close)
        self.show()
        self.raise_()

    def close(self):
        try:
            dialogList.remove(self)
        except Exception as E:
            print(("Error removing from list, ",E))

        super(MessageDialog, self).close()



















































