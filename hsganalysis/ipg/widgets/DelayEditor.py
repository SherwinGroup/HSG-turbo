__author__ = 'Home'
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtCore import Qt


# Map the position of the text string to the exponential
# it represents
posMapper = {
     0:   2,
     1:   1,
     2:   0,

     4:  -1,
     5:  -2,
     6:  -3,

     8:  -4,
     9:  -5,
    10:  -6,

    12:  -7,
    13:  -8,
    14:  -9,

    16: -10,
    17: -11,
    18: -12
}

class DelayEditor(QtWidgets.QLineEdit):

    # emits signal when the user has updated the vaalue in the
    # text box. Emits the new values
    sigNewValue = QtCore.pyqtSignal(object)



    MAX_VALUE = 999.999999999995
    def __init__(self, *args, **kwargs):
        super(DelayEditor, self).__init__(*args, **kwargs)

        val = kwargs.pop("value", 0)
        self.setValue(val)
        self.setHighlight(0)
        self.emitOnReturnOnly = False

        self.setFixedWidth(110)
        self.setFixedHeight(17)

        self.setContextMenuPolicy(Qt.NoContextMenu)

    def keyPressEvent(self, QKeyEvent):
        """
        :type QKeyEvent: QtGui.QKeyEvent
        """
        cur = self.value()
        if QKeyEvent.key() in [Qt.Key_Left, Qt.Key_Right,
                               Qt.Key_Up, Qt.Key_Down]:
            QKeyEvent.accept()
            if QKeyEvent.key() == Qt.Key_Left:
                self.setHighlight(dx=-1)
            elif QKeyEvent.key() == Qt.Key_Right:
                self.setHighlight(dx=1)
            elif QKeyEvent.key() == Qt.Key_Up:
                self.updateValue(1)
                if not self.emitOnReturnOnly:
                    self.sigNewValue.emit(self.value())
            elif QKeyEvent.key() == Qt.Key_Down:
                self.updateValue(-1)
                if not self.emitOnReturnOnly:
                    self.sigNewValue.emit(self.value())
        elif QKeyEvent.key() == Qt.Key_Home:
            self.setHighlight(0)
            return
        elif QKeyEvent.key() == Qt.Key_End:
            self.setHighlight(len(self.text()))
            return
        elif QKeyEvent.key() in [Qt.Key_0, Qt.Key_1, Qt.Key_2,
                                 Qt.Key_3, Qt.Key_4, Qt.Key_5,
                                 Qt.Key_6, Qt.Key_7, Qt.Key_8,
                                           Qt.Key_9          ]:
            QKeyEvent.accept()
            pos = self.selectionStart()
            super(DelayEditor, self).keyPressEvent(QKeyEvent)
            self.setHighlight(pos+1)
            if not self.emitOnReturnOnly:
                self.sigNewValue.emit(self.value())
        elif QKeyEvent.key() in [Qt.Key_Enter, Qt.Key_Return]:
            self.sigNewValue.emit(self.value())

        # super(DelayEditor, self).keyPressEvent(QKeyEvent)

    def mousePressEvent(self, QMouseEvent):
        """
        :type QMouseEvent: QtGui.QMouseEvent
        :param QMouseEvent:
        :return:
        """
        if QMouseEvent.button() == Qt.RightButton:
            return
        super(DelayEditor, self).mousePressEvent(QMouseEvent)
        self.setHighlight(self.cursorPosition())

    def mouseDoubleClickEvent(self, QMouseEvent):
        QMouseEvent.ignore()

    def mouseReleaseEvent(self, QMouseEvent):
        QMouseEvent.ignore()

    def mouseMoveEvent(self, QMouseEvent):
        QMouseEvent.ignore()

    def setHighlight(self, pos=None, dx=None):
        if dx is not None:
            pos = self.selectionStart()+dx

        if pos>=len(self.text()):
            pos = len(self.text())-1
        if pos < 0:
            pos = 0
        self.setSelection(pos, 1)
        if str(self.text())[pos] in ['.', ' ']:
            if dx is None: dx=1
            self.setHighlight(dx=dx)

    def updateValue(self, dy=0):
        # pos = self.selectionStart()
        # curString = list(str(self.text()))
        # curVal = int(curString[pos])
        # newVal = (curVal + dy + 10) % 10
        # curString[pos] = str(newVal)
        # self.setText(''.join(curString))
        # self.setHighlight(pos)
        pos = self.selectionStart()
        curVal = self.value()
        exp = posMapper[pos]
        self.setValue(curVal + dy * 10**exp)
        # print "want to update, cur pos", pos,
        # print "suggest", posMapper[pos]

        self.setHighlight(pos)

    def setValue(self, val):
        val = float(val)
        if val>DelayEditor.MAX_VALUE:
            val = DelayEditor.MAX_VALUE
        if val<0:
            val = 0
        st = "{:016.12f}".format(val)

        point = st.find('.')

        newSt = st[:point+1]
        subSt = st[point+1:]
        newSt += ' '.join([subSt[ii:ii+3] for ii in range(0, len(subSt), 3)])

        self.setText(newSt)

    def value(self):
        return float(str(self.text()).replace(' ',''))

class DelayTimeEditor(QtWidgets.QWidget):
    def __init__(self, *args, **kwargs):
        super(DelayTimeEditor, self).__init__(*args, **kwargs)
        self.layout = QtWidgets.QGridLayout()
        self.editor = DelayEditor(*args, **kwargs)
        self.layout.addWidget(self.editor, 0, 0, 1, 5)

        pref = ['', 'm', 'u', 'n', 'p']

        # self.layout.addItem(QtGui.QSpacerItem(10, 0), 1, 0)
        for ii in range(1, 5):
            label = QtWidgets.QLabel(self)
            label.setText("{}s".format(pref[ii]))
            label.setMargin(0)
            self.layout.addWidget(label, 1, ii, Qt.AlignCenter)

        self.layout.setSpacing(0)

        self.setFixedWidth(130)
        self.setFixedHeight(53)
        self.setLayout(self.layout)

if __name__ == '__main__':
    app = QtWidgets.QApplication([])
    # wid = DelayEditor()
    wid = DelayTimeEditor()
    wid.show()
    print("made")

    app.exec_()
