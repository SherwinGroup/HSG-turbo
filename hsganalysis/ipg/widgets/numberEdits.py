from PyQt5 import QtWidgets, QtGui, QtCore

class QFNumberEdit(QtGui.QLineEdit):
    #a signal to emit the new, approved number. Will emit False if the
    # inputted value is not accepted. Intended for float inputs
    textAccepted = QtCore.pyqtSignal(object)
    def __init__(self, parent = None, contents = ''):
        super(QFNumberEdit, self).__init__(parent)
        self.editingFinished.connect(self._handleEditingFinished)
        self.textChanged.connect(lambda: self._handleEditingFinished())
        self.returnPressed.connect(lambda: self._handleEditingFinished(True))
        self._before = contents



    def _handleEditingFinished(self, _return = False):
        before, after = self._before, str(self.text())
        if (not self.hasFocus() or _return) and before != after:
            val = self.parseInp(after)
            #if the return is False, need to catch that. Otherwise, may take
            #if val to be false when val=0, which is a valid input
            if type(val) is bool:
                self.setText(str(before))
                self.textAccepted.emit(False)
            else:
                self.setText(str(val))
                self._before = str(val)
                self.textAccepted.emit(val)


    def value(self):
        ret = -1
        if str(self.text()) == '':
            return float(self._before)
        try:
            ret = float(self.text())
        except:
            self._handleEditingFinished()
            ret = float(self.text())
        return ret

    def parseInp(self, inp):
        ret = None
        #see if we can just turn it into a number and leave if we can
        try:
            ret = float(inp)
            return ret
        except:
            pass

        toMatch = re.compile('(\d+\.?\d*|\d*\.\d+)\*(\d+\.?\d*|\d*\.\d+)')
        if re.match(toMatch, inp):
            print(("it's a command! {}".format(inp)))
            try:
                ret = eval(inp)
                return ret
            except Exception as e:
                print("Can't parse command", inp, e)
        #tests to see whether digit is whole number or decimal, and if it has
        #some modifier at the end
        toMatch = re.compile('-?(\d+\.?\d*|\d*\.\d+)(m|u|n|M|k)?\Z')
        if re.match(toMatch, inp):
            convDict = {'m': 1e-3, 'u':1e-6, 'n':1e-9, 'M':1e6, 'k':1e3}
            try:
                ret = (float(inp[:-1]) * #convert first part to number
                   convDict[[b for b in list(convDict.keys()) if b in inp][0]]) #and multiply by the exponential
                return ret
            except Exception as e:
                print("Can't parse float string")
                print(inp, type(inp))
                print(e)
                print('')
        else:
            return False

class QINumberEdit(QtGui.QLineEdit):
    #a signal to emit the new, approved number. Will emit False if the
    # inputted value is not accepted. Intended for integer inputs
    textAccepted = QtCore.pyqtSignal(object)
    def __init__(self, parent = None, contents = ''):
        super(QINumberEdit, self).__init__(parent)
        self.editingFinished.connect(self._handleEditingFinished)
        self.textChanged.connect(lambda: self._handleEditingFinished())
        self.returnPressed.connect(lambda: self._handleEditingFinished(True))
        self._before = contents

    def __add__(self, other):
        ret = copy.deepcopy(self)
        ret.setText(str(self.value()+other))
        return ret

    def __sub__(self, other):
        ret = copy.deepcopy(self)
        ret.setText(str(self.value()-other))
        return ret

    def _handleEditingFinished(self, _return = False):
        before, after = self._before, str(self.text())
        if (not self.hasFocus() or _return) and before != after:
            val = self.parseInp(after)
            #if the return is False, need to catch that. Otherwise, may take
            #if val to be false when val=0, which is a valid input
            if type(val) is bool:
                self.setText(str(before))
                self.textAccepted.emit(False)
            else:
                self.setText(str(val))
                self._before = str(val)
                self.textAccepted.emit(val)
    def value(self):
        ret = -1
        if str(self.text())=='':
            return int(self._before)
        try:
            ret = int(self.text())
        except:
            self._handleEditingFinished()
            ret = int(self.text())
        return ret


    def parseInp(self, inp):
        ret = None
        #see if we can just turn it into a number and leave if we can
        try:
            ret = int(inp)
            return ret
        except:
            return False
