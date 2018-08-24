import numpy as np
from scipy.optimize import minimize

from joblib import Parallel, delayed, cpu_count
from hsganalysis.jones import JonesVector as JV
from .expFanCompiler import FanCompiler
from .extractMatrices import *

def generateMC(alpha, gamma=None, angle=0, niter = 4000,
               ncores = cpu_count(),
               returnFull=False):
    """
    Provided alpha/gamma data, calculate T matrices with error by implementing a
    Monte Carlo error propagation. Takes the input alpha/gamma angles, add Gaussian
    noise, and calculate the T matrix. Repeat this process many times, and return the
    matrix of all calculated Ts, for potential statistics purposes. Parallelizes over
    multiple processors.
    :param alpha: Alpha angles as Nx3 [sb, alpha, alpha err]
    :param gamma: Gamma angles as Nx3 [sb, gamma, gamma err]
    :param angle: The angle of the THz WITH RESPECT TO [011]
    :param niter: How many iterations/repititions of the MC to run
    :param ncores: How many cores to run it on. Defaults to all
    :param returnFull: If true, returns an Nx2x2xniter of all the T matrices. If
        false, returns (mean, std) of the dataset
    :return:
    """
    if gamma is None and isinstance(alpha, FanCompiler):
        alpha, gamma, _ = alpha.build()
    alphaErr = np.column_stack((np.zeros_like(alpha[:,0]), alpha[:, 2::2]))
    alpha = np.column_stack((alpha[:, 0], alpha[:, 1::2]))

    gammaErr = np.column_stack((np.zeros_like(gamma[:, 0]), gamma[:, 2::2]))
    gamma = np.column_stack((gamma[:, 0], gamma[:, 1::2]))

    def doSingleMC():
        global failureRate
        thisAlpha = alpha + np.random.normal(scale=alphaErr)
        thisGamma = gamma + np.random.normal(scale=gammaErr)
        try:
            J = findJ(thisAlpha, thisGamma)
            T = makeT(J, angle)
            return T
        except Exception as e:
            print("Exception occurred")
            print(e)
            failureRate += 1


    allTs = Parallel(n_jobs=ncores, verbose=5)(delayed(doSingleMC)() for _ in range(
        niter))

    # Transpose it so the various iterations are on the final axis,
    # meanT=allTs.mean(axis=-1)
    allTs = np.array(allTs).transpose(1, 2, 3, 0)

    if returnFull:
        return allTs
    else:
        meanT = allTs.mean(axis=-1)

        # This is maybe a difficult point. The standard deviation is defined as sqrt(A*A)
        # (A* being the conjugate transpose), which means the std is a real number. Normally, this
        # makes sense, but I don't know if that's what's wanted here. The real and imaginary
        # parts are independently calculated in the fitting routine, which means each has its own
        #  error associated with it., so I'd rather keep those around like this.
        # stdT = allTs.std(axis=-1)
        stdT = np.real(allTs).std(axis=-1) + 1j*np.imag(allTs).std(axis=-1)
        return meanT, stdT


def errorPropagator(num, den, dnum, dden):
    """
    Do error propagation for two calculations on two complex numbers. Calculates
    |num/den| and angle(num/den), and provides errors on them. Error propagation
    determiend by calculating |num/den| = sqrt((a+ib/x+iy) * (a-ib/x-iy)), and doing
    error propagation on that. Similar for the angle, where angle=arctan(num/den).
    Propagations were done in Mathematica because that's too much algebra for me to
    mess up.
    :param num:
    :param den:
    :param dnum:
    :param dden:
    :return:
    """
    a, b = np.real(num), np.imag(num)
    x, y = np.real(den), np.imag(den)

    da, db = np.real(dnum), np.imag(dnum)
    dx, dy = np.real(dden), np.imag(dden)

    ratio = np.abs(num/den)
    dratio = np.sqrt(((a**2+b**2)**2*dx**2*x**2+(a**2+b**2)**2*dy**2*y**2+a**2*da**2*(x**2+y**2)**2+b**2*db**2*(x**2+y**2)**2)/((a**2+b**2)*(x**2+y**2)**3))

    ang = np.angle(num/den, deg=True)
    dang = np.sqrt(a**4*(dy**2*x**2+dx**2*y**2)+b**2*(da**2*(x**2+y**2)**2+b**2*(
            dy**2*x**2+dx**2*y**2))+a**2*(db**2*(x**2+y**2)**2+2*b**2*(dy**2*x**2+dx**2*y**2)))/(
            (a**2+b**2)*(x**2+y**2))
    dang = np.rad2deg(np.abs(dang))


    return ratio, dratio, ang, dang


def errorPropagatorDifference(num, den, dnum, dden):
    """
    Do error propagation for two calculations on two complex numbers. Calculates
    |num-den| and angle(num-den), and provides errors on them. Error propagation
    determiend by calculating |num-den| = sqrt((a+ib-(x+iy)) * (a-ib-(x-iy))),
    and doing error propagation on that. Similar for the angle, where angle=arctan(
    num/den). Propagations were done in Mathematica because that's too much algebra
    for me to mess up.
    :param num:
    :param den:
    :param dnum:
    :param dden:
    :return:
    """
    a, b = np.real(num), np.imag(num)
    x, y = np.real(den), np.imag(den)

    da, db = np.real(dnum), np.imag(dnum)
    dx, dy = np.real(dden), np.imag(dden)

    mag = np.abs(num-den)
    dmag = np.sqrt(((da**2+dx**2)*(a-x)**2+(db**2+dy**2)*(b-y)**2)/((a-x)**2+(b-y)**2))

    ang = -np.angle(num-den, deg=True)
    dang = np.sqrt((db**2+dy**2)*(a-x)**2+(da**2+dx**2)*(b-y)**2)/((a-x)**2+(b-y)**2)
    dang = np.rad2deg(np.abs(dang))


    return mag, dmag, ang, dang



try:
    from PyQt5 import QtWidgets, QtGui, QtCore
    class Ui_Form(object):
        def setupUi(self, Form):
            Form.setObjectName("Form")
            Form.resize(582, 555)
            self.verticalLayout_3 = QtWidgets.QVBoxLayout(Form)
            self.verticalLayout_3.setObjectName("verticalLayout_3")
            self.horizontalLayout = QtWidgets.QHBoxLayout()
            self.horizontalLayout.setObjectName("horizontalLayout")
            self.groupBox_2 = QtWidgets.QGroupBox(Form)
            self.groupBox_2.setObjectName("groupBox_2")
            self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.groupBox_2)
            self.verticalLayout_2.setObjectName("verticalLayout_2")
            self.rbxreal = QtWidgets.QRadioButton(self.groupBox_2)
            self.rbxreal.setObjectName("rbxreal")
            self.verticalLayout_2.addWidget(self.rbxreal)
            self.rbximag = QtWidgets.QRadioButton(self.groupBox_2)
            self.rbximag.setObjectName("rbximag")
            self.verticalLayout_2.addWidget(self.rbximag)
            self.lwx = QtWidgets.QListWidget(self.groupBox_2)
            self.lwx.setObjectName("lwx")
            item = QtWidgets.QListWidgetItem()
            self.lwx.addItem(item)
            item = QtWidgets.QListWidgetItem()
            self.lwx.addItem(item)
            item = QtWidgets.QListWidgetItem()
            self.lwx.addItem(item)
            item = QtWidgets.QListWidgetItem()
            self.lwx.addItem(item)
            self.verticalLayout_2.addWidget(self.lwx)
            self.horizontalLayout.addWidget(self.groupBox_2)
            self.groupBox = QtWidgets.QGroupBox(Form)
            self.groupBox.setObjectName("groupBox")
            self.verticalLayout = QtWidgets.QVBoxLayout(self.groupBox)
            self.verticalLayout.setObjectName("verticalLayout")
            self.rbyreal = QtWidgets.QRadioButton(self.groupBox)
            self.rbyreal.setObjectName("rbyreal")
            self.buttonGroup = QtWidgets.QButtonGroup(Form)
            self.buttonGroup.setObjectName("buttonGroup")
            self.buttonGroup.addButton(self.rbyreal)
            self.verticalLayout.addWidget(self.rbyreal)
            self.rbyimag = QtWidgets.QRadioButton(self.groupBox)
            self.rbyimag.setObjectName("rbyimag")
            self.buttonGroup.addButton(self.rbyimag)
            self.verticalLayout.addWidget(self.rbyimag)
            self.lwy = QtWidgets.QListWidget(self.groupBox)
            self.lwy.setObjectName("lwy")
            item = QtWidgets.QListWidgetItem()
            self.lwy.addItem(item)
            item = QtWidgets.QListWidgetItem()
            self.lwy.addItem(item)
            item = QtWidgets.QListWidgetItem()
            self.lwy.addItem(item)
            item = QtWidgets.QListWidgetItem()
            self.lwy.addItem(item)
            self.verticalLayout.addWidget(self.lwy)
            self.horizontalLayout.addWidget(self.groupBox)
            self.lwSBs = QtWidgets.QListWidget(Form)
            self.lwSBs.setObjectName("lwSBs")
            self.horizontalLayout.addWidget(self.lwSBs)
            spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
            self.horizontalLayout.addItem(spacerItem)
            self.horizontalLayout.setStretch(0, 1)
            self.horizontalLayout.setStretch(1, 1)
            self.horizontalLayout.setStretch(2, 1)
            self.horizontalLayout.setStretch(3, 10)
            self.verticalLayout_3.addLayout(self.horizontalLayout)
            self.gPlot = PlotWidget(Form)
            self.gPlot.setObjectName("gPlot")
            self.verticalLayout_3.addWidget(self.gPlot)
            self.verticalLayout_3.setStretch(0, 1)
            self.verticalLayout_3.setStretch(1, 100)

            self.retranslateUi(Form)
            QtCore.QMetaObject.connectSlotsByName(Form)

        def retranslateUi(self, Form):
            _translate = QtCore.QCoreApplication.translate
            Form.setWindowTitle(_translate("Form", "Form"))
            self.groupBox_2.setTitle(_translate("Form", "x"))
            self.rbxreal.setText(_translate("Form", "Real"))
            self.rbximag.setText(_translate("Form", "Imag"))
            __sortingEnabled = self.lwx.isSortingEnabled()
            self.lwx.setSortingEnabled(False)
            item = self.lwx.item(0)
            item.setText(_translate("Form", "T++"))
            item = self.lwx.item(1)
            item.setText(_translate("Form", "T+-"))
            item = self.lwx.item(2)
            item.setText(_translate("Form", "T-+"))
            item = self.lwx.item(3)
            item.setText(_translate("Form", "T--"))
            self.lwx.setSortingEnabled(__sortingEnabled)
            self.groupBox.setTitle(_translate("Form", "y"))
            self.rbyreal.setText(_translate("Form", "Real"))
            self.rbyimag.setText(_translate("Form", "Imag"))
            __sortingEnabled = self.lwy.isSortingEnabled()
            self.lwy.setSortingEnabled(False)
            item = self.lwy.item(0)
            item.setText(_translate("Form", "T++"))
            item = self.lwy.item(1)
            item.setText(_translate("Form", "T+-"))
            item = self.lwy.item(2)
            item.setText(_translate("Form", "T-+"))
            item = self.lwy.item(3)
            item.setText(_translate("Form", "T--"))
            self.lwy.setSortingEnabled(__sortingEnabled)

    class MonteCarloInteracter(QtWidgets.QWidget):
        def __init__(self, data):
            super(Window, self).__init__()
            self.data = data
            self.ui = Ui_Form()
            self.ui.setupUi(self)

            self.ui.lwSBs.addItems(["{}".format(ii) for ii in wantedSBs])

            self.ui.lwSBs.itemSelectionChanged.connect(self.updatePlot)
            self.ui.lwSBs.setSelectionMode(self.ui.lwSBs.ExtendedSelection)
            self.ui.lwx.itemSelectionChanged.connect(self.updatePlot)
            self.ui.lwy.itemSelectionChanged.connect(self.updatePlot)
            self.ui.rbxreal.clicked.connect(self.updatePlot)
            self.ui.rbximag.clicked.connect(self.updatePlot)
            self.ui.buttonGroup.buttonClicked.connect(self.updatePlot)


        def updatePlot(self):
            wantSBs = self.ui.lwSBs.selectedItems()
            print("want sbs", wantSBs)
            print([str(ii.text()) for ii in wantSBs])
            wantIdxs = [np.where(wantedSBs==int(str(ii.text())))[0][0] for ii in wantSBs]

            self.ui.gPlot.plotItem.clearPlots()

            if self.ui.rbxreal.isChecked():
                fx = np.real
            else:
                fx = np.imag

            xp = str(self.ui.lwx.selectedItems()[0].text())
            # T++
            xp = [int(xp[1]=="-"), int(xp[2]=="-")]

            if self.ui.rbyreal.isChecked():
                fy = np.real
            else:
                fy = np.imag

            yp = str(self.ui.lwy.selectedItems()[0].text())
            # T++
            yp = [int(yp[1] == "-"), int(yp[2] == "-")]

            print("xp, yp", xp, yp)

            for idx in wantIdxs:
                self.ui.gPlot.plot(
                    fx(self.data[xp[0], xp[1], idx, :]),
                    fy(self.data[yp[0], yp[1], idx, :]),
                    'o_'
                )
except ImportError:
    print("You'll need to install PyQt5 for interacting the MC data")

