import numpy as np
import hsganalysis.ipg as pg
import matplotlib.pylab as plt
import glob
import os
import json
import hsganalysis as hsg
from hsganalysis.ipg import PolarImagePlot
from PyQt5 import QtCore, QtGui, QtWidgets
from scipy.optimize import minimize
from hsganalysis import newhsganalysis
from hsganalysis import JonesVector as JV
newhsganalysis.plt = pg
np.set_printoptions(linewidth=400)

from .createFanDiagram import createFan
from .extractMatrices import *
from .expFanCompiler import *
from .polarPlot import *

## For the +54deg data
# alpha54FullDatas = np.genfromtxt(alpha54File, delimiter=',')
# alpha54SBs = alpha54FullDatas[1:,0]
# alpha54AlphaNIRS = alpha54FullDatas[0,1:]
# alpha54Data = alpha54FullDatas[1:, 1:]
#
# gamma54FullDatas = np.genfromtxt(gamma54File, delimiter=',')
# gamma54SBs = gamma54FullDatas[1:,0]
# gamma54AlphaNIRS = gamma54FullDatas[0,1:]
# gamma54Data = gamma54FullDatas[1:, 1:]
#
# sbs54 = np.array(alpha54SBs, dtype=int)
# nirAlphas54 = np.array(gamma54AlphaNIRS, dtype=int)
#
# para54Data = np.genfromtxt(para54Spec, delimiter=',')[3:]
# perp54Data = np.genfromtxt(perp54Spec, delimiter=',')[3:]
#
# para54Data = para54Data[(8<=para54Data[:,0]) & (para54Data[:,0]<=36)]
# perp54Data = perp54Data[(8<=perp54Data[:,0]) & (perp54Data[:,0]<=36)]
# intensity54Data = np.column_stack((para54Data[:,0], perp54Data[:, 3]/para54Data[:, 3]))
#
#
# ## For the -3deg data
# alpham3FullDatas = np.genfromtxt(alpham3File, delimiter=',')
# alpham3SBs = alpham3FullDatas[2:,0] # cut out 6th order with [2:...] slice
# alpham3AlphaNIRS = alpham3FullDatas[0,1:]
# alpham3Data = alpham3FullDatas[2:, 1:] # cut out 6th order with [2:...] slice
#
# gammam3FullDatas = np.genfromtxt(gammam3File, delimiter=',')
# gammam3SBs = gammam3FullDatas[2:,0] # cut out 6th order with [2:...] slice
# gammam3AlphaNIRS = gammam3FullDatas[0,1:]
# gammam3Data = gammam3FullDatas[2:, 1:] # cut out 6th order with [2:...] slice
#
#
# sbsm3 = np.array(alpham3SBs, dtype=int)
# nirAlphasm3 = np.array(gammam3AlphaNIRS, dtype=int)
#
# param3Data = np.genfromtxt(param3Spec, delimiter=',')[3:]
# perpm3Data = np.genfromtxt(perpm3Spec, delimiter=',')[3:]
#
# param3Data = param3Data[(8<=param3Data[:,0]) & (param3Data[:,0]<=36)]
# perpm3Data = perpm3Data[(8<=perpm3Data[:,0]) & (perpm3Data[:,0]<=36)]
# intensitym3Data = np.column_stack((param3Data[:,0], perpm3Data[:, 3]/param3Data[:, 3]))


## This attempt is going to try and reuse the code that was used for the DBR paper
## It uses some system of equations that come from doing the matrix multipliations
## out and stuff
##

totalcount = 0

def unflattenJ(jVec):
    j =  (jVec[:3]+1j*jVec[3:])
    j = np.append([1], j)
    return j.reshape(2,2)

def solver(r, J):
    global totalcount
    totalcount+=1
    # print "total runs", totalcount
    # print r.shape
    # print J.shape
    # N = len(r)//2
    # rVec = (r[:N]+1j*r[N:]).reshape(-1, 2)
    # print "in vector", rVec
    Jxy, Jyx, Jyy = (J[:3]+1j*J[3:])
    nir, sb = r
    eN = np.exp(1j*np.deg2rad(nir.delta))
    eH = np.exp(-1j*np.deg2rad(sb.delta))
    cotH = 1./np.tan(np.deg2rad(sb.phi))
    cN = np.cos(np.deg2rad(nir.phi))
    sN = np.sin(np.deg2rad(nir.phi))

    return cotH*eH*(Jyx*cN + Jyy*sN*eN)-cN-Jxy*sN*eN

def findJBad(obj):
    mod = QtWidgets.QApplication.keyboardModifiers()
    if mod & QtCore.Qt.ShiftModifier:
        print("Skipping")
        return
    alphaFullDatas = obj.alphaFullDatas
    alphaSBs = obj.sbs
    alphaAlphaNIRS = obj.nirAlphas
    gammaSBs = obj.sbs
    gammaAlphaNIRS = obj.nirAlphas
    sbs = obj.sbs
    nirAlphas = obj.nirAlphas
    cboxGroup = obj.cboxGroup
    sbGetter = obj.sbGetter
    palp = obj.palp
    pgam = obj.pgam
    intRatioCurve = obj.intRatioCurve


    outputAlphaData = np.empty((alphaSBs.shape[0]+1, alphaAlphaNIRS.shape[0]+1)) * \
                      np.nan
    outputAlphaData[1:, 0] = alphaSBs
    outputAlphaData[0, 1:] = alphaAlphaNIRS

    outputGammaData = np.empty((gammaSBs.shape[0] + 1, gammaAlphaNIRS.shape[0] + 1)) * \
                      np.nan
    outputGammaData[1:, 0] = gammaSBs
    outputGammaData[0, 1:] = gammaAlphaNIRS
    outputJMatrix = np.empty((len(sbs),9))

    wantNIRAlphas = [ii for ii in nirAlphas if cboxGroup.button(int(ii)).isChecked()]
    wantNIRIndices = np.array([nirAlphas.tolist().index(ii) for ii in wantNIRAlphas])

    for idx, sb in enumerate(sbs):
        als, gms = zip(*[sbGetter(sb, ii) for ii in wantNIRAlphas])
        sbJones = JV(alpha=als, gamma=gms)
        nirJones = JV(alpha=wantNIRAlphas, gamma=0)

        costfunc = lambda jmatrix: np.linalg.norm(solver([nirJones, sbJones], jmatrix))

        p = minimize(costfunc, x0=np.ones(6))
        J = unflattenJ(p.x)


        ### FOR ONLY SHOWING WHAT WAS USED IN FIT
        # nirJones.apply_transformation(J)
        # outputAlphaData[idx + 1, wantNIRIndices + 1] = nirJones.alpha
        # outputGammaData[idx + 1, wantNIRIndices + 1] = nirJones.gamma

        ### FOR ONLY SHOWING EVERYTHING
        nirJones = JV(alpha=nirAlphas, gamma=0)
        nirJones.apply_transformation(J)

        outputAlphaData[idx+1, 1:] = nirJones.alpha
        outputGammaData[idx+1, 1:] = nirJones.gamma

        np.set_printoptions(formatter={"float_kind":lambda x:"{: 6.2f}".format(x)})
        outputJMatrix[idx] = np.array([sb, 1] + p.x[:3].tolist() + [0] + p.x[3:].tolist())

    palp.setImage(outputAlphaData[1:, 1:])
    palp.setLevels(-90, 90)
    pgam.setImage(outputGammaData[1:, 1:])
    pgam.setLevels(-45, 45)
    palp.imageItem.render()

    # header = "#{}\n".format(
    #      [ii for ii in nirAlphas if cboxGroup.button(int(ii)).isChecked()]
    # )+"#\n"*99 + "\n\n"
    # np.savetxt("BackedOutFromJones_alpha.txt", outputAlphaData, delimiter=',',
    #            header=header, comments='', fmt="%0.6e")
    # np.savetxt("BackedOutFromJones_gamma.txt", outputGammaData, delimiter=',',
    #            header=header, comments='', fmt="%0.6e")
    # np.savetxt("JonesMatrix.txt", outputJMatrix, delimiter=',',
    #            header=header, comments='', fmt="%0.6e")

    intRatioCurve.setData(alphaSBs, np.sqrt(
        np.abs(outputJMatrix[:, 4]**2+1j*outputJMatrix[:, -1]**2))
                          )

    obj.outputJMatrix = outputJMatrix

def updateTCurves():
    try:
        J54 = mainwid54.outputJMatrix
        Jm3 = mainwidm3.outputJMatrix
    except (AttributeError, NameError):
        return

    T54 = makeT(J54, 54)
    Tm3 = makeT(Jm3, -3)

    mainwid54.TMatrix = T54
    mainwidm3.TMatrix = Tm3

    Tm3pm = 180 / 3.14159 * np.angle(Tm3[0, 1, :] / Tm3[1, 0, :])
    T54pm = 180 / 3.14159 * np.angle(T54[0, 1, :] / T54[1, 0, :])

    print("T-3 +- {}".format(Tm3pm[:5].mean()))
    print("T54 +- {}".format(T54pm[:5].mean()))
    m3TpmTmp.setData(sbsm3, Tm3pm)
    p54TpmTmp.setData(sbsm3, T54pm)
    diffTpmTmp.setData(sbsm3, Tm3pm-T54pm)

    Tm3pp = 180 / 3.14159 * np.angle(Tm3[0, 0, :] / Tm3[1, 1, :])
    T54pp = 180 / 3.14159 * np.angle(T54[0, 0, :] / T54[1, 1, :])
    m3TppTmm.setData(sbsm3, Tm3pp)
    p54TppTmm.setData(sbsm3, T54pp)

    Tm3pp = np.abs(Tm3[0, 0, :])
    Tm3mm = np.abs(Tm3[1, 1, :])
    T54pp = np.abs(T54[0, 0, :])
    T54mm = np.abs(T54[1, 1, :])
    m3TppMag.setData(sbsm3, Tm3pp)
    m3TmmMag.setData(sbsm3, Tm3mm)
    p54TppMag.setData(sbsm3, T54pp)
    p54TmmMag.setData(sbsm3, T54mm)


def makeInteractiveFanWidget(compiledAlpha, compiledGamma, crystalOrientation,
               intensityData = [],
               plotFanNIRs = np.arange(-90, 90, 5),
               sliceInExp = True,
               calculationCallback = lambda J, T: T,
                **kwargs):
    """

    :param compiledAlpha: saved files created from a FanCompiler.buildAndSave()
    :param compiledGamma:
    :param calculationCallback: function to be called when T/J matrices are recalculated
    :return:
    """
    NMax = kwargs.get("NMax", 12)

    nirAlphas = compiledAlpha[0, 1:]
    sbs = compiledAlpha[1:, 0]
    print(f"Found sbs: {sbs}")
    print(f"found nira: {nirAlphas}")
    alphaData = compiledAlpha[1:, 1:]
    gammaData = compiledGamma[1:, 1:]

    mainwid = QtWidgets.QWidget()
    layout = QtWidgets.QVBoxLayout()
    cboxesLayout = QtWidgets.QGridLayout()
    cboxGroup = QtWidgets.QButtonGroup()
    cboxGroup.setExclusive(False)

    for idx, ang in enumerate(nirAlphas):
        cb = QtWidgets.QCheckBox(str(ang))
        cb.setChecked(True)
        cboxGroup.addButton(cb, id=int(ang))
        cboxesLayout.addWidget(cb, idx // 12, idx % 12)

    layout.addLayout(cboxesLayout)
    tabWid = QtWidgets.QTabWidget(mainwid)

    palp, pgam = createFan(plotFanNIRs, sbs)
    # from hsganalysis.QWPProcessing.fanDiagram import FanDiagram


    def updateJ():

        mod = QtWidgets.QApplication.keyboardModifiers()
        if mod & QtCore.Qt.ShiftModifier:
            print("Skipping")
            return
        # Which NIR alpha angles you've selected to include in the J matrix calculation
        wantNIRAlphas = [ii for ii in nirAlphas if
                         cboxGroup.button(int(ii)).isChecked()]
        # And converting those to indices. Include the 0th order to include the SB
        wantNIRIndices = [0] + \
                         [nirAlphas.tolist().index(ii)+1 for ii in wantNIRAlphas]

        toFitAlphas = compiledAlpha[:, wantNIRIndices]
        # print(f"Fitting to: \n{toFitAlphas}")
        toFitGammas = compiledGamma[:, wantNIRIndices]

        J = findJ(toFitAlphas, toFitGammas)
        reconstructedAlpha, reconstructedGamma = jonesToFans(
            sbs, J, wantNIR=plotFanNIRs)
        # print(f"Recon to: \n{reconstructedAlpha}")

        if sliceInExp:
            for idx, _ in enumerate(nirAlphas):

                niralpha = compiledAlpha[0, idx+1]
                if np.abs(compiledGamma[0, idx+1])>1: continue # finite gamma don't make
                # sense in this plot
                try:
                    reconstructedAlpha[1:,
                        np.argwhere(reconstructedAlpha[0, :].astype(int) == niralpha)[0][0]
                    ] = compiledAlpha[1:, idx+1]

                    reconstructedGamma[1:,
                        np.argwhere(reconstructedGamma[0, :].astype(int) == niralpha)[0][0]
                    ] = compiledGamma[1:, idx+1]
                except IndexError:
                    print("Warning! Unable to slice in NIR alpha = {}".format(niralpha))

        palp.setImage(reconstructedAlpha[1:, 1:])
        pgam.setImage(reconstructedGamma[1:, 1:])

        # palp.ui.histogram.setHistogramRange(-90, 90)
        palp.setLevels(-90, 90)
        # pgam.ui.histogram.setHistogramRange(-45, 45)
        pgam.setLevels(-45, 45)


        T = makeT(J, crystalOrientation)



        TppPolar.setData(sbs[:NMax]+np.abs(T[0, 0, :NMax]), np.angle(T[0, 0, :NMax],
                                                                 deg=False))
        TpmPolar.setData(sbs[:NMax]+np.abs(T[0, 1, :NMax]), np.angle(T[0, 1, :NMax],
                                                                 deg=False))
        TmpPolar.setData(sbs[:NMax]+np.abs(T[1, 0, :NMax]), np.angle(T[1, 0, :NMax],
                                                                 deg=False))
        TmmPolar.setData(sbs[:NMax]+np.abs(T[1, 1, :NMax]), np.angle(T[1, 1, :NMax],
                                                                deg=False))

        TppLinear.setData(sbs[:NMax], np.abs(T[0, 0, :NMax]))
        TppALinear.setData(sbs[:NMax], np.angle(T[0, 0, :NMax],deg=True))
        TpmLinear.setData(sbs[:NMax], np.abs(T[0, 1, :NMax]))
        TpmALinear.setData(sbs[:NMax], np.angle(T[0, 1, :NMax],deg=True))
        TmpLinear.setData(sbs[:NMax], np.abs(T[1, 0, :NMax]))
        TmpALinear.setData(sbs[:NMax], np.angle(T[1, 0, :NMax],deg=True))
        TmmLinear.setData(sbs[:NMax], np.abs(T[1, 1, :NMax]))
        TmmALinear.setData(sbs[:NMax], np.angle(T[1, 1, :NMax],deg=True))

        TppoTmmLinear.setData(sbs[:NMax], np.abs(T[0, 0, :NMax] / T[1, 1, :NMax]))
        TppoTmmALinear.setData(sbs[:NMax], np.angle(T[0, 0, :NMax] / T[1, 1, :NMax],
                               deg=True))

        TpmoTmpLinear.setData(sbs[:NMax], np.abs(T[0, 1, :NMax] / T[1, 0, :NMax]))
        TpmoTmpALinear.setData(sbs[:NMax], np.angle(T[0, 1, :NMax] / T[1, 0, :NMax],
                               deg=True))



        calculationCallback(J, T)
    mainwid.updateJ = updateJ # need to keep a reference


    tabWid.addTab(palp, "Fan Diagram")

    # intRatioWid = pg.plot(intRatioWidtensityData, fmt="ko-")
    # intRatioWid = pg.PlotContainerWindow()
    # intRatioCurve = intRatioWid.plot(intensityData, fmt="ko-")
    # intRatioCurve = intRatioWid.plot(fmt="ro-")
    #
    # tabWid.addTab(intRatioWid, "Ratio")


    TPolars = makePolarPlot()
    pg.legend()
    tabWid.addTab(TPolars, "T Matrices (Polar)")

    TppPolar = polarPlot('ko-', name="T++")
    TpmPolar = polarPlot('ro-', name="T+-")
    TmpPolar = polarPlot('bo-', name="T-+")
    TmmPolar = polarPlot('mo-', name="T--")

    double = pg.DoubleYPlot()
    TLinears = pg.PlotContainerWindow(plotWidget=double)
    TLinears.addLegend()
    TppLinear = TLinears.plot('ko-', name="T++")
    TpmLinear = TLinears.plot('ro-', name="T+-")
    TmpLinear = TLinears.plot('bo-', name="T-+")
    TmmLinear = TLinears.plot('mo-', name="T--")

    TppALinear = TLinears.plot('ko--', name="T++")
    double.p2.addItem(TppALinear)
    TpmALinear = TLinears.plot('ro--', name="T+-")
    double.p2.addItem(TpmALinear)
    TmpALinear = TLinears.plot('bo--', name="T-+")
    double.p2.addItem(TmpALinear)
    TmmALinear = TLinears.plot('mo--', name="T--")
    double.p2.addItem(TmmALinear)

    tabWid.addTab(TLinears, "T Matrices (Linear)")






    double2 = pg.DoubleYPlot()
    TLinears2 = pg.PlotContainerWindow(plotWidget=double2)
    TLinears2.addLegend()
    TppoTmmLinear = TLinears2.plot('ko-', name="|T++/T--|")
    TpmoTmpLinear = TLinears2.plot('ro-', name="|T+-/T+-|")

    TppoTmmALinear = TLinears2.plot('ko--', name="ph(T++/T--)")
    double2.p2.addItem(TppoTmmALinear)
    TpmoTmpALinear = TLinears2.plot('ro--', name="ph(T+-/T-+)")
    double2.p2.addItem(TpmoTmpALinear)


    tabWid.addTab(TLinears2, "T Matrices (Ratios)")





    # mainwid.findJ = lambda :findJ(mainwid)
    cboxGroup.buttonToggled.connect(updateJ)


    layout.addWidget(tabWid)
    mainwid.setLayout(layout)

    return mainwid
