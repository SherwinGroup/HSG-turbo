import numpy as np
import os
import hsganalysis as hsg
from hsganalysis.jones import JonesVector as JV

class FanCompiler(object):
    """
    Helper class for compiling the data of a polarimetry NIR alpha sweep

    Typical use scenario:

    datasets = [ ... list of folders, each being a dataset of different NIR alphas ...]
    outputs = FanComilier(<whatever sideband orders you want compiled>)
    for data in datasets:
        laserParams, rawData = hsg.hsg_combine_qwp_sweep(folder, save=False, verbose=False)
        _, fitDict = hsg.proc_n_fit_qwp_data(rawData, laserParams, vertAnaDir="VAna" in folder,
                                        series=folder)
        outputs.addSet(nira, fitDict)
    outputs.buildAndSave(fname)

    """
    def __init__(self, wantedSBs, keepErrors = False, negateNIR = True):
        """

        :param wantedSBs:
        :param keepErrors:
        :param negateNIR: flag for whether to negate the NIR alpha value. Currently,
        this is done because the PAX views -z direction, while home-built views +z (
        with NIR)
        """
        self.want = np.array(wantedSBs)
        self.arrA = wantedSBs.reshape(-1, 1) # so I can stack in the opposite direction
        self.arrG = wantedSBs.reshape(-1, 1)  # so I can stack in the opposite direction
        self.arrS = wantedSBs.reshape(-1, 1)  # so I can stack in the opposite direction
        self.nirAlphas = []
        self.nirGammas = []
        self._e = keepErrors
        self._n = negateNIR

    @staticmethod
    def fromDataFolder(folder, wantedSBs, keepErrors = False, negateNIR = True):
        """
        Create a fan compiler by passing the data path. Handles looping through the
        folder's sub-folders to find
        :param folder:
        :return:
        """
        comp = FanCompiler(wantedSBs, keepErrors, negateNIR)
        for nirFolder in hsg.natural_glob(folder, "*"):
            if "skip" in nirFolder.lower(): continue
            laserParams, rawData = hsg.hsg_combine_qwp_sweep(nirFolder, save=False,
                                                             verbose=False,
                                                             loadNorm=False)

            _, fitDict = hsg.proc_n_fit_qwp_data(rawData, laserParams,
                                                 vertAnaDir="VAna" in nirFolder,
                                                 series=nirFolder)
            comp.addSet(fitDict)
        return comp

    def addSet(self, dataSet):
        """ Assume it's passed from  proc_n_fit_qwp_data"""
        newAData = []
        newGData = []
        newSData = []
        alphaSet = dataSet["alpha"]
        gammaSet = dataSet["gamma"]
        s0Set = dataSet["S0"]

        nirAlpha = dataSet["alpha"][0][1]
        nirGamma = dataSet["gamma"][0][1]
        if self._e:
            for sb in self.want:
                newAData.append([ii[1:] for ii in alphaSet if ii[0] == sb])
                newGData.append([ii[1:] for ii in gammaSet if ii[0] == sb])
                newSData.append([ii[1:] for ii in s0Set if ii[0] == sb])
                if not newAData[-1]: # no data was found.
                    newAData[-1] = [np.nan, np.nan]
                    newGData[-1] = [np.nan, np.nan]
                    newSData[-1] = [np.nan, np.nan]
        else:
            for sb in self.want:
                newAData.append([ii[1] for ii in alphaSet if ii[0] == sb])
                newGData.append([ii[1] for ii in gammaSet if ii[0] == sb])
                newSData.append([ii[1] for ii in s0Set if ii[0] == sb])
                if not newAData[-1]: # no data was found.
                    newAData[-1] = [np.nan]
                    newGData[-1] = [np.nan]
                    newSData[-1] = [np.nan]

        try:
            self.arrA = np.column_stack((self.arrA, newAData))
            self.arrG = np.column_stack((self.arrG, newGData))
            self.arrS = np.column_stack((self.arrS, newSData))
        except:
            print(self.arrA.shape)
            print(np.array(newAData).shape)
            raise


        if self._n:
            print("WARNING: Negating NIR Alpha")
            self.nirAlphas.append(-nirAlpha)
        else:
            self.nirAlphas.append(nirAlpha)
        self.nirGammas.append(nirGamma)

    def build(self):
        fullDataA = np.append([-1], self.nirAlphas)
        fullDataA = np.row_stack((fullDataA, self.arrA))

        fullDataG = np.append([-1], self.nirGammas)
        fullDataG = np.row_stack((fullDataG, self.arrG))

        fullDataS = np.append([-1], self.nirAlphas)
        fullDataS = np.row_stack((fullDataS, self.arrS))

        return fullDataA, fullDataG, fullDataS

    def buildAndSave(self, fname, *args):
        """
        fname: filename to save to. Must have a least one string formatter position
        to allow for saving separate alpha/gamma/s0 files. *args are passed to any
        other formatting positions.
        :param fname:
        :param args:
        :return:
        """

        if not os.path.exists(os.path.dirname(fname)):
            os.mkdir(os.path.dirname(fname))

        oh = "#\n" * 100
        oh += "\n\n"

        fullDataA, fullDataG, fullDataS = self.build()

        # fullData = np.append([-1], self.nirAlphas)
        # fullData = np.row_stack((fullData, self.arrA))
        np.savetxt(fname.format("alpha", *args), fullDataA, header=oh, delimiter=',',
                   comments='')

        # fullData = np.append([-1], self.nirAlphas)
        # fullData = np.row_stack((fullData, self.arrG))
        np.savetxt(fname.format("gamma", *args), fullDataG, header=oh, delimiter=',',
                   comments='')

        # fullData = np.append([-1], self.nirAlphas)
        # fullData = np.row_stack((fullData, self.arrS))
        np.savetxt(fname.format("S0", *args), fullDataS, header=oh, delimiter=',',
                   comments='')

class SbStateGetter(object):
    """
    sister function to FanCompiler. Useful for taking the arrays out of a FanCompiler
    allowing indexing based on sb number and nir alpha angles.

    Example creation:
    fc = FanCompiler(wantSBs)
    // initialize fc

    getter = SbStateGetter(
        fc.arrA[:, 1:],
        fc.arrG[:, 1:],
        fc.want,
        fc.nirAlphas
    )
    """
    def __init__(self, alphas, gammas, sbNum, niralphas):
        self.alphas = alphas
        self.gammas = gammas
        self.sbNum = sbNum
        self.niralphas = niralphas

        self.invS = {kk: ii for ii, kk in enumerate(sbNum)}

        self.invA = {kk: ii for ii, kk in enumerate(niralphas)}

    def getSBState(self, sb, nirA):
        alpha = self.alphas[self.invS[sb], self.invA[nirA]]
        gamma = self.gammas[self.invS[sb], self.invA[nirA]]

        return alpha, gamma

    def __call__(self, sb, nirA):
        return self.getSBState(sb, nirA)

    def getStateDict(self, sb, nirA):
        a, g = self.getSBState(sb, nirA)
        return {"alpha": a, "gamma": g}

def jonesToFans(sbs, J, wantNIR = np.arange(-90, 90, 5)):
    """
    Take a Jones matrix as an Nx2x2
    :param J:
    :return:
    """
    vec = JV(alpha=wantNIR, gamma=0)
    vec.apply_transformation(J)


    alphas = np.column_stack((sbs, vec.alpha))
    alphas = np.row_stack(([-1] + wantNIR.tolist(), alphas ))

    gammas = alphas.copy()
    gammas[1:, 1:] = vec.gamma


    return alphas, gammas