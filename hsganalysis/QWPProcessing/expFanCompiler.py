import numpy as np
import os
import hsganalysis as hsg
from hsganalysis.jones import JonesVector as JV

class FanCompilerWithoutStokes(object):
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
        self._n = 1
        if negateNIR:
            print("WARNING: NEGATING NIR ALPHA")
            self._n = -1



    @staticmethod
    def fromDataFolder(folder, wantedSBs, keepErrors = False, negateNIR = True):
        """
        Create a fan compiler by passing the data path. Handles looping through the
        folder's sub-folders to find
        :param folder: The folder to search through. Alternatively, if it's a
        list/iterable, iterate through that instead. Useful if external code is
        directly removing sets of data.
        :return:
        """
        comp = FanCompiler(wantedSBs, keepErrors, negateNIR)
        # If it's a string, assume it's a single path that wants to be searached
        if isinstance(folder, str):
            wantFolders = hsg.natural_glob(folder, "*")
        else:
            # Otherwise, assume they've passed an iterable to search through
            wantFolders = folder


        for nirFolder in wantFolders:
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

        # nirAlpha = dataSet["alpha"][0][1]
        # nirGamma = dataSet["gamma"][0][1]
        if self._e:
            # Need to double to account for
            nirAlpha = [self._n*dataSet["alpha"][0][1], dataSet["alpha"][0][2]]
            nirGamma = [dataSet["gamma"][0][1], dataSet["gamma"][0][2]]
            for sb in self.want:
                # the list(*[]) bullshit is to make sure a list gets appended,
                # not a numpy array. Further complicated because if the list
                # comprehension returns nothing, it doesn't append anything,
                # hence casting to a list.
                newAData.append(list(*[ii[1:] for ii in alphaSet if ii[0] == sb]))
                newGData.append(list(*[ii[1:] for ii in gammaSet if ii[0] == sb]))
                newSData.append(list(*[ii[1:] for ii in s0Set if ii[0] == sb]))
                if not newAData[-1]: # no data was found.
                    newAData[-1] = [np.nan, np.nan]
                    newGData[-1] = [np.nan, np.nan]
                    newSData[-1] = [np.nan, np.nan]
        else:
            nirAlpha = [self._n*dataSet["alpha"][0][1]]
            nirGamma = [dataSet["gamma"][0][1]]
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

        # extending created lists accounts for keeping errors r not
        self.nirAlphas.extend(nirAlpha)
        self.nirGammas.extend(nirGamma)

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
        # self.arrA = wantedSBs.reshape(-1, 1) # so I can stack in the opposite direction
        # self.arrG = wantedSBs.reshape(-1, 1)  # so I can stack in the opposite direction
        # self.arrS = wantedSBs.reshape(-1, 1)  # so I can stack in the opposite direction

        keys = [ "S0", "S1", "S2", "S3", "alpha", "gamma", "DOP"]
        self.outputArrays = {ii: wantedSBs.reshape(-1, 1) for ii in keys}


        self.nirAlphas = []
        self.nirGammas = []
        self._e = keepErrors
        self._n = 1
        if negateNIR:
            print("WARNING: NEGATING NIR ALPHA")
            self._n = -1



    @staticmethod
    def fromDataFolder(folder, wantedSBs, keepErrors = False, negateNIR = True):
        """
        Create a fan compiler by passing the data path. Handles looping through the
        folder's sub-folders to find
        :param folder: The folder to search through. Alternatively, if it's a
        list/iterable, iterate through that instead. Useful if external code is
        directly removing sets of data.
        :return:
        """
        comp = FanCompiler(wantedSBs, keepErrors, negateNIR)
        # If it's a string, assume it's a single path that wants to be searached
        if isinstance(folder, str):
            wantFolders = hsg.natural_glob(folder, "*")
        else:
            # Otherwise, assume they've passed an iterable to search through
            wantFolders = folder


        for nirFolder in wantFolders:
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
        # newAData = []
        # newGData = []
        # newSData = []
        # alphaSet = dataSet["alpha"]
        # gammaSet = dataSet["gamma"]
        # s0Set = dataSet["S0"]

        newData = {ii: [] for ii in self.outputArrays}

        # nirAlpha = dataSet["alpha"][0][1]
        # nirGamma = dataSet["gamma"][0][1]
        if self._e:
            # Need to double to account for
            nirAlpha = [self._n*dataSet["alpha"][0][1], dataSet["alpha"][0][2]]
            nirGamma = [dataSet["gamma"][0][1], dataSet["gamma"][0][2]]
            for sb in self.want:
                # the list(*[]) bullshit is to make sure a list gets appended,
                # not a numpy array. Further complicated because if the list
                # comprehension returns nothing, it doesn't append anything,
                # hence casting to a list.
                [jj.append( list(*[ii[1:] for ii in dataSet[kk] if ii[0] == sb]) )
                 for kk, jj in newData.items()]
                if not newData["alpha"]: # no data was found.
                    for jj in newData.values():
                        jj[-1] = [np.nan, np.nan]
        else:
            nirAlpha = [self._n*dataSet["alpha"][0][1]]
            nirGamma = [dataSet["gamma"][0][1]]
            for sb in self.want:
                # print("dataSet", [list([ii[1] for ii in dataSet[kk] if ii[0] == sb] for\
                #         kk,
                #                                                    jj in newData.items()] )
                [jj.append( list([ii[1] for ii in dataSet[kk] if ii[0] == sb]) )
                 for kk, jj in newData.items()]
                if not newData["alpha"][-1]: # no data was found.
                    for jj in newData.values():
                        jj[-1] = [np.nan]

        try:
            for k, v in newData.items():
                self.outputArrays[k] = np.column_stack((self.outputArrays[k], v))
        except:
            # print(self.arrA.shape)
            # print(np.array(newAData).shape)
            raise

        # extending created lists accounts for keeping errors r not
        self.nirAlphas.extend(nirAlpha)
        self.nirGammas.extend(nirGamma)

    def build(self):
        data = self.buildAll()
        return data["alpha"], data["gamma"], data["S0"]


    def buildAll(self):
        fullData = {ii: np.append([-1], self.nirAlphas) for ii in
                    self.outputArrays.keys()}

        for k, v in self.outputArrays.items():
            fullData[k] = np.row_stack((fullData[k], v))

        # insert the gamma values into that array for the NIR laser
        fullData["gamma"][0, 1:] = self.nirGammas
        return fullData

    def buildAndSave(self, fname, *args, saveStokes=False):
        """
        fname: filename to save to. Must have a least one string formatter position
        to allow for saving separate alpha/gamma/s0 files. *args are passed to any
        other formatting positions.
        :param fname:
        :param args:
        :param saveStokes: Pass true if you want to save the stokes files
        :return:
        """

        if not os.path.exists(os.path.dirname(fname)):
            os.mkdir(os.path.dirname(fname))

        oh = "#\n" * 100
        oh += "\n\n"

        fullDataA, fullDataG, fullDataS = self.build()

        outputs = self.buildAll()

        if saveStokes:
            saveEms = [ii for ii in outputs.keys()]
        else:
            saveEms = ["alpha", "gamma", "S0"]

        for saveEm in saveEms:
            np.savetxt(fname.format(saveEm, *args), outputs[saveEm], header=oh,
                       delimiter=',',
                       comments='')

        # # fullData = np.append([-1], self.nirAlphas)
        # # fullData = np.row_stack((fullData, self.arrA))
        # np.savetxt(fname.format("alpha", *args), fullDataA, header=oh, delimiter=',',
        #            comments='')
        #
        # # fullData = np.append([-1], self.nirAlphas)
        # # fullData = np.row_stack((fullData, self.arrG))
        # np.savetxt(fname.format("gamma", *args), fullDataG, header=oh, delimiter=',',
        #            comments='')
        #
        # # fullData = np.append([-1], self.nirAlphas)
        # # fullData = np.row_stack((fullData, self.arrS))
        # np.savetxt(fname.format("S0", *args), fullDataS, header=oh, delimiter=',',
        #            comments='')

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