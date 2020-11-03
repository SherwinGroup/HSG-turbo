import numpy as np
from scipy.optimize import minimize
import hsganalysis as hsg
import hsganalysis.ipg as pg
from hsganalysis.jones import JonesVector as JV
from hsganalysis.QWPProcessing.expFanCompiler import *
from hsganalysis.QWPProcessing.fanDiagram import *


def rtd(a):
    # radians to degrees
    return a*180/np.pi

def dtr(a):
    # degrees to radians
    return a*np.pi/180

def makeU(t):
    # unitary rotation matrix between x/y and sigma+- bases
    # t is the crystal rotation angle. By note below, it should be relative to [011]
    # t = dtr(t-45) # +45 because I'm relative to [011], Matrix wants [010]
    #             # I SHOULD DOUBLE CHECK THAT IT'S +45 NOT -45

    # no, fuck it, I'm tired of thinking of angles right now. Pass the damn angle relative
    # to [010], sort it out yourself.
    t = dtr(t)
    a =    np.exp(-1j * t)
    b =   -np.exp( 1j * t)
    c = 1j*np.exp(-1j * t)
    d = 1j*np.exp( 1j * t)
    mat = np.array([[a, b],[c,d]])/np.sqrt(2)
    return mat

def makeT(J, ang):
    """take a J matrix as saved from processing,
    return a fully shaped T matrix"""

    if J.ndim == 3:
        ...
    else:
        sbs = J[:, 0].astype(int)
        J = J[:, 1:]
        # Shape it correctly
        J = J[:, :4] + 1j * J[:, 4:]
        J = J.reshape(-1, 2, 2)
    # J = J.transpose(1, 2, 0)

    U = makeU(ang)
    T = np.einsum("ij,jlx,lm->imx", U.conj().T, J, U)

    return T

def makeJfromT(T, ang):
    """take a T matrix as saved from processing,
    return a fully shaped J matrix"""

    U = makeU(ang)
    J = np.einsum("ij,jlx,lm->imx", U, T, U.conj().T)

    return J


def unflattenJ(jVec):
    jVec = np.array(jVec)
    j =  (jVec[:3]+1j*jVec[3:])
    j = np.append([1], j)
    return j.reshape(2,2)

def solver(r, J):
    """ Reminder:
    a Jones Vector is represented by
    (cos(phi), sin(phi)exp(i delta)).T
    """
    Jxy, Jyx, Jyy = (J[:3]+1j*J[3:])
    nir, sb = r
    eN = np.exp(1j*np.deg2rad(nir.delta))
    eH = np.exp(-1j*np.deg2rad(sb.delta))
    cH = np.cos(np.deg2rad(sb.phi))
    sH = np.sin(np.deg2rad(sb.phi))
    cN = np.cos(np.deg2rad(nir.phi))
    sN = np.sin(np.deg2rad(nir.phi))

    # cotH = 1./np.tan(np.deg2rad(sb.phi))
    # return cotH*eH*(Jyx*cN + Jyy*sN*eN)-cN-Jxy*sN*eN

    return cH*eH*(Jyx*cN + Jyy*sN*eN) - sH*(cN + sN*eN*Jxy)

def findJ(alphas, gammas=None, **kwargs):
    """
    Extract the Jones matrix (x/y basis) from given data.
    alphas/gammas should be the matrices saved from the FanCompiler, of form:

    arb     | niralpha1 | niralpha2 | niralpha3 | niralpha4 | ...
    1st sb  | 1sb alpha | 1sb alpha |     .
    2nd sb  | 2sb alpha | 2sb alpha |     .
    3rd sb  | 3sb alpha | 3sb alpha |     .
      .
      .
      .

    Assumes both alphas/gammas are the same shape

    kwarg options:
       returnFlat-- return a flattened (Nx9) Jones matrix, of form
         [[sb#, Re(Jxx), Re(Jxy), Re(Jyx), Re(Jyy), Im(Jxx), Im(Jxy), Im(Jyx), Im(Jyy)],
          [  .. ]
          ]
        useful for saving to file.
        If false, return an 2x2xN,
          [[[Jxx, Jxy],
            [Jyx, Jyy]],
            [[ .. ]]
          ]
        useful for continuing into processing (and JonesVector maths).

        NOTE: You probably shouldn't use the "Return Flat" argument for saving.
        Instead, get the J matrix back and use saveT() to avoid accidentally
        introducing errors from difference in the order of the columns of the flattened
        matrix in this function vs saveT/loadT

    You can also just pass a FanCompiler object and it'll pull the alpha/gamma from
    that.

    :return:
    """

    defaults = {
        "returnFlat": False
    }
    defaults.update(kwargs)

    if gammas is None and isinstance(alphas, FanCompiler):
        alphas, gammas, _ = alphas.build(withErrors=False)

    sbs = alphas[1:,0]
    nirAlphas = alphas[0, 1:]
    nirGammas = gammas[0, 1:]
    # This SbStateGetter makes it more convenient to get the alpha and gamma
    # angles for a specific sideband and NIR alpha
    sbGetter = SbStateGetter(alphas[1:, 1:], gammas[1:, 1:], sbs, nirAlphas)

    ## Initialize the matrices
    outputFlatJMatrix = np.empty((len(sbs),9))
    outputJMatrix = np.empty((2, 2, len(sbs)), dtype=complex)

    # There's one J matrix for each sideband order, so naturally have to iterate over
    # each sideband
    for idx, sb in enumerate(sbs):
        # Get the list of alpha and gamma angles for each of the NIR Alphas used
        als, gms = zip(*[sbGetter(sb, ii) for ii in nirAlphas])
        # Check to make sure all of the data is reasonable (not nan, meaning the sideband
        # wasn't observed for all NIRalpha, or infinite when the fits fucked up)
        # Leaving these in causes issues for the minimizer, so they have to be skipped
        if not any(np.isfinite(als)) or not any(np.isfinite(gms)):
            # Do some python magic so I can still use p.x further and not have to
            # wrap everything in a try/except
            p = type("_", (object, ), {"x": np.array([np.nan]*3 + [0]*3)})
        else:
            sbJones = JV(alpha=als, gamma=gms)
            nirJones = JV(alpha=nirAlphas, gamma=nirGammas)
            # Note: We current'ly don't do error propagation at this step
            costfunc = lambda jmatrix: np.linalg.norm(solver([nirJones, sbJones], jmatrix))

            p = minimize(costfunc, x0=np.ones(6))
        J = unflattenJ(p.x)
        outputJMatrix[..., idx] = J

        outputFlatJMatrix[idx] = np.array([sb, 1] # Re(Jxx) === 1
                                          + p.x[:3].tolist() # Re(Jxy-Jyy)
                                          + [0]  # Im(Jxx) === 0
                                          + p.x[3:].tolist() # Im(Jxy-Jyy)
                                          )

    if defaults["returnFlat"]: return outputFlatJMatrix
    return outputJMatrix


def saveT(T, sbs, out):
    """
    Save a complex T matrix, input as an Nx2x2, into a text file. Dumps it as a CSV
    where the first four columns are the real components, the last four are imaginary
    :param T:
    :param out:
    :return:
    """
    T = np.array(T.transpose(2, 0, 1))

    ## I'm nervous of trusting how numpy handles .view() on complex types. I feel like
    # I've seen it swap orders or something, where I've had to change the loadT function
    # to compensate. I guess when in doubt, process data from scratch, save it and
    # reload it and make sure the memory and disk matrices agree.

    # 01/04/19 No, fuck this. I don't trust view at all. I'm looking at two different
    # T matrices, and in one instance this gets reordered as
    #     ReT++,ReT+-,ReT-+,ReT--,ImT++,ImT+-,ImT-+,ImT--
    # while in another, it does it as
    #     ReT++,ImT++,ReT+-,ImT+-,ReT-+,ImT-+,ReT--,ImT--
    #
    # I have no fucking clue why it does it that way, but I'm sick and fucking tired of it
    # So no more
    # 12/02/19 - How do you really feel Darren?
    #
    # flatT = T.reshape(-1, 4).view(float).reshape(-1, 8)

    flatT = T.reshape(-1, 4)
    flatT = np.column_stack((flatT.real, flatT.imag))

    # I'm also going to complicate this, because I want to save it like qile's matlab
    # code save his files, so that we can use the same loading.
    # As of 12/19/18, I believe the above code should be ordering columns as,

    ###   0    1     2     3      4     5    6     7
    ### ReT++,ReT+-,ReT-+,ReT--,ImT++,ImT+-,ImT-+,ImT--

    # Qile saves as
    ###   0    1     2     3      4     5    6     7
    ### ReT--,ImT--,ReT+-,ImT+-,ReT-+,ImT-+,ReT++,ImT++

    reorder = [ 3, 7, 1, 5, 2, 6, 0, 4 ]

    flatT = flatT[:, reorder]



    flatT = np.column_stack((sbs, flatT))

    header = "SB,ReT++,ImT++,ReT+-,ImT+-,ReT-+,ImT-+,ReT--,ImT--"
    header = "SB,ReT++,ReT+-,ReT-+,ReT--,Im++,Im+-,Im-+,Im--"

    header = "SB,ReT--,ImT--,ReT+-,ImT+-,ReT-+,ImT-+,ReT++,ImT++"
    np.savetxt(out,flatT, header=header, comments='', delimiter=',',)
    print("saved {}\n".format(out))

def loadT(name):
    """
    Load the file saved by saveT
    :param name:
    :return:
    """
    d = np.genfromtxt(name, delimiter=',')[1:] # label line

    sbs = d[:,0]

    T = d[:, 1:5]+1j*d[:, 5:]
    # T = d[:, 1::2] + 1j * d[:, 2::2]
    T = T.reshape(-1, 2, 2)
    T = T.transpose(1, 2, 0)


    T = d[:, 1::2] + 1j * d[:, 2::2]

    ## Reversing the axis gets back into
    # T++, T+-, T-+, T--, which makes it easier to then
    #             reshape it to appropriate 2x2 matrices
    #  and then                      put them in the correct indices.
    T = T[:, ::-1].reshape(-1, 2, 2).transpose(2, 1, 0)
    # T = T.reshape(-1, 2, 2)


    return T, sbs


def fan_n_Tmat(file,observedSidebands,crystalAngle,saveFileName,save_results=False,Plot=True):
    """
    This function will take a folder of polarimetry scans and produce
    DOP, alpha/gamma, J and T matrices, Fan, and Matrix Relations

    :param file: String of folder name containing 4 polarimetry scans
    :param observedSidebands: np array of observed sidebands. Data will be
        cropped such that these sidebands are included in everything.
    :param crystalAngle: (Float) Angle of the sample from the 010 crystal face
    :saveFileName: Str of what you want to call the text files to be saved
    :save_results: Boolean controls if things are saved to txt files.
        Currently saves DOP, alpha, gamma, J matrix, T matrix, Fan, and Matrix
        Relations
    :Plot: Boolean controls if plots are displayed

    :return: DOP,alpha,gamma,J matrix, T matrix, Matrix Relations
    """

    # Make a list of folders for the fan data
    datasets = hsg.natural_glob(file, "*")

    # FanCompiler class keeps alpha and gamma data together for different pols
    fanData = FanCompiler(observedSidebands)

    # Initialize arrays for DOP
    dop00 = np.array([])
    dop45 = np.array([])
    dop90 = np.array([])
    dopn45 = np.array([])

    for data in datasets:
        laserParams, rawData = hsg.hsg_combine_qwp_sweep(data)
        _, fitDict = hsg.proc_n_fit_qwp_data(rawData,vertAnaDir = False, laserParams = laserParams,wantedSbs = observedSidebands)
        # Add the new alpha and gamma sets to the fan class
        fanData.addSet(fitDict)

        # Get Stoke Parameters
        s0 = fitDict['S0']
        s1 = fitDict['S1']
        s2 = fitDict['S2']
        s3 = fitDict['S3']

        # Trims the Stoke Parameters down to remove any sidebands beyond observedSidebands
        # This does mean that it won't calculate the DOP for any extra sidebands, even if
        # the software detected them.

        while s0[-1,0] > observedSidebands[-1]:
            s0 = s0[:-1,:]
            s1 = s1[:-1,:]
            s2 = s2[:-1,:]
            s3 = s3[:-1,:]

        # This actually calculates the DOP and DOP error
        dop = s0
        dop[1:,1] = np.sqrt((s1[1:,1])**2+(s2[1:,1])**2+(s3[1:,1])**2)/s0[1:,1]
        dop[1:,2] = np.sqrt( (s1[1:,1]**2)*(s1[1:,2]**2)/(s0[1:,1]**2)/((s1[1:,1])**2+(s2[1:,1])**2+(s3[1:,1])**2) + (s2[1:,1]**2)*(s2[1:,2]**2)/(s0[1:,1]**2)/((s1[1:,1])**2+(s2[1:,1])**2+(s3[1:,1])**2) + (s3[1:,1]**2)*(s3[1:,2]**2)/(s0[1:,1]**2)/((s1[1:,1])**2+(s2[1:,1])**2+(s3[1:,1])**2)) + ((s1[1:,1])**2+(s2[1:,1])**2+(s3[1:,1])**2)*s0[1:,2]**2/s0[1:,1]**4

        # Save according to alpha. This will break if you use an alpha not in the usual
        #   0,45,90,-45 range.
        if laserParams['lAlpha'] == 00:
            dop00 = dop[1:,:]
        if laserParams['lAlpha'] == 45:
            dop45 = dop[1:,:]
        if laserParams['lAlpha'] == 90:
            dop90 = dop[1:,:]
        if laserParams['lAlpha'] == -45:
            dopn45 = dop[1:,:]

    # Put into appropriate form for saving
    totdop = [dop00,dop45[:,1:],dop90[:,1:],dopn45[:,1:]]
    totdop = np.hstack(totdop)

    # Saves as txt file with columns, SB Order, 00 DOP and error, 45 DOP and error,
    # 90 DOP and error, -45 DOP and error
    if save_results:
        np.savetxt(saveFileName + "_DOP.txt", totdop, delimiter=',', header=
            'SB Order, 00 DOP, 00 DOP Error, 45 DOP, 45 DOP Error, 90 DOP, 90 DOP Error, -45 DOP, -45 DOP Error')


    # Building the fan compiler causes it to make  2D np arrays with the alpha and
    # gamma angles where each column is a different alpha_NIR excitation
    alphaData, gammaData, _ = fanData.build()

    # save the alpha and gamma results
    if save_results:
        fanData.buildAndSave(saveFileName + "_{}.txt")

    # Now we need to calculate the Jones matrices.
    J = findJ(alphaData,gammaData)

    # Get the T matrix:
    T = makeT(J,crystalAngle)

    # and save the matrices
    if save_results:
        saveT(J, observedSidebands, "{}_JMatrix.txt".format(saveFileName))
        saveT(T, observedSidebands, "{}_TMatrix.txt".format(saveFileName))

    #Now make a Fan Diagram:

    from PyQt5 import QtWidgets
    app = QtWidgets.QApplication([])

    # interpolate so we have alpha and gamma for all excitation pols
    resampledAlpha, resampledGamma = jonesToFans(observedSidebands, J)

    # And create the fan
    f = FanDiagram(resampledAlpha, resampledGamma)
    # And show the fann
    if Plot:
        f.show()
        app.exec_()

    # ToDo: Add code for making the fans look nicer

    # You can save the fan with
    if save_results:
        f.export(saveFileName+"_fanDiagram.png")

    # And to look at the ratios of the T matrix directly:
    if Plot:
        pg.figure("Magnitudes")
        pg.plot(observedSidebands,
                np.abs(T[0,0,:]/T[1,1,:]),
                'o-',
                label="T++/T--")
        pg.plot(observedSidebands,
                np.abs(T[0,1,:]/T[1,0,:]),
                'o-',
                label="T+-/T-+")

        pg.figure("Angles")
        pg.plot(observedSidebands,
                np.angle(T[0,0,:]/T[1,1,:], deg=True),
                'o-',
                label="T++/T--")
        pg.plot(observedSidebands,
                np.angle(T[0,1,:]/T[1,0,:], deg=True),
                'o-',
                label="T+-/T-+")

        pg.show()

    # Ok this is sort of a quick way to get what I want for Origin plotting
    # the relevant T matrix values
    #
    # ToDo: Format the text files to fit with the standards of other txt files

    tmag = np.transpose(np.array([observedSidebands,np.abs(T[0,0,:]/T[1,1,:]),np.abs(T[0,1,:]/T[1,0,:])]))
    tang = np.transpose(np.array([observedSidebands,np.angle(T[0,0,:]/T[1,1,:], deg = True),np.angle(T[0,1,:]/T[1,0,:], deg = True)]))

    if save_results:
        np.savetxt(saveFileName + "_{}.txt".format("TmatrixMag"), tmag, delimiter=',', header='SB Order, |T++/T--|, |T+-/T-+|')
        np.savetxt(saveFileName + "_{}.txt".format("TmatrixAng"), tang, delimiter=',', header='SB Order, Angle(T++/T--), Angle(T+-/T-+)')

    return totdop,alphaData,gammaData,J,T,tmag,tang

def fan_generator(file,observedSidebands,crystalAngle,saveFileName,save_results=False,Plot=True):
    '''
    This is a split off from fan_n_Tmat that just handles the fan. Apparently
    some of the fan creation stops scripts from running so to avoid that we're
    seperating the creation of the fan into this.

    :param file: String of folder name containing 4 polarimetry scans
    :param observedSidebands: np array of observed sidebands. Data will be
        cropped such that these sidebands are included in everything.
    :param crystalAngle: (Float) Angle of the sample from the 010 crystal face
    :saveFileName: Str of what you want to call the text files to be saved
    :save_results: Boolean controls if the fan diagram png is saved.

    Returns:
    :f: This is a fan object. I'm not really sure how to mess with these in
        greater depth.
    '''
    # Make a list of folders for the fan data
    datasets = hsg.natural_glob(file, "*")

    # FanCompiler class keeps alpha and gamma data together for different pols
    fanData = FanCompiler(observedSidebands)

    for data in datasets:
        laserParams, rawData = hsg.hsg_combine_qwp_sweep(data)
        _, fitDict = hsg.proc_n_fit_qwp_data(rawData,vertAnaDir = False, laserParams = laserParams,wantedSbs = observedSidebands)
        # Add the new alpha and gamma sets to the fan class
        fanData.addSet(fitDict)

    alphaData, gammaData, _ = fanData.build()

    from PyQt5 import QtWidgets
    app = QtWidgets.QApplication([])

    # interpolate so we have alpha and gamma for all excitation pols
    resampledAlpha, resampledGamma = jonesToFans(observedSidebands, J)

    # And create the fan
    f = FanDiagram(resampledAlpha, resampledGamma)
    # And show the fann
    if Plot:
        f.show()
        app.exec_()

    # ToDo: Add code for making the fans look nicer

    # You can save the fan with
    if save_results:
        f.export(saveFileName+"_fanDiagram.png")

    return f

def J_T_proc(file,observedSidebands,crystalAngle,saveFileName,save_results=False,Plot=False):
    """
    This function will take a folder of polarimetry scans and produce
    DOP, alpha/gamma, J and T matrices, and Matrix Relations This is the same
    as fan_n_Tmat but doesn't create the fan itself. Otherwise creates pretty
    much everything one would need.

    :param file: String of folder name containing 4 polarimetry scans
    :param observedSidebands: np array of observed sidebands. Data will be
        cropped such that these sidebands are included in everything.
    :param crystalAngle: (Float) Angle of the sample from the 010 crystal face
    :saveFileName: Str of what you want to call the text files to be saved
    :save_results: Boolean controls if things are saved to txt files.
        Currently saves DOP, alpha, gamma, J matrix, T matrix, Fan, and Matrix
        Relations
    :Plot: Boolean controls if plots of matrix relations are displayed

    :return: DOP,alpha,gamma,J matrix, T matrix, Matrix Relations
    """

    # Make a list of folders for the fan data
    datasets = hsg.natural_glob(file, "*")

    # FanCompiler class keeps alpha and gamma data together for different pols
    fanData = FanCompiler(observedSidebands)

    # Initialize arrays for DOP
    dop00 = np.array([])
    dop45 = np.array([])
    dop90 = np.array([])
    dopn45 = np.array([])

    for data in datasets:
        laserParams, rawData = hsg.hsg_combine_qwp_sweep(data)
        try:
            _, fitDict = hsg.proc_n_fit_qwp_data(rawData,vertAnaDir = False,
            laserParams = laserParams,wantedSbs = observedSidebands)
            # Add the new alpha and gamma sets to the fan class
        except IndexError:
            print('incorrect number of files in data folder ',data)
            print('proceeding without fourier analysis')
            _, fitDict = hsg.proc_n_fit_qwp_data(rawData,vertAnaDir = False,
                laserParams = laserParams,
                wantedSbs = observedSidebands,fourier = False)
        # Add the new alpha and gamma sets to the fan class
        fanData.addSet(fitDict)

        # Get Stoke Parameters
        s0 = fitDict['S0']
        s1 = fitDict['S1']
        s2 = fitDict['S2']
        s3 = fitDict['S3']

        # Trims the Stoke Parameters down to remove any sidebands beyond observedSidebands
        # This does mean that it won't calculate the DOP for any extra sidebands, even if
        # the software detected them.

        while s0[-1,0] > observedSidebands[-1]:
            s0 = s0[:-1,:]
            s1 = s1[:-1,:]
            s2 = s2[:-1,:]
            s3 = s3[:-1,:]

        # This actually calculates the DOP and DOP error
        dop = s0
        dop[1:,1] = np.sqrt((s1[1:,1])**2+(s2[1:,1])**2+(s3[1:,1])**2)/s0[1:,1]
        dop[1:,2] = np.sqrt( (s1[1:,1]**2)*(s1[1:,2]**2)/(s0[1:,1]**2)/((s1[1:,1])**2+(s2[1:,1])**2+(s3[1:,1])**2) + (s2[1:,1]**2)*(s2[1:,2]**2)/(s0[1:,1]**2)/((s1[1:,1])**2+(s2[1:,1])**2+(s3[1:,1])**2) + (s3[1:,1]**2)*(s3[1:,2]**2)/(s0[1:,1]**2)/((s1[1:,1])**2+(s2[1:,1])**2+(s3[1:,1])**2)) + ((s1[1:,1])**2+(s2[1:,1])**2+(s3[1:,1])**2)*s0[1:,2]**2/s0[1:,1]**4

        # Save according to alpha. This will break if you use an alpha not in the usual
        #   0,45,90,-45 range.
        if laserParams['lAlpha'] == 00:
            dop00 = dop[1:,:]
        if laserParams['lAlpha'] == 45:
            dop45 = dop[1:,:]
        if laserParams['lAlpha'] == 90:
            dop90 = dop[1:,:]
        if laserParams['lAlpha'] == -45:
            dopn45 = dop[1:,:]

    # Put into appropriate form for saving
    #totdop = [dop00,dop45[:,1:],dop90[:,1:],dopn45[:,1:]]
    #totdop = np.hstack(totdop)

    # Saves as txt file with columns, SB Order, 00 DOP and error, 45 DOP and error,
    # 90 DOP and error, -45 DOP and error
    # if save_results:
    #     np.savetxt(saveFileName + "_DOP.txt", totdop, delimiter=',', header=
    #         'SB Order, 00 DOP, 00 DOP Error, 45 DOP, 45 DOP Error, 90 DOP, 90 DOP Error, -45 DOP, -45 DOP Error')


    # Building the fan compiler causes it to make  2D np arrays with the alpha and
    # gamma angles where each column is a different alpha_NIR excitation
    alphaData, gammaData, _ = fanData.build()

    # save the alpha and gamma results
    if save_results:
        fanData.buildAndSave(saveFileName + "_{}.txt")

    # Now we need to calculate the Jones matrices.
    J = findJ(alphaData,gammaData)

    # Get the T matrix:
    T = makeT(J,crystalAngle)

    # and save the matrices
    if save_results:
        saveT(J, observedSidebands, "{}_JMatrix.txt".format(saveFileName))
        saveT(T, observedSidebands, "{}_TMatrix.txt".format(saveFileName))

    # And to look at the ratios of the T matrix directly:
    if Plot:
        pg.figure("Magnitudes")
        pg.plot(observedSidebands,
                np.abs(T[0,0,:]/T[1,1,:]),
                'o-',
                label="T++/T--")
        pg.plot(observedSidebands,
                np.abs(T[0,1,:]/T[1,0,:]),
                'o-',
                label="T+-/T-+")

        pg.figure("Angles")
        pg.plot(observedSidebands,
                np.angle(T[0,0,:]/T[1,1,:], deg=True),
                'o-',
                label="T++/T--")
        pg.plot(observedSidebands,
                np.angle(T[0,1,:]/T[1,0,:], deg=True),
                'o-',
                label="T+-/T-+")

        pg.show()

    # Ok this is sort of a quick way to get what I want for Origin plotting
    # the relevant T matrix values
    #
    # ToDo: Format the text files to fit with the standards of other txt files

    tmag = np.transpose(np.array([observedSidebands,np.abs(T[0,0,:]/T[1,1,:]),np.abs(T[0,1,:]/T[1,0,:])]))
    tang = np.transpose(np.array([observedSidebands,np.angle(T[0,0,:]/T[1,1,:], deg = True),np.angle(T[0,1,:]/T[1,0,:], deg = True)]))

    if save_results:
        np.savetxt(saveFileName + "_{}.txt".format("TmatrixMag"), tmag, delimiter=',', header='SB Order, |T++/T--|, |T+-/T-+|')
        np.savetxt(saveFileName + "_{}.txt".format("TmatrixAng"), tang, delimiter=',', header='SB Order, Angle(T++/T--), Angle(T+-/T-+)')

    return alphaData,gammaData,J,T,tmag,tang
