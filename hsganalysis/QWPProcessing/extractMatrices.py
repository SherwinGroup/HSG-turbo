import numpy as np
from scipy.optimize import minimize

from hsganalysis.jones import JonesVector as JV
from .expFanCompiler import *


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
            # Note: We current'y don't do error propagation at this step
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
    np.savetxt(out,
               flatT, header=header, comments='', delimiter=',',
               fmt="%.6f")
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


