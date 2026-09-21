# Numerical functions copied unchanged from BioXTAS RAW SASCalc.py.
# Copyright BioXTAS RAW contributors. GPL-3.0-or-later; see LICENSE_RAW.txt.
# Full upstream source is in ../bioxtasraw-source on the Desktop installation.
import numpy as np
import packaging.version as packv
pnp_version = packv.Version(np.__version__)

def runEFA(A, forward=True):
    """Runs the forward or backward evolving factor calculations."""
    slist = np.zeros_like(A)

    jmax = A.shape[1]

    if not forward:
        A = A[:,::-1]

    for j in range(jmax):
        s = np.linalg.svd(A[:, :j+1], full_matrices = False, compute_uv = False)
        slist[:s.size, j] = s

        if j > 0:
            slist[s.size-1:, j-1] = s[-1]

    return slist

def runRotation(D, intensity, err, ranges, force_positive, svd_v, previous_results=None,
    method='Hybrid', niter=1000, tol=1e-12):
    """
    Runs the full EFA rotation.

    :param numpy.array intensity: An array of scattering intensities with intensities
        corresponding to a single experiment in each column and intensities
        corresponding to a single q in each row.

    :param numpy.array err: An array of scattering intensity errors with each
        column being the errors associated with the corresponding column in the
        intensity matrix.

    :param numpy.array ranges: An array of the EFA ranges, which each row is
        the range of a given component, and contains two integers, the first is
        start of that component the second the end of the component.

    :param list force_positive: A list with as many entries as there are EFA
        components (i.e. rows in the ranges array). Each entry is either True or False,
        corresponding to whether the intensity is forced to remain positive or not.

    :param np.array svd_v: The right singular matrix obtained from SVD of the
        appropriate intensity matrix (either error normalized or otherwise). Note that
        for numpy.svd results, it initial returns V transpose. So if you had
        U, S, Vt = np.linalg.svd(...), this V = Vt.T.

    :param list previous_results: This list contains the previous results of
        EFA, if any. This is used to improve the starting guess. The first
        entry in the list is a boolean corresponding to whether the previous
        results converged, and the second is a dictionary of the previous results,
        corresponding to the rotation_data dictionary returned by this function.
        Defaults to None, which should be used if no previous results are available.

    :param str method: The method of EFA to be used. Options are 'Hybrid', 'Iterative',
        and 'Explicit'. Defaults to 'Hybrid'.

    :param int niter: Number of iterations to run in the iterative and hybrid methods.
        Defaults to 1000.

    :param float tol: Tolerance for convergence of the iterative and hybrid methods.
        Defaults to 1e-12.
    """

    init_dict = {'Hybrid'       : initHybridEFA,
                'Iterative'     : initIterativeEFA,
                'Explicit'      : initExplicitEFA}

    run_dict = {'Hybrid'        : runIterativeEFARotation,
                'Iterative'     : runIterativeEFARotation,
                'Explicit'      : runExplicitEFARotation}

    #Calculate the initial matrices
    num_sv = ranges.shape[0]

    M = np.zeros_like(svd_v[:,:num_sv])

    if previous_results is not None:
        converged, rotation_data = previous_results
    else:
        converged = False
        rotation_data = {}

    for j in range(num_sv):
        M[ranges[j][0]:ranges[j][1]+1, j] = 1

    if converged and M.shape[0] != rotation_data['C'].shape[0]:
        converged = False
        rotation_data = {}

    if converged:
        C_init = rotation_data['C']
    else:
        C_init = svd_v[:,:num_sv]

    V_bar = svd_v[:,:num_sv]

    failed, C, T = init_dict[method](M, num_sv, D, C_init, converged, V_bar) #Init takes M, num_sv, and D, C_init, and converged and returns failed, C, and T in that order. If a method doesn't use a particular variable, then it should return None for that result
    C, failed, converged, dc, k = run_dict[method](M, D, failed, C, V_bar, T, niter, tol, force_positive) #Takes M, D, failed, C, V_bar, T in that order. If a method doesn't use a particular variable, then it should be passed None for that variable.

    if not failed:
        if method != 'Explicit':
            conv_data = {'steps'   : dc,
                'iterations': k,
                'final_step': dc[-1],
                'options'   : {'niter': niter, 'tol': tol, 'method': method},
                'failed'    : False,
                }
        else:
            conv_data = {'steps'   : None,
                'iterations': None,
                'final_step': None,
                'options'   : {'niter': niter, 'tol': tol, 'method': method},
                'failed'    : False,
                }

    else:
        conv_data = {'steps'   : None,
                'iterations': None,
                'final_step': None,
                'options'   : {'niter': niter, 'tol': tol, 'method': method},
                'failed'    : True,
                }

    #Check whether the calculation converged

    if method != 'Explicit':
        if k == niter and dc[-1] > tol:
            converged = False
        elif failed:
            converged = False
        else:
            converged = True

    else:
        if failed:
            converged = False
        else:
            converged = True

    if converged:
        #Calculate SAXS basis vectors
        mult = np.linalg.pinv(np.transpose(M*C))
        intensity = np.dot(intensity, mult)
        err = np.sqrt(np.dot(np.square(err), np.square(mult)))

        int_norm = np.dot(D, mult)
        resid = D - np.dot(int_norm, np.transpose(M*C))

        chisq = np.mean(np.square(resid),0)

        #Save the results
        rotation_data = {'M'   : M,
            'C'     : C,
            'int'   : intensity,
            'err'   : err,
            'chisq' : chisq}

    return converged, conv_data, rotation_data

def runExplicitEFARotation(M, D, failed, C, V_bar, T, niter, tol, force_pos):
    num_sv = M.shape[1]

    for i in range(num_sv):
        V_i_0 = V_bar[np.logical_not(M[:,i]),:]

        T[i,1:num_sv] = -np.dot(V_i_0[:,0].T, np.linalg.pinv(V_i_0[:,1:num_sv].T))

    C = np.dot(T, V_bar.T)

    C = C.T

    if -1*C.min() > C.max():
        C = C*-1

    converged = True

    csum = np.sum(M*C, axis = 0)
    if pnp_version >= packv.Version('1.10'):
        C = C/np.broadcast_to(csum, C.shape) #normalizes by the sum of each column
    else:
        norm = np.array([csum for i in range(C.shape[0])])

        C = C/norm #normalizes by the sum of each column

    return C, failed, converged, None, None

def runIterativeEFARotation(M, D, failed, C, V_bar, T, niter, tol, force_pos):
    #Carry out the calculation to convergence
    k = 0
    converged = False

    dc = []

    while k < niter and not converged and not failed:
        k = k+1
        try:
            Cnew = EFAUpdateRotation(M, C, D, force_pos)
        except np.linalg.linalg.LinAlgError:
           failed = True

        dck = np.sum(np.abs(Cnew - C))

        dc.append(dck)

        C = Cnew

        if dck < tol:
            converged = True

    return C, failed, converged, dc, k

def EFAUpdateRotation(M,C,D, force_pos):
    S = np.dot(D, np.linalg.pinv(np.transpose(M*C)))

    Cnew = np.transpose(np.dot(np.linalg.pinv(S), D))

    for i, fp in enumerate(force_pos):
        if fp:
            Cnew[Cnew[:,i] < 0,i] = 0

    csum = np.sum(M*Cnew, axis = 0)

    if pnp_version >= packv.Version('1.10'):
        Cnew = Cnew/np.broadcast_to(csum, Cnew.shape) #normalizes by the sum of each column
    else:
        norm = np.array([csum for i in range(Cnew.shape[0])])

        Cnew = Cnew/norm #normalizes by the sum of each column

    return Cnew

def EFAFirstRotation(M,C,D):
    #Have to run an initial rotation without forcing C>=0 or things typically fail to converge (usually the SVD fails)
    S = np.dot(D, np.linalg.pinv(np.transpose(M*C)))

    Cnew = np.transpose(np.dot(np.linalg.pinv(S), D))

    csum = np.sum(M*Cnew, axis = 0)
    if pnp_version >= packv.Version('1.10'):
        Cnew = Cnew/np.broadcast_to(csum, Cnew.shape) #normalizes by the sum of each column
    else:
        norm = np.array([csum for i in range(Cnew.shape[0])])

        Cnew = Cnew/norm #normalizes by the sum of each column

    return Cnew

def initIterativeEFA(M, num_sv, D, C, converged, V_bar):

    #Set a variable to test whether the rotation fails for a numerical reason
    failed = False

    #Do an initial rotation
    try:
        C = EFAFirstRotation(M, C, D)
    except np.linalg.linalg.LinAlgError:
        failed = True

    return failed, C, None

def initExplicitEFA(M, num_sv, D, C, converged, V_bar):

    T = np.ones((num_sv, num_sv))

    failed = False

    return failed, None, T

def initHybridEFA(M, num_sv, D, C, converged, V_bar):
    failed = False

    if not converged:
        failed, temp, T = initExplicitEFA(M, num_sv, D, C, converged, V_bar)
        C, failed, temp1, temp2, temp3 = runExplicitEFARotation(M, None, None, None, V_bar, T, None, None, None)

    return failed, C, None
