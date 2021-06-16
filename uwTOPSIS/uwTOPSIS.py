import pandas as pd
from numpy import *
import scipy.optimize as opt
from scipy.spatial.distance import minkowski, chebyshev

def requirements(data, directions, L, U, norm, p, w0, forceideal, display):
    
    # data requirements
    if len(data.shape) < 2:
        raise ValueError('[!] data must be matrix-shaped.')
    elif data.shape[0] < 2 or data.shape[1] < 2:
        raise ValueError('[!] data must be matrix-shaped.')
    
    # directions requirements
    if len(directions) != data.shape[1]:
        raise ValueError('[!] Number of optimal directions must be equal to the number of criteria.')
    if not all([d == 'min' or d == 'max' for d in directions]):
        raise ValueError('[!] Optimal directions must be either "max" or "min".')
    
    # norm requirements
    if norm not in ["euclidean", "gauss", "minmax", "none"]:
        raise ValueError('[!] The normalization method must be either "euclidean", "gauss", "minmax", or "none".')
    
    # p requirements
    if not any([p >= 1, p == "max"]):
        raise ValueError('[!] The p values of the L-p metric must be either a number equal or greater than 1 or "max".')

    # L requirements
    if len(L) != data.shape[1]:
        raise ValueError('[!] Number of lower bounds must be equal to the number of criteria.')
    if any([l < 0 or l > 1 for l in L]):
        raise ValueError('[!] Lower bounds must belong to [0,1].')
    if sum(L) > 1:
        raise ValueError('[!] The sum of lower bounds must be less than 1.')
    
    # U requirements
    if len(U) != data.shape[1]:
        raise ValueError('[!] Number of upper bounds must be equal to the number of criteria.')
    if any([u < 0 or u > 1 for u in U]):
        raise ValueError('[!] Upper bounds must belong to [0,1].')
    if sum(U) < 1:
        raise ValueError('[!] The sum of upper bounds must be greater than 1.')
    
    # w0 requirements
    if len(w0) > 0:
        if len(w0) != data.shape[1]:
            raise ValueError('[!] Length of initial weight must be equal to the number of criteria.')
        if all([w < 0 for w in w0]) and all([w > 1 for w in w0]):
            raise ValueError('[!] Initial weights must belong to [0,1].')
        if abs(sum(w0) - 1) > 10**(-6):
            raise ValueError('[!] Initial weights must sum 1.')
    
    # forceideal requirements
    if int(forceideal) not in [0,1]:
        raise ValueError('[!] "forceideal" must be boolean.')
    
    # display requirements
    if int(display) not in [0,1]:
        raise ValueError('[!] "display" must be boolean.')
    
    return

def normalization(data, norm):
    if norm == 'euclidean':
        data_norm = data.apply(lambda x: x/linalg.norm(x))
    elif norm == 'minmax':
        data_norm = data.apply(lambda x: (x-min(x))/(max(x)-min(x)))
    elif norm == 'gauss':
        data_norm = data.apply(lambda x: 1/(std(x)*sqrt(2*pi))*exp(-1/2*((x-mean(x))/std(x))**2))
    elif norm == 'none':
        data_norm = data
    return data_norm

def get_ideals(data, directions, forceideal, J):
    if forceideal:
        Ideal = [int(i == 'max') for i in directions]
        Antiideal = [int(i == 'min') for i in directions]
    else:
        pos_max = [int(i == 'max') for i in directions]
        pos_min = [int(i == 'min') for i in directions]
        col_max = data.apply(lambda z: max(z))
        col_min = data.apply(lambda z: min(z))
        Ideal = array(col_max)*array(pos_max) + array(col_min)*array(pos_min)
        Antiideal = array(col_max)*(1-array(pos_max)) + array(col_min)*(1-array(pos_min))
    return Ideal, Antiideal

def distance(x, v, p):
    if p == 'max':
        d = chebyshev(x, v)
    else:
        d = minkowski(x, v, p)
    return d

def R(data, Ideal, Antiideal, p):
    # ????? Introduce data[i] instead of the loop ????? 
    d_Ideal = array([distance(x, Ideal, p) for x in data])
    d_Antiideal = array([distance(x, Antiideal, p) for x in data])
    R = d_Antiideal/(d_Ideal + d_Antiideal)
    return R

def R_i(w, data, Ideal, Antiideal, p, i, optimal_mode, J):
    data_norm = w*data
    Ideal = array(w)*Ideal
    Antiideal = array(w)*Antiideal
    r = R(data_norm, Ideal, Antiideal, p)
    if optimal_mode == 'min':
        r_i = r[i]
    else:
        r_i = -r[i]
    return r_i

def R_gradient(w, data, Ideal, Antiideal, p, i, optimal_mode, J):
    d_Ideal = array(distance(data[i], Ideal, p))
    d_Antiideal = array(distance(data[i], Antiideal, p))
    gradient = []
    for k in range(J):
        d_Ideal_partial = w[k]*(Ideal[k]-data[i][k])**2/d_Ideal
        d_Antiideal_partial = w[k]*(Antiideal[k]-data[i][k])**2/d_Antiideal
        gradient.append((d_Antiideal_partial*(d_Ideal+d_Antiideal)-d_Antiideal*(d_Ideal_partial+d_Antiideal_partial))/((d_Ideal+d_Antiideal)**2))
    return gradient

def optimize_TOPSIS(data, Ideal, Antiideal, p, L, U, optimal_mode, display, I, J):
    
    bounds = [(l,u) for l, u in zip(L, U)]
    constraints = ({'type': 'ineq', 'fun': lambda w: 1-sum(w)},
                   {'type': 'ineq', 'fun': lambda w: sum(w)-1},)
    
    # Optimizing the R-score according to the optimal_mode
    r = []
    w = []
    for i in range(I):
        id_max = argmax(abs(data[i]-Antiideal))
        id_min = argmin(abs(data[i]-Antiideal))
        w0 = 1/(J-2)*(1-L-U)
        if optimal_mode == 'max':
            w0[id_max] = U[id_max]
            w0[id_min] = L[id_min]
        if optimal_mode == 'min':
            w0[id_max] = L[id_max]
            w0[id_min] = U[id_min]
        # For Bounded-Constrained problems, we may apply either L-BFGS_B, Powell or TNC methods
        opt_i = opt.minimize(fun = R_i,
                            x0 = w0,
                            jac = R_gradient,
                            args = (data, Ideal, Antiideal, p, i, optimal_mode, J),
                            method = 'SLSQP',
                            bounds = bounds,
                            constraints =  constraints,
                            tol = 10**(-8),
                            options = {'disp': display})
        if optimal_mode == 'max':
            opt_i.fun = -opt_i.fun
        r.append(opt_i.fun)
        w.append(opt_i.x)
    return r, w

def uwTOPSIS(data, directions, L, U, norm = "euclidean", p = 2, w0=[], alpha = 1/2, forceideal = False, display = False):
    """
    uwTOPSIS method
    Input:
        data: dataframe which contains the alternatives and the criteria.
        directions: array with the optimal direction of the criteria.
        L: array with the lower bounds of the weigths.
        U: array with the upper bounds of the weigths.
        norm: normalization method for the data, whether "euclidean", "gauss", "minmax", "none".
        p: integer value for the L-p distance.
        w0: array with the initial guess of the weights.
        alpha: value of the convex lineal combination of the uwTOPSIS score.
        forceideal: logical argument to indicate whether to force the ideal solution. If true, ideal solution is 1-array and antiideal is 0-array.
        display: logical argument to indicate whether to show print convergence messages or not.
    Output:
        Dictionary which contains three keys.
            Ranking: List with R_min and R_max scores in regard of the optimal weights, plus the uwTOPSIS score.
            Weights_min: List with the weights that minimizes the R score.
            Weights_max: List with the weights that maximizes the R score.
    """
    # Check whether the data input verifies the basic requirements
    requirements(data, directions, L, U, norm, p, w0, forceideal, display)
    
    # 1st step: Normalize data
    I = len(data.index)
    J = len(data.columns)
    data_norm = normalization(data, norm)

    # 2nd step: Compute the Ideal, Antiideal elements
    Ideal, Antiideal = get_ideals(data_norm, directions, forceideal, J)
    
    # 3rd step: Optimize R score
    r_min, w_min = optimize_TOPSIS(array(data_norm), Ideal, Antiideal, p, L, U, 'min', display, I, J)
    r_max, w_max = optimize_TOPSIS(array(data_norm), Ideal, Antiideal, p, L, U, 'max', display, I, J)
    uwTOPSIS = [(1-alpha)*m + alpha*M for m, M in zip(r_min,r_max)]
    
    # Output prepation
    scores = {'R_min': r_min, 'R_max': r_max, 'uwTOPSIS': uwTOPSIS}
    output_uwTOPSIS = {'Ranking': scores, 'Weights_min': w_min, 'Weights_max': w_max}

    return output_uwTOPSIS
