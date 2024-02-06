import numpy as np
import scipy.optimize as opt
from scipy.spatial.distance import minkowski, chebyshev

def uwTOPSIS_requirements(data, directions, L, U, norm, p, forceideal, display):
    '''
    UW-TOPSIS requirements
    ---------------------

    Checks whether the input parameters satisfy the UW-TOPSIS hypothesis.
    '''
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
    if norm not in ["euclidean", "minmax", "none"]:
        raise ValueError('[!] The normalization method must be either "euclidean", "minmax", or "none".')
    
    # p requirements
    if not any([p >= 1, p == "max"]):
        raise ValueError('[!] The p values of the L-p metric must be either a number equal or greater than 1 or "max".')

    # L requirements
    if len(L) != data.shape[1]:
        raise ValueError('[!] Number of lower bounds must be equal to the number of criteria.')
    if any([l < 0 or l > 1 for l in L]):
        raise ValueError('[!] Lower bounds must belong to [0,1].')
    if np.sum(L) > 1:
        raise ValueError('[!] The sum of lower bounds must be less than 1.')
    
    # U requirements
    if len(U) != data.shape[1]:
        raise ValueError('[!] Number of upper bounds must be equal to the number of criteria.')
    if any([u < 0 or u > 1 for u in U]):
        raise ValueError('[!] Upper bounds must belong to [0,1].')
    if np.sum(U) < 1:
        raise ValueError('[!] The sum of upper bounds must be greater than 1.')
    
    # forceideal requirements
    if int(forceideal) not in [0,1]:
        raise ValueError('[!] "forceideal" must be boolean.')
    
    # display requirements
    if int(display) not in [0,1]:
        raise ValueError('[!] "display" must be boolean.')
    return

def data_normalization(data, norm):
    '''
    Normalization of the decision matrix
    ------------------

    Normalize the decision matrix
    '''
    if norm == 'euclidean':
        data_normalized = np.apply_along_axis(lambda x: x/np.sum(x**2)**0.5,
                                        axis = 0,
                                        arr = data)
    elif norm == 'minmax':
        data_normalized = np.apply_along_axis(lambda x: (x - np.min(x))/(np.max(x) - np.min(x)),
                                        axis = 0,
                                        arr = data)
    elif norm == 'none':
        data_normalized = data
    return data_normalized

def get_ideal_solutions(data, directions, forceideal, J):
    '''
    Extraction of the ideal solutions
    ----------------

    Get the Positive (PIS) and Negative (NIS) ideal solutions.
    '''
    if forceideal:
        # Solutions are binary vectors
        PIS = np.array([i == 'max' for i in directions]).astype("int")
        NIS = np.array([i == 'min' for i in directions]).astype("int")
    else:
        # The ideal solutions are composed of the "best possible" attributes in data
        # Position of max/min directions of criteria
        column_idx_max = np.array([i == 'max' for i in directions]).astype("int")
        column_idx_min = np.array([i == 'min' for i in directions]).astype("int")
        # Values associated to max/min positions
        column_value_max = data.max(axis = 0)
        column_value_min = data.min(axis = 0)
        # PIS and NIS solutions
        PIS = column_value_max * column_idx_max + column_value_min * column_idx_min
        NIS = column_value_max * (1 - column_idx_max) + column_value_min * (1 - column_idx_min)
    return PIS, NIS

def separation_measure(x, v, p):
    '''
    Distance function that defines the relative proximity index
    -----------------------------------------------------------
    '''
    if p == 'max':
        d = chebyshev(x, v)
    else:
        d = minkowski(x, v, p)
    return d

def R_score(data, PIS, NIS, p):
    '''
    Relative proximity index (R-score) of TOPSIS
    -----------------

    Relative proximity index to the ideal solutions
    '''
    # Compute separation measures
    distance_PIS = np.array([separation_measure(x, PIS, p) for x in data])
    distance_NIS = np.array([separation_measure(x, NIS, p) for x in data])
    # Compute the relative proximity index
    R_index = distance_NIS / (distance_PIS + distance_NIS)
    return R_index

def uwTOPSIS_objective_function(weights, data, PIS, NIS, p, i, optimal_mode, J):
    '''
    Objective function of the R-score for the UW-TOPSIS
    -------------

    R-score for the ith alternative in order to optimize it for the UW-TOPSIS method.
    '''
    #
    weights = np.array(weights)
    # Weight the decision matrix and the ideal solutions
    data_weighted = weights * data
    PIS_weighted  = weights * PIS
    NIS_weighted  = weights * NIS
    # Compute the relative proximity index
    r = R_score(data_weighted, PIS_weighted, NIS_weighted, p)
    # Determine the optimality
    if optimal_mode == 'min':
        r_i = r[i]
    else:
        r_i = -r[i]
    return r_i

def create_initial_guess(x, L, U, PIS, NIS, optimal_mode):
    '''
    Create initial weighting scheme  to initialize the optimization
    ---------

    Create a synthetic weighting scheme (w0) for the initialization of `opt.minimize` function
    '''
    # Initialize as lower bounds
    w0 = np.copy(L)
    # Absolute distances regarding PIS and NIS
    abs_dist_PIS = np.abs(x - PIS)
    abs_dist_NIS = np.abs(x - NIS)
    relative_proximity = abs_dist_NIS / (abs_dist_PIS + abs_dist_NIS)
    idx_relative_proximity = np.argsort(relative_proximity)
    # Get max/min indices regarding relative distance
    idx_max_0 = idx_relative_proximity[-1]
    idx_max_1 = idx_relative_proximity[-2]
    idx_min_0 = idx_relative_proximity[0]
    idx_min_1 = idx_relative_proximity[1]
    if optimal_mode == 'max':
        w0[idx_max_0] = U[idx_max_0]
        w0[idx_max_1] = 1 - w0.sum() + w0[idx_max_1]
    if optimal_mode == 'min':
        w0[idx_min_0] = U[idx_min_0]
        w0[idx_min_1] = 1 - w0.sum() + w0[idx_min_1]
    # Normalize the initial guess (just in case)
    w0 /= w0.sum()
    return w0

def optimize_TOPSIS(data, PIS, NIS, p, L, U, optimal_mode, display, I, J):
    '''
    Optimize the R-score of the TOPSIS technique
    --------------------------------------------

    Optimization (min/max) of the R-score per each alternative.
    '''
    # Define bounds and constraints of the optimization problem
    bounds = [(l,u) for l, u in zip(L, U)]
    constraints = ({'type': 'ineq', 'fun': lambda w: 1-sum(w)},
                   {'type': 'ineq', 'fun': lambda w: sum(w)-1},)
    # Optimizing the R-score according to the optimal_mode
    r = []
    w = []
    for i in range(I):
        # Initial guess for the i-th alternative
        w0 = create_initial_guess(data[i], L, U, PIS, NIS, optimal_mode)
        # Optimize the R[i] scores
        opt_i = opt.minimize(fun = uwTOPSIS_objective_function,
                            x0 = w0,
                            args = (data, PIS, NIS, p, i, optimal_mode, J),
                            method = 'SLSQP', # SLSQP, COBYLA, trust-constr
                            bounds = bounds,
                            constraints =  constraints,
                            # tol = 10**(-8),
                            options = {'disp': display},
                            )
        if optimal_mode == 'max':
            opt_i.fun = -opt_i.fun
        r.append(opt_i.fun)
        w.append(opt_i.x)
    return r, w

def uwTOPSIS(data,
             directions,
             L,
             U,
             norm = "euclidean",
             p = 2,
             alpha = 1/2,
             forceideal = False,
             display = False):
    """
    UW-TOPSIS: UnWeighted TOPSIS technique
    ==============================

    Input:
    -------
        data: dataframe which contains the alternatives and the criteria.
        directions: array with the optimal direction of the criteria, whether "max" or "min".
        L: array with the lower bounds of the weigths.
        U: array with the upper bounds of the weigths.
        norm: normalization method for the data, whether "euclidean", "gauss", "minmax", "none".
        p: integer value for the L-p distance.
        alpha: value of the convex lineal combination of the uwTOPSIS score.
        forceideal: logical argument to indicate whether to force the ideal solution. If true, ideal solution (PIS) is a np.ones-array and antiideal is np.zeros-array.
        display: logical argument to indicate whether to show print convergence messages or not.
        
    Output:
    -------
        Dictionary which contains three keys.
            Ranking: List with R_min and R_max scores in regard of the optimal weights, plus the uwTOPSIS score.
            Weights_min: List with the weights that minimizes the R score.
            Weights_max: List with the weights that maximizes the R score.
    """
    # Check whether the data input verifies the basic requirements
    data = np.array(data)
    uwTOPSIS_requirements(data, directions, L, U, norm, p, forceideal, display)
    
    # 1st step: Normalize data
    I, J = data.shape
    data_normalized = data_normalization(data, norm)

    # 2nd step: Compute the PIS and NIS elements
    PIS, NIS = get_ideal_solutions(data_normalized, directions, forceideal, J)
    
    # 3rd step: Optimize the relative proximity index (R-score)
    r_min, w_min = optimize_TOPSIS(data_normalized, PIS, NIS, p, L, U, 'min', display, I, J)
    r_max, w_max = optimize_TOPSIS(data_normalized, PIS, NIS, p, L, U, 'max', display, I, J)
    uwTOPSIS = [(1-alpha) * m + alpha * M for m, M in zip(r_min, r_max)]
    
    # Output of UW-TOPSIS
    scores = {'R_min': r_min, 'R_max': r_max, 'uwTOPSIS': uwTOPSIS}
    output_uwTOPSIS = {'Ranking': scores, 'Weights_min': w_min, 'Weights_max': w_max}

    return output_uwTOPSIS

def TOPSIS(data,
            directions,
            weights,
            norm = "euclidean",
            p = 2,
            forceideal = False):
    """
    TOPSIS technique
    ========================

    Technique for Order of Preference by Similarity to Ideal Solution (TOPSIS) is a
    multi-criteria decision analysis method developed by Hwang and Yoon in 1981.

    Input:
    -------
        data: dataframe which contains the alternatives and the criteria.
        directions: array with the optimal direction of the criteria, whether "max" or "min".
        weights: array with the weighting scheme.
        norm: normalization method for the data, whether "euclidean", "gauss", "minmax", "none".
        p: integer value for the L-p distance.
        forceideal: logical argument to indicate whether to force the ideal solution. If true, ideal solution (PIS) is a np.ones-array and antiideal is np.zeros-array.

    Output:
    -------
        Relative proximity index vector of the TOPSIS.
    """
    # 1st step: Normalize and weight the data
    data = np.array(data)
    I, J = data.shape
    data_normalized = data_normalization(data, norm)
    data_normalized_weighted = weights * data_normalized

    # 2nd step: Compute the PIS and NIS elements
    PIS, NIS = get_ideal_solutions(data_normalized_weighted, directions, forceideal, J)

    # 3rd step: Compute the relative proximity index (R-score)
    topsis = R_score(data_normalized_weighted, PIS, NIS, p)

    return topsis
