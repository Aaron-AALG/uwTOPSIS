# Unweighted TOPSIS method

The Un-Weighted Technique for Order Preference by Similarity to Ideal Solution (uwTOPSIS) ranks decision alternatives based on the classical TOPSIS approach, however this method does not require the introduction of a priori weights.

As a consequence of working with unknown weights, the method does not take into account the relative importance of criteria. Then, the positive ideal solution (PIS) and a negative ideal solution (NIS) varies depending on the conditions of problem. Hence, the function of relative proximity (_R_) is an operator which are optimized as two mathematical programming problems of maximize (_R<sub>L_) and minimize (_R<sub>U_), considering weights as variables. Finally, per each alternative, we get the intervals [_R<sub>L_, _R<sub>U_] so we can rank them in accordance with a determined comparison method.

For a better understanding about either the algorithm or the method, please check:

[V. Liern and B. Pérez-Gladish, “Multiple criteria ranking method based on functional proximity index: un-weighted topsis,” _Annals of Operations Research_, 2020.](https://doi.org/10.1007/s10479-020-03718-1)

## Installation

You can install the uwTOPSIS library from GitHub:

```terminal
git clone https://github.com/Aaron-AALG/uwTOPSIS.git
python3 -m pip install -e uwTOPSIS
```  
  
You can also install it directly from PyPI:
```terminal
pip install uwTOPSIS
```  
  
## Input-Output arguments

**Input**:
>    data: dataframe which contains the alternatives and the criteria.
>
>    directions: array with the optimal direction of the criteria.
>    
>    L: array with the lower bounds of the weigths.
>    
>    U: array with the upper bounds of the weigths.
>    
>    norm: normalization method for the data, whether "euclidean", "gauss", "minmax", "none". (By default norm = "euclidean").
>    
>    p: integer value for the L-p distance. (By default p=2).
>    
>    w0: array with the initial guess of the weights. (By default w0=[]).
>    
>    alpha: value of the convex lineal combination of the uwTOPSIS score. (By default alpha=1/2).
>    
>    forceideal: logical argument to indicate whether to force the ideal solution. If true, ideal solution 
>    is 1-array and antiideal is 0-array. (By default forceideal = False).
>    
>    display: logical argument to indicate whether to show print convergence messages or not. (By default display = False).

**Output**:

Dictionary which contains three keys.
>    Ranking: List with R_min and R_max scores in regard of the optimal weights, plus the uwTOPSIS score.
>
>    Weights_min: List with the weights that minimizes the R score.
>
>    Weights_max: List with the weights that maximizes the R score.
  
## Usage

uwTOPSIS is implemented in order to manage **Pandas** DataFrames as input data which will be converted to **NumPy** arrays. Here is an example based on the paper abovementioned (V. Liern et. al 2020), in which we only use three alternatives and four criteria:

```python
import pandas as pd
import numpy as np
from uwTOPSIS.uwTOPSIS import *

data = pd.DataFrame({"c1":[173, 176, 142],
                     "c2":[10, 11, 5],
                     "c3":[11.4, 12.3, 8.2],
                     "c4":[10.01, 10.48, 7.3]})
directions = ["max", "max", "min", "min"]
L = np.array([0.1 for _ in range(data.shape[1])])
U = np.array([0.4 for _ in range(data.shape[1])])
norm = "euclidean"
p = 2

x = uwTOPSIS(data, directions, L, U, norm, p)
```

The output of the function is a dictionary whose entries are Ranking, Weights_min, and, Weights_max. Besides, Ranking entry is another dictionary with the arguments R_min, R_max, and, uwTOPSIS. The Weights_min and Weights_max output contains the arrays with the optimal solution of each alternative as minimize and maximize respectively.

## Optimization in Python

This library uses the [minimize](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html) function of the scipy.optimize module to carry out the optimization problems. In particular, _R<sub>L_ and _R<sub>U_ are obtained one by one, thus we can compute the gradient and apply the __SLSQP__ method.
