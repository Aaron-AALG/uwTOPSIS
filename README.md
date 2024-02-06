# Un-Weighted TOPSIS method

![python-version](https://img.shields.io/badge/python->=3.8-green.svg)
![license](https://img.shields.io/pypi/l/uwtopsis.svg)
[![pypi-version](https://img.shields.io/pypi/v/uwtopsis.svg)](https://pypi.python.org/pypi/uwtopsis/)
[![Downloads](https://static.pepy.tech/personalized-badge/uwtopsis?period=total&left_color=grey&right_color=orange&left_text=Downloads)](https://pepy.tech/project/uwtopsis)

The Un-Weighted Technique for Order Preference by Similarity to Ideal Solution (UW-TOPSIS) ranks decision alternatives based on the classical TOPSIS approach, however, this method does not require the introduction of a priori weights. Instead, it makes use of lower and upper bounds to create a weighted decision space that determines the domain of the.

As a consequence of working with unknown weights, the method does not take into account the relative importance of criteria. Then, the positive ideal solution $( PIS )$ and a negative ideal solution $( NIS )$ vary depending on the conditions of the problem. Hence, the function of relative proximity $( R )$ is an operator which is optimized as two mathematical programming problems of maximize $( R^{U} )$ and minimize $( R^{L} )$, considering weights as decision variables. Finally, per each alternative, we get the intervals $[ R^{L}$, $R^{U} ]$ so we can rank them following a determined comparison method.

For a better understanding of either the algorithm or the method, please check:

[V. Liern and B. Pérez-Gladish (2020), “**Multiple criteria ranking method based on functional proximity index: un-weighted TOPSIS**”, _Annals of Operations Research_.](https://doi.org/10.1007/s10479-020-03718-1)

The motivation of this repository is the application of UW-TOPSIS to relatively large datasets as we discussed in the following paper:

[O. Blasco-Blasco, M. Liern-García, A. López-García, S.E. Parada Rico (2021), "**An Academic Performance Indicator Using Flexible Multi-Criteria Methods**";  _Mathematics, Applications of Quantitative Methods in Business and Economics Research_.](https://doi.org/10.3390/math9192396)

# Installation

You can install the uwTOPSIS library from GitHub:

```terminal
git clone https://github.com/Aaron-AALG/uwTOPSIS.git
python3 -m pip install -e uwTOPSIS
```

You can also install it directly from [PyPI](https://pypi.org/project/uwTOPSIS/):

```terminal
pip install uwTOPSIS
```

# Input-Output 

## Input

**data**: dataframe which contains the alternatives and the criteria.

**directions**: array with the optimal direction of the criteria.

**L**: array with the lower bounds of the weights.

**U**: array with the upper bounds of the weights.

**norm**: normalization method for the data, whether "euclidean", "minmax", or "none" (By default norm = "euclidean").

**p**: integer value for the L-p distance (By default p=2).

**alpha**: value of the convex linear combination of the uwTOPSIS score (By default alpha=1/2).

**forceideal**: logical argument to indicate whether to force the ideal solution. If true, the ideal solutions
 are boolean arrays regarding the `directions` (By default forceideal = False).

**display**: logical argument to indicate whether to show print convergence messages or not (By default display = False).

## Output

Dictionary which contains three keys.

**Ranking**: List with R_min and R_max scores in regard to the optimal weights, plus the uwTOPSIS score.

**Weights_min**: List with the weights that minimize the R score.

**Weights_max**: List with the weights that maximize the R score.

## Example

UW-TOPSIS is implemented in order to manage **Pandas** DataFrames as input data which will be converted to **NumPy** arrays. Here is an example based on the paper of V. Liern and B. Pérez-Gladish (2020), in which we only use three alternatives and four criteria:

```python
import pandas as pd
import numpy as np
from uwTOPSIS.uwTOPSIS import *

data = pd.DataFrame({"c1":[173, 176, 142],
                     "c2":[10, 11, 5],
                     "c3":[11.4, 12.3, 8.2],
                     "c4":[10.01, 10.48, 7.3]})
directions = ["max", "max", "min", "min"]
L = np.repeat(0.1, data.shape[1])
U = np.repeat(0.4, data.shape[1])
norm = "euclidean"
p = 2

x = uwTOPSIS(data, directions, L, U, norm, p)
```

The output of the function is a dictionary whose entries are `Ranking`, `Weights_min`, and `Weights_max`. Besides, `Ranking` entry is another dictionary with the arguments `R_min`, `R_max`, and, `uwTOPSIS`. The `Weights_min` and `Weights_max` output contains the arrays with the optimal solution of each alternative as minimize and maximize respectively.

## Generalization to classic TOPSIS

Given that UW-TOPSIS generalizes TOPSIS, we can also compute it by limiting the amplitude of the boundaries. The user can utilize the Numpy numerical epsilon as the difference between lower and upper bounds. Here is an example:

```python
weights = np.array([0.25, 0.2, 0.2, 0.35])
epsilon = np.finfo(float).eps

try:
  x = uwTOPSIS(data,
               directions, 
               weights, 
               weights + epsilon, 
               norm,
               p)
except:
  x = uwTOPSIS(data,
               directions, 
               weights - epsilon, 
               weights, 
               norm,
               p)
```

However, it is strongly recommended to use the TOPSIS function included in our package instead:

```python
x = TOPSIS(data, directions, weights, norm, p)
```


## Optimization in Python

This library uses the [minimize](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html) function of the `scipy.optimize` module to carry out the optimization problems. In particular, $R^{L}$ and $R^{U}$ are obtained one at a time. Thus, we can apply the SLSQP optimization method.


# Literature review of UW-TOPSIS

Since the first implementation of UW-TOPSIS in MCDA in 2020, several researchers in the field have shown interest in this technique. The following table shows the works in which UW-TOPSIS has been used as a method for the case study or experimental part.

| Year | Title    |
|------|----------|
| 2020 | [V. Liern and B. Pérez-Gladish, Multiple criteria ranking method based on functional proximity index: un-weighted TOPSIS, Annals of Operations Research](https://doi.org/10.1007/s10479-020-03718-1) |
| 2021 | [R. Benítez and V. Liern, Unweighted TOPSIS: a new multi-criteria tool for sustainability analysis, International Journal of Sustainable Development & World Ecology](https://doi.org/10.1080/13504509.2020.1778583) |
| 2021 | [V. Liern, B. Pérez-Gladish, F. Rubiera-Morollón and B. M'Zali, Residential choice from a multiple criteria sustainable perspective, Annals of Operations Research](https://doi.org/10.1007/s10479-021-04480-8) |
| 2021 | [V. Liern and B. Pérez-Gladish, Building Composite Indicators With Unweighted-TOPSIS, IEEE Transactions on Engineering Management](https://doi.org/10.1109/TEM.2021.3090155) |
| 2021 | [B. Pérez-Gladish and F. A. Ferreira, MCDM/A studies for economic development, social cohesion and environmental sustainability: introduction](https://doi.org/10.1080/13504509.2020.1821257) |
| 2021 | [O. Blasco-Blasco, M. Liern-García, A. López-García, S.E. Parada Rico, An Academic Performance Indicator Using Flexible Multi-Criteria Methods,  Mathematics, Applications of Quantitative Methods in Business and Economics Research](https://doi.org/10.3390/math9192396)|
| 2021 | [J. Vicens-Colom, J. Holles and V. Liern, Measuring Sustainability with Unweighted TOPSIS: An Application to Sustainable Tourism in Spain, Sustainability](https://doi.org/10.3390/math9192396)|
| 2022 | [V. Liern and B. Pérez-Gladish, Multiple criteria ranking method based on functional proximity index: Un-weighted TOPSIS, Annals of Operations Research](https://doi.org/10.1007/s10479-020-03718-1) |
| 2022 | [T. Fernández-García, V. Liern, B. Pérez-Gladish and F. Rubiera-Morollón, Measuring the territorial effort in research, development, and innovation from a multiple criteria approach: Application to the Spanish regions case, Technology in Society](https://doi.org/10.1016/j.techsoc.2022.101975) |
| 2022 | [K. Bouslah, V. Liern, J. Ouenniche and B. Pérez-Gladish. Ranking firms based on their financial and diversity performance using multiple-stage unweighted TOPSIS, International Transactions in Operational Research](https://doi.org/10.1111/itor.13143) |
| 2023 | [V. Liern and B. Pérez-Gladish. Measuring Corporate Gender Diversity and Inclusion with UW-TOPSIS and Linguistic Intervals, Operational Research Methods in Business, Finance and Economics](https://doi.org/10.1007/978-3-031-31241-0_4) |
| 2023 | [A. López-García, V. Liern and B. Pérez-Gladish. Determining the underlying role of corporate sustainability criteria in a ranking problem using UW-TOPSIS, Annals of Operations Research](https://doi.org/10.1007/s10479-023-05543-8)|
| 2023 | [A. López-García, O. Blasco-Blasco, M. Liern-García and S.E. Parada Rico. Early detection of students’ failure using Machine Learning techniques, Operations Research Perspectives](https://doi.org/10.1016/j.orp.2023.100292)|
| 2023 | [A. López-García, Evaluation of optimal solutions in multicriteria models for intelligent decision support, Ph.D. Thesis](https://doi.org/10.13140/RG.2.2.19136.30724)|
