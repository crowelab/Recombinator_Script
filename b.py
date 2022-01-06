import sys
import numpy as np
from scipy.stats import truncnorm
import matplotlib.pyplot as plt

X1 = np.random.normal(5, 1, 10000000)
plt.hist(X1, 10,normed=True)
plt.show()
X=X1.round().astype(int)
plt.hist(X1, 10,normed=True)
plt.show()
#print X
#for i in X:
#   print i
