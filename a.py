import numpy as np
from scipy.stats import truncnorm
import matplotlib.pyplot as plt


def get_truncated_normal(mean=0, sd=1, low=0, upp=10):
    return truncnorm(
        (low - mean) / sd, (upp - mean) / sd, loc=mean, scale=sd)
#X1 = get_truncated_normal(mean=2, sd=1, low=-10, upp=10)
#X2 = get_truncated_normal(mean=5.5, sd=1, low=-10, upp=10)
X3 = get_truncated_normal(mean=8, sd=1, low=0, upp=10)

import matplotlib.pyplot as plt
fig, ax = plt.subplots(3, sharex=True)
#ax[0].hist(X1.rvs(10000), normed=True)
#ax[1].hist(X2.rvs(10000), normed=True)
ax[2].hist(X3.rvs(10000), normed=True)
plt.show()
X=X3.rvs(10000)
X = X.round().astype(int)
print X
