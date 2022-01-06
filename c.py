import matplotlib.pyplot as plt
import scipy.stats as stats

#lower, upper = 3.5, 6
#mu, sigma = 5, 0.7
#lower, upper = 0, 5
mu, sigma = 2, 1
lower, upper = 0, mu+3*sigma

X = stats.truncnorm(
    (lower - mu) / sigma, (upper - mu) / sigma, loc=mu, scale=sigma)
N = stats.norm(loc=mu, scale=sigma)

fig, ax = plt.subplots(2, sharex=True)
ax[0].hist(X.rvs(10000), normed=True)
ax[1].hist(N.rvs(10000), normed=True)
plt.show()
