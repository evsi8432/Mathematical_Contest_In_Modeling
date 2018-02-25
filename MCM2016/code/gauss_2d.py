def gauss_2d(x0, y0):
	mean = [x0, y0]
	cov = [[1, 0], [0, 1]]

	x,y = np.random.multivariate_normal(mean, cov, 500).T
	return zip(x,y)


# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

def was_infected(lam,dt):
    import random
    from scipy import exp
    rnd = random.random()
    if rnd <= 1-exp(-dt*lam):
        return True
    else:
        return False