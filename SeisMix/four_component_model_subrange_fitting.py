import numpy as np
from scipy import optimize


###### Functions for four-parameter model fitting ######
def inst_model(x, params):
    iw = params[3] * (x-params[0]) + params[1]
    turb = 1./3. * (x-params[0]) + params[1]
    noise = 2. * (x-(params[0]+params[2])) + (params[1] + 1./3. * params[2])
    return np.fmax(iw, np.fmax(turb,noise))

def misfit(params, x, y):
    return np.nansum(np.power((y-inst_model(x,params)),2))



### Fit four-component model with initial parameters par.
### par[0] = iw turb crossover x, par[1] = iw turb crossover y, par[2] = turb length x, par[3] = slope iw. Calls function misfit(par) to perform fitting
def subrange_fitting(logkx, logphi, par):
    best_par = optimize.fmin(misfit, par, args=(logkx, logphi), disp=0)
    return(best_par)

