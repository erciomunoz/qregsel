import pandas as pd
import io
import requests
import numpy as np
from cvxopt import matrix, solvers
from sfi import Data, Scalar, Matrix


def rotated_qr(varnames, varname, cdf, touse, coefs, nname, kname, dfname):
	X = Data.get(varnames, None, touse)
	X = pd.DataFrame(data=X)
	K = X.shape[1]
	X.insert(K,"intercept",1,True)
	K = X.shape[1]
	N = X.shape[0]
	y = Data.get(varname, None, touse)
	y = pd.DataFrame(data=y)
	tau = Data.get(cdf, None, touse)
	tau = np.array(np.transpose(tau))
	tau_neg = 1-tau
	# equality constraints - left hand side
	A1 = X.to_numpy() # intercepts & data points - positive weights
	A2 = X.to_numpy() * - 1 # intercept & data points - negative weights
	A3 = np.identity(N) # error - positive
	A4 = np.identity(N)*-1 # error - negative
	A = np.concatenate((A1,A2,A3,A4 ), axis= 1) #all the equality constraints
	# equality constraints - right hand side
	b = y.to_numpy()
	#goal function - intercept & data points have 0 weights
	#positive error has tau weight, negative error has 1-tau weight
	c = np.concatenate((np.repeat(0, 2 * K), tau, tau_neg))
	#converting from numpy types to cvxopt matrix
	Am = matrix(A)
	bm = matrix(b)
	cm = matrix(c)
	# all variables must be greater than zero
	# adding inequality constraints - left hand side
	n = Am.size[1]
	G = matrix(0.0, (n,n))
	G[::n+1] = -1.0
	# adding inequality constraints - right hand side (all zeros)
	h = matrix(0.0, (n,1))
	# Solving the model
	sol = solvers.lp(cm,G,h,Am,bm, solver='glpk')
	x = sol['x']
	# both negative and positive components get values above zero, this gets fixed here
	beta = x[0:K] - x[K :2*K]
	Matrix.create(coefs, K, 1, 0)
	Matrix.store(coefs, beta)
	Scalar.setValue(nname, N)
	Scalar.setValue(kname, K)
	Scalar.setValue(dfname, N-K)