import numpy as np

def finDiff(fun, x0, h):
	return (fun(x0-2*h) - 8*fun(x0-h) + 8*fun(x0 + h) - fun(x0+2*h))/(12*h)

def NR(fun, x0, tol, h):
	if(abs(fun(x0)) < tol):
		return x0
	else:
		x = x0 - fun(x0)/finDiff(fun, x0, h)
		return NR(fun, x, tol,h)

if __name__ == '__main__':
	def fun(x):
		return np.sin(x)

	x0 = 2*np.pi/3
	tol = 1e-10
	h = 0.1
	print(abs(NR(fun, x0, tol, h) - np.pi))