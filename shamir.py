'''

Shamir's secret sharing scheme -- Lagrange interpolation over finite fields GF(p)

first we must make the extended euclidean algorithm to find inverses in GF(p)

now we must make the lagrange interpolation formula:
	given a dictionary
		f = {a : f(a), b : f(b), ... }
	we take 
		C = [a, b, ... ] = list(f.keys())
	then define
		delta(/, x, *, i, C) = prod( [(x-j) * inv_mod_p(i-j) for j in C if j != i else 1] )

'''
from functools import reduce

# product operator
prod = lambda X: reduce(lambda x, y: x*y, X)


# ---- NUMBER THEORY & ALGEBRA ---- #

# p defining our field - GF(p)
p = 5575621

# extended euclidean algorithm on hcf(a, p)
def eea_p(a : int, *, p : int = p) -> int:
	if a == 0: return p, 0, 1
	hcf, s, t = eea_p(p % a, p=a)
	return hcf, t - (p // a) * s, s

# multiplicative inverses modulo p
def inv_mod_p(a : int, *, p : int = p) -> int:
	_, s, _ = eea_p(a, p=p)
	return s




# ---- LAGRANGE INTERPOLATION ---- #

# delta function for interpolation
def delta(x, *, i : int, C : list[int], p : int = p): # there are '% p's everywhere because this is done in the finite field GF(p)
	terms = [ (((x-j) % p) * inv_mod_p((i-j) % p, p=p)) % p for j in C if j != i ]
	return prod(terms) % p

def lagrange_interp(x, *, f_dict : dict, p : int = p): # resultant polynomial
	C = list(f_dict.keys())
	return sum([ (f_dict[i]*delta(x, i=i, C=C, p=p)) % p for i in C ]) % p



# ---- RUNTIME ---- #
if __name__ == "__main__":

	# secrets
	F = {870193: -23613404754021249939363940813,
		 485592: -2289717337456309501708473607,
		 3994760: -10487199360175451308104343835783,
		 4325261: -14412723039002678222346852964541,
		 3730509: -7975705554298882208355190391485} # generated using the random_polynomial library with y-intercept 1935737

	# lagrange interpolate
	func = lambda x : lagrange_interp(x, f_dict = F)
	
	# find f(0) -> secret
	print(func(0))
