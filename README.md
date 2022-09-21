# shamir.py
A Python 3.10 implementation of Shamir's Secret Sharing (S3) scheme.

# How the S3 scheme works
The scheme works based off of Galois fields of prime order, we will denote by $\mathrm{GF}(p)$ the Galois field of prime order $p$. We omit the square brackets denoting equivalence classes modulo $p$ to aid with readability.

The basic idea behind the scheme is that any polynomial of degree $n$ is uniquely defined by $(n+1)$ points. If we wish to encode some message $m$ we can simply pick some $n$-th polynomial that somehow encodes $m$ and then share $(n+1)$ distinct points on this curve, say to $(n+1)$ trustworthy friends. To retrieve our secret, a recipient must then find the curve passing through our points _et voilà_!


## Encoding $m$ in a polynomial
There are many ways for a curve to encode a message - e.g. the function's value at a certain point, the curvature at a certain point etc. - but by far the simplest is to take its value at $x=0$. This is simple both graphically and algebraically, as this is just its leading constant term.

Hence the curve we initially work with to define our points will be of the form 
$$f(x) = m + a_1 x + \cdots + a_n x^n.$$ 

Since we are working over $\mathrm{GF}(p)$ we pick $a_1, \cdots, a_n \in \mathrm{GF}(p)$ to define our function $f(x)$, and take values of $f$ modulo $p$.


## Intermediate step
We now pick $(n+1)$ values $x_1, \cdots, x_{n+1} \in \mathrm{GF}(p)$, producing our points $P_i := (x_i, f(x_i))$. We then give out these points to our trustworthy friends.


## Retrieving $m$ from $\\{P_i\\}$
To reconstruct the curve from these points we use a version of [Lagrange interpolation](https://en.wikipedia.org/wiki/Lagrange_polynomial) adjusted for Galois fields.

For each point $P_i$ we define a polynomial $\delta_i$ using the Lagrange interpolation formula:
$$\delta_i(x) := \prod_{1 \leq j \leq n+1 \atop j\neq i} (x-x_j)(x_i-x_j)^{-1},$$
where $x^{-1}$ represents the multiplicative inverse of $x$ in $\mathrm{GF}(p)$. We then define the reconstructed polynomial $\hat{f}(x)$ by the following formula:
$$\hat{f}(x) = \sum_{i=1}^{n+1}f(x_i)\delta_i(x).$$


# Implementation in Python 3.10
We have 2 main parts of our implementation:
- Multiplicative inverses in $\mathrm{GF}(p)$, and
- Lagrange interpolation in $\mathrm{GF}(p)$.

## Multiplicative inverses in $\mathrm{GF}(p)$
To find the multiplicative inverse of the element $a \in \mathrm{GF}(p)$ we us the Extended Euclidean Algorithm (EEA) for $\mathrm{hcf}(a, p)$. This returns three integers, namely the highest common factor of $a$ and $p$, and two integers $s, t \in \mathbb{Z}$ such that $as +pt = \mathrm{hcf}(a, p)$.

Since $p$ is prime we have that by definition $\mathrm{hcf}(a,p)=1$, and so $s, t$ satisfy $as+pt=1$. Taking this equation modulo $p$ we find that 
$$as \equiv 1 \quad (\mathrm{mod}\\,p) \implies s \equiv a^{-1}\quad (\mathrm{mod}\\,p).$$

We implement the EEA in Python as follows:
```python
def eea_p(a : int, *, p : int = p) -> int:
	if a == 0: return p, 0, 1
	hcf, s, t = eea_p(p % a, p=a)
	return hcf, t - (p // a) * s, s
```

## Lagrange interpolation in $\mathrm{GF}(p)$
To undergo Lagrange interpolation in $\mathrm{GF}(p)$ we use the following code to find $\delta_i (x)$:
```python
# delta function for interpolation
def delta(x, *, i : int, C : list[int], p : int = p): 
	# there are '% p's everywhere because this is done in the finite field GF(p)
	terms = [ (((x-j) % p) * inv_mod_p((i-j) % p, p=p)) % p for j in C if j != i ]
	return prod(terms) % p
```

We then apply the last step of interpolation to yield $\hat{f}$:
```python
def lagrange_interp(x, *, f_dict : dict, p : int = p): 
	# resultant polynomial from interpolation
	C = list(f_dict.keys())
	return sum([ (f_dict[i]*delta(x, i=i, C=C, p=p)) % p for i in C ]) % p
```
