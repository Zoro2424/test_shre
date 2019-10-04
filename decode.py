#====Python example====
#<div style="overflow:auto;">
#<syntaxhighlight lang="python">
'''
The following Python implementation of Shamir's Secret Sharing is
released into the Public Domain under the terms of CC0 and OWFa:
https://creativecommons.org/publicdomain/zero/1.0/
http://www.openwebfoundation.org/legal/the-owf-1-0-agreements/owfa-1-0

See the bottom few lines for usage. Tested on Python 2 and 3.
'''

from __future__ import division
from __future__ import print_function

import random
import functools

# 12th Mersenne Prime
# (for this application we want a known prime number as close as
# possible to our security level; e.g.  desired security level of 128
# bits -- too large and all the ciphertext is large; too small and
# security is compromised)
_PRIME = 6864797660130609714981900799081393217269435300143305409394463459185543183397656052122559640661454554977296311391480858037121987999716643812574028291115057151
# 13th Mersenne Prime is 2**521 - 1

_RINT = functools.partial(random.SystemRandom().randint, 0)

def _eval_at(poly, x, prime):
    '''evaluates polynomial (coefficient tuple) at x, used to generate a
    shamir pool in make_random_shares below.
    '''
    accum = 0
    for coeff in reversed(poly):
        accum *= x
        accum += coeff
        accum %= prime
    return accum

def make_random_shares(minimum, shares, prime=_PRIME):
    '''
    Generates a random shamir pool, returns the secret and the share
    points.
    '''
    if minimum > shares:
        raise ValueError("pool secret would be irrecoverable")
    poly = [_RINT(prime) for i in range(minimum)]
    points = [(i, _eval_at(poly, i, prime))
              for i in range(1, shares + 1)]
    return poly[0], points

def _extended_gcd(a, b):
    '''
    division in integers modulus p means finding the inverse of the
    denominator modulo p and then multiplying the numerator by this
    inverse (Note: inverse of A is B such that A*B % p == 1) this can
    be computed via extended Euclidean algorithm
    http://en.wikipedia.org/wiki/Modular_multiplicative_inverse#Computation
    '''
    x = 0
    last_x = 1
    y = 1
    last_y = 0
    while b != 0:
        quot = a // b
        a, b = b, a%b
        x, last_x = last_x - quot * x, x
        y, last_y = last_y - quot * y, y
    return last_x, last_y

def _divmod(num, den, p):
    '''compute num / den modulo prime p

    To explain what this means, the return value will be such that
    the following is true: den * _divmod(num, den, p) % p == num
    '''
    inv, _ = _extended_gcd(den, p)
    return num * inv

def _lagrange_interpolate(x, x_s, y_s, p):
    '''
    Find the y-value for the given x, given n (x, y) points;
    k points will define a polynomial of up to kth order
    '''
    k = len(x_s)
    assert k == len(set(x_s)), "points must be distinct"
    def PI(vals):  # upper-case PI -- product of inputs
        accum = 1
        for v in vals:
            accum *= v
        return accum
    nums = []  # avoid inexact division
    dens = []
    for i in range(k):
        others = list(x_s)
        cur = others.pop(i)
        nums.append(PI(x - o for o in others))
        dens.append(PI(cur - o for o in others))
    den = PI(dens)
    num = sum([_divmod(nums[i] * den * y_s[i] % p, dens[i], p)
               for i in range(k)])
    return (_divmod(num, den, p) + p) % p

def recover_secret(shares, 	   prime=6864797660130609714981900799081393217269435300143305409394463459185543183397656052122559640661454554977296311391480858037121987999716643812574028291115057151):
    '''
    Recover the secret from share points
    (x,y points on the polynomial)
    '''
    if len(shares) < 2:
        raise ValueError("need at least two shares")
    x_s, y_s = zip(*shares)
    return _lagrange_interpolate(0, x_s, y_s, prime)

def main():
    '''main function'''
    #secret, shares = make_random_shares(minimum=3, shares=6)
    secret = 3063559254092205807650747490288814162146139100196871271715080212491337958438231709169732944775743992268217969913962936401625573337510431854118472651453903892
    shares = [(168427778, 3063559254092205807650747490288814162146139100196871271715080212491337958438231709169732944775743992268217969913962936401625573337510431854118472651453903892),
	      (168427779, 5409716963002645780047888977554539643130004472674845561514278793595277775840949371047823231801115186288421069156147671806542063559106824754229501660102569117),
	      (168427780, 3470417805014014928539193600953069745295468460519400123413216448282568886399397762865850903490979219985409759029173215286027965072761905524796215584941811165),
	      (168427781, 2878419120972089606715980629617148184331032344463862875022201704319410719806902569554125690095484336015989524481635461552366860114539834115434583609389681147),
	      (168427782, 3156739409999502346852670849870538480859964394153791685926632033488865052561227643625234097937057218444031144425306647656841774667806225384071649237139857054),
	      (168427783, 3267304366804597800944739190404196863137781564086706182918988277609202049669835684409860490582533711880362420353771236115678417951344127536279980566682413626),
	      (168427784, 4763152168408032446855026439464797499484371418768390386478327605584682184404842887866753104718371397762758012541004985057388358773148252674754321881163317183),
	      (168427785, 2428215446241943969397222629087610739096730010989775314009911286517198472065495122766564376117450500866589381408269687918069344721086954895959813118548860280),
	      (168427786, 1346012099245027255137846727095389802665539710218591031924196107275702342801810675138475726268268109034725954792888077009226006685063243040114310221915741361),
	      (168427787, 3167103626579990164481551900631455406227649498339870931463659678762007414973185848385145925751739606156366598497723358368426602431521721329611646179289862914)]

     print('secret:                                                     ',
     secret)
     print('shares:')
      if shares:
   	for share in shares:
          print('  ', share)

     print('secret recovered from minimum subset of shares:             ',
        recover_secret(shares[:10]))
     print('secret recovered from a different minimum subset of shares: ',
        recover_secret(shares[-3:]))

if __name__ == '__main__':
    main()


