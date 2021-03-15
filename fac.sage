from fpylll import IntegerMatrix, SVP,FPLLL
import sys
from gmpy2 import gcd,mpz
#FPLLL.set_external_enumerator(None)
ccorn_variant = False

P=[]
def svp(B):
	A = IntegerMatrix.from_matrix(B)
	return SVP.shortest_vector(A)

def first_primes(n):
	p = 1
	P = []
	while len(P) < n:
		p = next_prime(p)
		P += [p]
	return P

def is_smooth(x, P):
	y = abs(x)
	for p in P:
		while p.divides(y):
			y /= p
	return abs(y) == 1


# This piece of code is borrowed from pollards p-1 algorithm
# It checks if a powersmooth number M is gcd((a**M)-1,N) > 1.
def try_factor(N,M,B=10000):
    found = False
    N = mpz(N)
    M = mpz(abs(M))
    for base in range(2,B):
        if gcd(base,N) == 1:
            p = gcd(pow(base,M,N)-1,N)
            if p > 1:
                q = N//p
                print("Base: %d ,M: %d" % (base,M))
                print("Found: %d = %d * %d" % (N,p,q))
                found = True
                break
    return found

# Test if a factoring relation was indeed found.
def test_Schnorr(N, n, prec=1000):
    global P
    if len(P) == 0:
        P = Primes()[:n]
    f = list(range(1, n+1))
    shuffle(f)

    # Scale up and round
    def sr(x):
        return round(x * 2^prec)

    if ccorn_variant: 
        N1 = round((N^(1/(n+1))) * (2^prec)) / (2^prec)
        diag = [sr(N1*f[i]) for i in range(n)] + [sr(N1*ln(N))]
        B = diagonal_matrix(diag, sparse=False)
        for i in range(n):
            B[i, n] = sr(N1*ln(P[i]))
    else:
        diag = [sr(N*f[i]) for i in range(n)] + [sr(N*ln(N))]
        B = diagonal_matrix(diag, sparse=False)
        for i in range(n):
            B[i, n] = sr(N*ln(P[i]))

    b = svp(B)
    if ccorn_variant:
        e = [b[i] / diag[i] for i in range(n)]
        en = (b[n] - sum(e[i]*B[i, n] for i in range(n))) / B[n, n]
        assert en in ZZ
        if en == 1:
            #print("\nFlipping sign of SVP solution with e_{n+1} == %d" % en)
            e = [-ei for ei in e]
            en = -en    # just for consistency
        elif en != -1:
            #print("\nSkipping SVP solution with e_{n+1} == %d" % en)
            return False,False
    else:
        e = [b[i] / sr(N*f[i]) for i in range(n)]


    u = 1
    v = 1
    for i in range(n):
        assert e[i] in ZZ
        if e[i] > 0:
            u *= P[i]^e[i]
        if e[i] < 0:
            v *= P[i]^(-e[i])

    r = (u - v*N, P)
    b  = is_smooth(r[0],r[1])
    f = False
    if b:
        f = try_factor(N,r[0])
    return (b,f)
            

try:
	bits = int(sys.argv[1])
except:
	bits = 400

try:
	n = int(sys.argv[2])
except:
	n = 47

try:
	trials = int(sys.argv[3])
except:
	trials = 100


print("Testing Schnorr's relation finding algorithm with n=%d on RSA-moduli of %d bits, %d trials"%(n, bits, trials))

fc=0
successes = 0
for i in range(trials):
    p = random_prime(2^(bits/2), false, 2^(bits/2-1))
    q = random_prime(2^(bits/2), false, 2^(bits/2-1))
    N = p*q
    success,f = test_Schnorr(N, n)
    successes += success
    if success:
        fc +=f

        print("Trial: %d success: %d" % (i,successes), end="\t")
        print("factored count:",fc)
    sys.stdout.flush()

print("\n%d Factoring Relation found out of %d trials"%(successes, trials))
