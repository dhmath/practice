"""
Created on Tue Jul 08 12:17:31 2014

@author: DH
"""

###################################################################
#   A probabilistic primality test
#
#   First checks against a list of primes (default: primes < 100)
#   Then runs Rabin-Miller
#   Returns:
#             0 if composite, 
#             1 if definitely prime,
#             2 if probably prime
###################################################################
    
def modSqr(a, n):
    ''' modular squaring '''
    return (a ** 2) % n

def binary(n):
    ''' get n in binary '''
    a = n
    b = a % 2
    bits = []
    while a > 0:
        #print a, b
        bits.insert(0, b)
        a = (a - b) / 2
        b = a % 2

    return bits

def modExp(a, e, n):
    ''' get a ** e mod n by repeated squaring '''
    bits = binary(e)
    
    pow = 1
    for b in bits:
        pow = ((pow ** 2) * (a ** b)) % n
    
    return pow

def sieveOfErat(N):
    ''' sieve of Eratosthenes (get list of primes <= N) '''
    primes = []
    proper_divisors = [0 for i in range(N + 1)]
    
    # Loop through proper_divisors and if nth element is 0, 
    # add n to list of primes and increment every nth element of proper_divisors. 
    # Since we don't consider 0 or 1 as primes, we'll start loop at 2.
    # (This could also be adapted to calculate the divisor function or sum-of-divisor-powers)
    
    for n in range(2, N + 1):    
        if proper_divisors[n] == 0:
            primes.append(n)
            for i in range(2 * n, N + 1, n):
                proper_divisors[i] += 1
                
    return primes


# when we build a list of primes less than ub (upper bound),
# I want to keep track of the value of ub that we used
class PrimeListSieved:
    def __init__(self, ub, show_output = 1):
        self.ub = ub
        self.primes = sieveOfErat(ub)
        self.length = len(self.primes)
        if show_output:
            print("{0} primes up to {1}".format(self.length, self.ub))
            
            
def isPrimeCheckAgainstList(n, ub = 100, show_output = 0):
    ''' check n against list of primes
    return 0 if composite
    return 1 if prime
    return 2 if inconclusive (i.e., for all p in list, n%p != 0 and n > p ** 2)
    '''
    
    # make sure prime list has been defined
    global prime_list_global
    
    try:
        prime_list_global
    except NameError:
        if show_output:
            print("prime_list_global not yet defined. Creating it now:\n")
        prime_list_global = PrimeListSieved(ub, show_output)   
    
    # if necessary, redo sieve to make longer list of primes
    if prime_list_global.ub < ub:
        if show_output:
            print("Length of prime_list_global is less than {0}. Rerunning sieve.\n".format(ub))
        prime_list_global = PrimeListSieved(ub, show_output)
    
    for p in prime_list_global.primes:
        if p == n:
            return 1
        elif n % p == 0:
            return 0
        
    # if sqrt(n) < largest prime in list and 
    # n passed through previous step, then n is prime
    # otherwise test is inconclusive
    
    if prime_list_global.ub ** 2 < n:
        return 2
    else:
        return 1
    
    
def isPrimeRabinMiller(n, p):
    ''' test if n is a base-p pseudoprime
    returns 0 if composite, 2 if it passes test
    returns -1 for n = 1
    '''
    test_pass = 2
    if n == 1:
        return(-1)
    # treat n as composite until it passes test
    test = 0
    
    # find largest power of 2 dividing n - 1,
    # and put n - 1 = (2 ^ s) * m
    s = 0
    m = n - 1
    while m % 2 == 0:
        s += 1
        m = m / 2
    
    # look at the sequence {p ^ (m * (2 ^ j)) mod n}
    # if it's 1 for j = 0 , or if some other term is n - 1,
    # then n passes the test
    power = modExp(p, m, n)
    if power == 1:
        # n passes test if p ^ m == 1 mod n
        return(test_pass)
    for j in range(s):
        if power == n - 1:
            # n passes test if some p ^ (m * (2 ^ j)) == -1 mod n
            return(test_pass)
        power = modSqr(power, n)
    return(test)


def isPrime(n, ub = 100, prime_bases = 9, show_output = 0):
    ''' check n against list of primes <= ub,
    then run Miller-Rabin on prime bases
    '''
    
    if n == 1:
        return(-1)
    
    # make sure prime list has been defined
    global prime_list_global
    
    try:
        prime_list_global
    except NameError:
        if show_output:
            print("prime_list_global not yet defined. Creating it now:\n")
        prime_list_global = PrimeListSieved(ub, show_output)   
    
    # if necessary, redo sieve to make longer list of primes
    if prime_list_global.ub < ub:
        if show_output:
            print("Length of prime_list_global is less than {0}. Rerunning sieve.\n".format(ub))
        prime_list_global = PrimeListSieved(ub, show_output)
    
    number_of_prime_bases = min(prime_list_global.length, prime_bases)
    
    # run check against list of primes
    test = isPrimeCheckAgainstList(n, ub)
    
    # Run Rabin-Miller on prime base.
    # The test has been verified for n < 3,825,123,056,546,413,051
    # when prime_bases <= 9
    # (see http://en.wikipedia.org/wiki/Miller%E2%80%93Rabin_primality_test)
    j = 0
    while test == 2 and j < number_of_prime_bases:
        p = prime_list_global.primes[j]
        test = isPrimeRabinMiller(n, p)
        j += 1
    return(test)


def primePi(X):
    ''' # of (definite and probable) primes less than X '''
    count = 0
    n = 1         # we don't want to count 1 as prime
    while n < X:
        n += 1
        if isPrime(n):
            count += 1
    return(count)
    
    

######################################################################
###########                SOME EXAMPLES           ###################
######################################################################

# this should output 9592
X = 100000
print "primePi({0}) = {1}\n".format(X, primePi(X))

##################       some Mersenne numbers      ##################

# M_607 is prime
N = 607
print "M_{0} : {1}".format(N, isPrime(2 ** N - 1))

# M_613 is not prime
N = 613
print "M_{0} : {1}".format(N, isPrime(2 ** N - 1))

# M_1279 is prime
N = 1279
print "M_{0} : {1}".format(N, isPrime(2 ** N - 1))

# not prime
N = 1283
print "M_{0} : {1}".format(N, isPrime(2 ** N - 1))

# prime
N = 4253
print "M_{0} : {1}".format(N, isPrime(2 ** N - 1))

print 
################    some more large composite numbers    ############

# RSA-100
m = 100
n = 1522605027922533360535618378132637429718068114961380688657908494580122963258952897654000350692006139
print "RSA-{0} : {1}".format(m, isPrime(n))

# RSA-2048
m = 2048
n = 25195908475657893494027183240048398571429282126204032027777137836043662020707595556264018525880784406918290641249515082189298559149176184502808489120072844992687392807287776735971418347270261896375014971824691165077613379859095700097330459748808428401797429100642458691817195118746121515172654632282216869987549182422433637259085141865462043576798423387184774447920739934236584823824281198163815010674810451660377306056201619676256133844143603833904414952634432190114657544454178424020924616515723350778707749817125772467962926386356373289912154831438167899885040445364023527381951378636564391212010397122822120720357
print "RSA-{0} : {1}".format(m, isPrime(n))
