import math
import itertools
import time
import fractions
from itertools import count

def quadraticRoots(a, b, c):
    return([(-b + math.sqrt(b ** 2 - 4 * a * c)) / (2 * a), (-b - math.sqrt(b ** 2 - 4 * a * c)) / (2 * a)])

def intLen(n):
    return(len(str(n)))

def digitsSimulation():
    total1 = 0
    total2 = 0
    total3 = 0
    results = []
    for n in range(1, 100000):
        for i in range(10):
            start1 = time.clock()
            len(str(n))
            end1 = time.clock()
            start2 = time.clock()
            math.floor(math.log10(n)) + 1
            end2 = time.clock()
            start3 = time.clock()
            num = n
            digits = 0
            while num > 0:
                num //= 10
                digits += 1
            digits
            end3 = time.clock()
            total1 += (end1 - start1)
            total2 += (end2 - start2)
            total3 += (end3 - start3)
            results.append([total1 / 1000000, total2 / 1000000, total3 / 1000000])
            total1 = 0
            total2 = 0
            total3 = 0
    print(results)

izip = itertools.zip_longest
chain = itertools.chain.from_iterable
compress = itertools.compress
# Input n >= 6
# Returns primes 2 <= p < n
def generatePrimes(n):
    zero = bytearray([False])
    size = n // 3 + (n % 6 == 2)
    sieve = bytearray([True]) * size
    sieve[0] = False
    for i in range(int(n ** 0.5) // 3 + 1):
        if sieve[i]:
            k = 3 * i + 1 | 1
            start = (k * k + 4 * k - 2 * k * (i & 1)) // 3
            sieve[(k * k) // 3 :: 2 * k] = zero * ((size - (k * k) // 3 - 1) // (2 * k) + 1)
            sieve[start :: 2 * k] = zero * ((size - start - 1) // (2 * k) + 1)
    ans = [2, 3]
    poss = chain(izip(*[range(i, n, 6) for i in (1, 5)]))
    ans.extend(compress(poss, sieve))
    return ans

primes = generatePrimes(7654322)

# Returns index of nth instance of x in arr
def nthIndex(n, x, arr):
    i = -1
    for j in range(n):
        i = arr.index(x, i + 1)
    return(i)
        
def isPalindrome(n):
    return(str(n) == str(n)[::-1])

def isPrime(n):
    num = int(n)
    if num <= 1:
        return False
    for i in range(2, math.floor(math.sqrt(num)) + 1):
        if num % i == 0:
            return False
    return True

# Multiples of 3 and 5
def euler1(n = 1000):
    total = 0
    for i in range(n):
        if i % 3 == 0 or i % 5 == 0:
            total += i
    print("Sum: " + str(total))

# Even Fibonacci numbers
def euler2(n = 4000000):
    total = 0
    back1 = 1
    back2 = 1
    while back1 + back2 < n:
        i = back1 + back2
        if i % 2 == 0:
            total += i
        back2 = back1
        back1 = i
    print("Sum: " + str(total))

# Largest prime factor   
def euler3(n = 600851475143):
    factor = 2
    candidate = n
    while factor < math.floor(math.sqrt(candidate)):
        if candidate % factor == 0:
            candidate /= factor
        else:
            factor += 1
    print("Largest prime factor: " + str(candidate))
         
# Largest palindrome product
def euler4(lb, ub):
    maxPalindrome = -1
    for x in range(lb, ub):
        for y in range(lb, x):
            if isPalindrome(x * y) and x * y > maxPalindrome:
                maxPalindrome = x * y
    print("Largest palindrome: " + str(maxPalindrome))

# Highly divisible triangular number
# Brute force. Assume i > 1
def euler12(minFactors):
    n = 1
    testFactor = 2
    factorCounter = 2
    while True:
        i = n * (1 + n) / 2
        while testFactor <= math.floor(math.sqrt(i)):
            if i % testFactor == 0:
                factorCounter += 2
            testFactor += 1
        if factorCounter >= minFactors:
            print(i)
            break
        n += 1
        factorCounter = 2

# 10001st prime
def euler7(n, ub):
    candidates = [True for i in range(ub)]
    candidates[0] = False
    for j in range(2, math.floor(math.sqrt(ub)) + 1):
        for k in range(2, ub // j + 1):
            candidates[j * k - 1] = False
    l = -1
    for m in range(n):
        l = candidates.index(True, l + 1)
    print(l + 1)

# Largest product in a series
def euler8(n):
    openFile = open("C:\\Users\\RR\\Dropbox\\CS\\Workspace\\Euler\\euler8.txt", "r")
    readFile = openFile.read().splitlines() 
    maxProduct = -1
    number = ""
    for line in readFile:
        number += line
    for i in range(len(number) - (n - 1)):
        product = 1
        for j in range(n):
            product *= int(number[i + j])
            if product > maxProduct:
                maxProduct = product
    print(maxProduct)

# Summation of primes
def euler10(ub):
    candidates = [True for i in range(ub - 1)]
    total = 0
    for j in range(2, math.floor(math.sqrt(ub - 1)) + 1):
        for k in range(2, (ub - 1) // j + 1):
            candidates[j * k - 1] = False
    for l in range(len(candidates)):
        if candidates[l]:
            total += l + 1
    print("Sum of all primes below " + str(ub) + " = " + str(total - 1)) # - 1 compensates for 1 being nonprime

# Maximum path sum I
def euler18(fileAddress = "C:\\Users\\RR\\Dropbox\\CS\\Workspace\\Euler\\euler18.txt"):
    with open(fileAddress) as file:
        arr = [line.split() for line in file]
        for r in range(len(arr) - 2, -1, -1):
            for c in range(0, len(arr[r])):
                if int(arr[r + 1][c]) >= int(arr[r + 1][c + 1]):
                    arr[r][c] = int(arr[r][c]) + int(arr[r + 1][c])
                else:
                    arr[r][c] = int(arr[r][c]) + int(arr[r + 1][c + 1])
        print(arr[0][0])

# Amicable numbers
def euler21(ub):
    total = 0
    for i in range(1, ub):
        if d(d(i)) == i and i != d(i):
            print("True")
            total += i
    print(total)

# euler21 helper method
def d(n):
    total = 1
    for i in range(2, math.floor(math.sqrt(n)) + 1):
        if n % i == 0:
            total += i
            if n / i != i:
                total += n / i
    return(total)

# Non-abundant sums
# Inefficient!
def euler23():
    total = 0
    for n in range(1, 28124):
        can = False
        for i in range(1, n // 2 + 1):
            if abundant(i) and abundant(n - i):
                can = True
                break
        if not can:
            total += n
        print(total)
    print(total)
    
# euler23 helper method
def abundant(n):
    return(d(n) > n)

# Lexicographic permutations
# Permutes numbers formed with digits from 0 to n
def euler24(n, pos):
    list = [i for i in range(0, n + 1)]
    arr = sorted(itertools.permutations(list))
    print(arr[pos - 1])
    
# Reciprocal cycles
def euler26(ub):
    maxRepeat = -1
    maxP = -1
    for p in primes[:ub]:
        n = intLen(p) - 1
        print("prime: " + str(p) + " n: " + str(n))
        strNum = str(1 / p)[n + 2: n + 2 + 2 * p]
        r = 1
        modStr = ""
        for i in range(int(strNum) // r):
            modStr + strNum[:r]
        print(modStr)
        print(strNum)

# Quadratic primes
def euler27(bound = 1000):
    maxLength = 0                   # Maximum number of primes for consecutive values of n
    maxProduct = 0
    for a in range(-bound + 1, bound):
        for b in range(-bound + 1, bound):
            length = 0
            n = 0
            while isPrime(n ** 2 + a * n + b):
                length += 1
                n += 1
            if length > maxLength:
                maxLength = length
                maxProduct = a * b
    print(maxProduct)
                
# Coin sums
def euler31(n):
    combos = 0
    for p200 in range(n // 200 + 1):
        for p100 in range((n - p200 * 200) // 100 + 1):
            for p50 in range((n - p200 * 200 - p100 * 100) // 50 + 1):
                for p20 in range((n - p200 * 200 - p100 * 100 - p50 * 50) // 20 + 1):
                    for p10 in range((n - p200 * 200 - p100 * 100 - p50 * 50 - p20 * 20) // 10 + 1):
                        for p5 in range((n - p200 * 200 - p100 * 100 - p50 * 50 - p20 * 20 - p10 * 10) // 5 + 1):
                            for p2 in range((n - p200 * 200 - p100 * 100 - p50 * 50 - p20 * 20 - p10 * 10 - p5 * 5) // 2 + 1):
                                for p1 in range(n - p200 * 200 - p100 * 100 - p50 * 50 - p20 * 20 - p10 * 10 - p5 * 5 - p2 * 2 + 1):
                                    if n == p200 * 200 + p100 * 100 + p50 * 50 + p20 * 20 + p10 * 10 + p5 * 5 + p2 * 2 + p1:
                                        combos += 1
    print(str(combos) + " ways to form " + str(n) + " pence.")

# Pandigital products
def euler32():
    pandigital = [str(n) for n in range(1, 10)]
    pandigitalProducts = set()
    for f1 in range(1, 100):
        lb = math.ceil(1000 / f1)
        ub = math.floor(10000 / f1)
        for f2 in range(lb, ub):
            product = f1 * f2
            if pandigital == sorted(str(f1) + str(f2) + str(product)):
                pandigitalProducts.add(product)
    print(sum(pandigitalProducts))

# Digit cancelling fractions
# Given: numerator = n1n2 and denominator = d1d2
def euler33():
    productNum = 1
    productDen = 1
    for n1 in range(1, 10):
        for n2 in range(10):
            for d1 in range(1, 10):
                for d2 in range(10):
                    if n1 * 10 + n2 < d1 * 10 + d2 and not (n2 == 0 and d2 == 0):
                        frac = (n1 * 10 + n2) / (d1 * 10 + d2)
                        num = n1 * 10 + n2
                        den = d1 * 10 + d2
                        strFrac = str(n1 * 10 + n2) + "/" + str(d1 * 10 + d2)
                        if n1 == d1 and d2 != 0:
                            if frac == n2 / d2:
                                print(strFrac)
                                productNum *= num
                                productDen *= den
                        if n1 == d2 and d1 != 0:
                            if frac == n2 / d1:
                                print(strFrac)
                                productNum *= num
                                productDen *= den
                        if n2 == d1 and d2 != 0:
                            if frac == n1 / d2:
                                print(strFrac)
                                productNum *= num
                                productDen *= den
                        if n2 == d2 and d1 != 0:
                            if frac == n1 / d1:
                                print(strFrac)
                                productNum *= num
                                productDen *= den
    factor = fractions.gcd(productNum, productDen)
    print(productDen / factor)
       
# Circular primes
def euler35(ub):
    count = 0
    for n in range(1, ub):
        num = str(n)
        allPrime = True
        for i in range(len(num)):
            if not isPrime(num[i:] + num[:i]):
                allPrime = False
                break
        if allPrime:
            count += 1
    print(count)

# Double-base palindromes
def euler36(ub):
    total = 0
    for n in range(1, ub):
        if isPalindrome(n) and isPalindrome("{0:b}".format(n)):
            total += n
    print(total)

# Truncatable primes
def euler37(ub):
    total = 0
    for p in primes[4:ub]:
        testLR = p % 10 ** (len(str(p)) - 1)
        testRL = p // 10
        lR = True
        rL = True
        # Test right to left
        while testRL > 0:
            if not isPrime(testRL):
                rL = False
                break
            testRL //= 10
        # Test left to right
        while testLR > 0:
            if not isPrime(testLR):
                lR = False
                break
            testLR = testLR % 10 ** (len(str(testLR)) - 1)
        if lR and rL:
            total += p
            print(p)
    print(total)

# Pandigital multiples
# Inefficient
def euler38():
    maxProduct = -1
    pandigital = [str(i) for i in range(1, 10)]
    for n in range(1, 10):
        for x in range(1, 10000):
            num = ""
            for i in range(1, n + 1):
                num += str(x * i)
            if sorted(list(num)) == pandigital and int(num) > maxProduct:
                maxProduct = int(num)
    print(maxProduct)
# Pandigital prime
def euler41():
    for i in range(len(primes)):
        l = len(str(primes[-(i + 1)]))
        arr = [str(n) for n in range(1, l + 1)]
        if sorted(list(str(primes[-(i + 1)]))) == arr:
            return(primes[-(i + 1)])
    return(-1)

# Triangular, pentagonal, and hexagonal
# Begins at lb = 144 since it is given that n = 143 works
def euler45(ub = 100000):
    for n in range(144, ub):
        h = n * (2 * n - 1)
        for rootP in quadraticRoots(3, -1, -2 * h):
            if rootP > 0 and rootP % 1 == 0:
                for rootT in quadraticRoots(1, 1, -2 * h):
                    if rootT > 0 and rootT % 1 == 0:
                        print(h)

# Consecutive prime sum
# 78499th prime is the first prime above one million
def euler50(ub):
    # Determine upper bound for number of terms to sum
    i = 0
    total = 0
    while total + primes[i + 1] < primes[ub]:
        total += primes[i]
        i += 1
    print(i + 1)
    maxLength = 1
    n = 1
    while n < ub:
        for start in range(ub - n):
            end = start + 1
            while 
            for end in range(start + 1, start + 1 + 546):
                print(sum(primes[start:end]))
                if sum(primes[start:end]) == primes[ub - n]:
                    if end - start > maxLength:
                        maxLength = end - start
        n += 1
    print(maxLength)
euler50(78499)
# Permuted multiples
# euler52()
def euler52():
    n = 123456
    while True:
        if sorted(str(n)) == sorted(str(2 * n)) and sorted(str(2 * n)) == sorted(str(3 * n)) and sorted(str(3 * n)) == sorted(str(4 * n)) and sorted(str(4 * n)) == sorted(str(5 * n)) and sorted(str(5 * n)) == sorted(str(6 * n)):
            break
        n += 1
    print(n)

# Maximum path sum II
def euler67():
    euler18("C:\\Users\\RR\\Dropbox\\CS\\Workspace\\Euler\\euler67.txt")
    
################################# 
