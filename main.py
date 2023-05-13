from decimal import *
import gmpy2, time, math
from gmpy2 import mpz, mpfr

out = open("results.txt", "a+")

def machin(digits):
    getcontext().prec = digits + 10     
    scale = 10**(digits + 10)           
        
    def atan(x):

        current_value = 0
        divisor = 1
        x_squared = x * x
        current_term = scale // x    

        while True:
            current_value += current_term // divisor
            divisor += 2
            current_term //= -x_squared

            if current_term == 0:
                break
        return Decimal(current_value) / scale
    
    pi_fourths = 4 * atan(5) - atan(239)
    return pi_fourths * 4
    
def chudnovsky(digits):
    scale = 10**(mpz(digits+20))  

    gmpy2.get_context().precision = int(math.log2(10) * (digits + 20)) 
    
    k = mpz(1)
    a_k = scale
    a_sum = scale
    b_sum = mpz(0)
    C = mpz(640320)
    C_cubed_over_24 = C**3 // 24
    
    while True:
        a_k *= -(6*k-5) * (2*k-1) * (6*k-1)
        a_k //= k**3 * C_cubed_over_24
        a_sum += a_k
        b_sum += k * a_k
        k += 1
        if a_k == 0:
            break
    
    total_sum = mpfr(13591409 * a_sum + 545140134 * b_sum)
    pi = (426880 * gmpy2.sqrt(mpfr(10005))) / total_sum
 
    return pi*scale

def gauss_legendre(digits):
    gmpy2.get_context().precision = int(math.log2(10) * (digits + 5))
    
    iterations = int(math.log2(digits)) 
    
    a = mpfr(1)
    b = 1 / gmpy2.sqrt(mpfr(2))
    t = mpfr(1)/4
    x = mpfr(1)

    for i in range(iterations):
        y = a
        a = (a + b) / 2
        b = gmpy2.sqrt(b * y)
        t = t - x * (y-a)**2
        x = 2 * x
        
    pi = ((a+b)**2) / (4*t)
    return pi
     

def main():
    with open("results.txt", "a+") as out:
        for x in range(1000, 100000, 500):
            start_time = time.time()      
            # chudnovsky(x)
            # machin(x)
            #gauss_legendre(x)
            #elapsed_time = (time.time() - start_time) * 1000
            out.write(str(x) + "\n")

    
main()
