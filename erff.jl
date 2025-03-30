
# Float32 erf(x) function 
#
# Julia implementation of https://github.com/ARM-software/optimized-routines/blob/master/math/erff.c
# 
# License: MIT



# Minimax approximation of erff. 
const A=(0x1.06eba6p-03, -0x1.8126e0p-02, 0x1.ce1a46p-04, -0x1.b68bd2p-06, 0x1.47348p-08, -0x1.3a1a82p-11)
const B=(0x1.079d0cp-3, 0x1.450aa0p-1, 0x1.b55cb0p-4, -0x1.8d6300p-6, 0x1.d1336p-9, -0x1.91d2ccp-12, 0x1.222900p-16)

const TwoOverSqrtPiMinusOne=0x1.06eba8p-3



#  Efficient implementation of erff
#  using either a pure polynomial approximation or
#  the exponential of a polynomial.
#  Worst-case error is 1.09ulps at 0x1.c111acp-1.  

function erff(x::Float32)

#   Reinterpreting float binary as unsigned Int32
    ix::UInt32=reinterpret(UInt32,x)
    sign::Bool=x<0
    ia12=(ix>>20) & 0x7ff


#    Limit of both intervals is 0.875 for performance reasons but coefficients
#    computed on [0.0, 0.921875] and [0.921875, 4.0], which brought accuracy
#    from 0.94 to 1.1ulps.  
  if (ia12 < 0x3f6)
    #   a = |x| < 0.875.  

    x2 = x * x

    #    Normalized cases (|x| < 0.921875). Use Horner scheme for x+x*P(x^2).  

        r = A[6]
        r = fma(r, x2, A[5])
        r = fma(r, x2, A[4])
        r = fma(r, x2, A[3])
        r = fma(r, x2, A[2])
        r = fma(r, x2, A[1])
        r = fma(r, x, x)
      
  elseif (ia12 < 0x408)
    #  |x| < 4.0 - Use a custom Estrin scheme.  

        a = abs(x)
    #  Start with Estrin scheme on high order (small magnitude) coefficients.  
        r = fma(B[7], a, B[6])
        u = fma(B[5], a, B[4])
        x2 = x * x
        r = fma(r, x2, u)

    #  Then switch to pure Horner scheme.  
        r = fma(r, a, B[3])
        r = fma(r, a, B[2])
        r = fma(r, a, B[1])
        r = fma(r, a, a)

    #  Single precision exponential with ~0.5ulps,
	#  ensures erff has max. rel. error
	#  < 1ulp on [0.921875, 4.0],
	#  < 1.1ulps on [0.875, 4.0].  
      r = exp(-r)

    #    Explicit copysign (calling copysignf increases latency).  
      if (sign)
	r = -1.0 + r
      else
	r = 1.0 - r
      end
  else
    #  |x| >= 4.0.  
        #    Explicit copysign (calling copysignf increases latency).  
        if (sign)
        r = -1.0
        else
        r = 1.0
        end 
    end
  return r



end












