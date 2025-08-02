# Float64 erf(x) function 
#
# Julia implementation of https://github.com/ARM-software/optimized-routines/blob/master/math/erf.c
# 
# License: MIT




#    Fast erf implementation using a mix of
#    rational and polynomial approximations.
#    Highest measured error is 1.01 ULPs at 0x1.39956ac43382fp+0.  
function erf(x::Float64)

    # # top 32 bits 
    ix::UInt32=reinterpret(UInt64,x)>>32
    # # top 32, without sign bit 
    ia::UInt32=ix & 0x7fffffff
    # # sign
    # sign::UInt32=ix>>31

    sign::Bool=x<0



    if (ia < 0x3feb0000)
    #  a = |x| < 0.84375.  

        x2 = x * x

        if (ia < 0x3fe00000)
        ## a < 0.5  - Use polynomial approximation.  

            # Minimax approximation of erf of the form x*P(x^2) approximately on the interval [0;0.5]
            PA=(0.12837916709551256, -0.3761263890318287, 0.11283791670896592, -0.026866170630075903, 0.005223977428649761, -0.0008548312229989974, 0.00012054654502151287, -1.4906315067498891e-5, 1.6126444349070824e-6, -1.3074259056679966e-7)
        
            r= fma(x,evalpoly(x2,PA),x) ## This fma is crucial for accuracy.  

            return r
        else
        ## 0.5 <= a < 0.84375 - Use rational approximation.  

            # Rational approximation on [0x1p-28, 0.84375] 
            NA=(0x1.06eba8214db68p-3, -0x1.4cd7d691cb913p-2, -0x1.d2a51dbd7194fp-6,-0x1.7a291236668e4p-8, -0x1.8ead6120016acp-16)
            DA=(1,0x1.97779cddadc09p-2, 0x1.0a54c5536cebap-4, 0x1.4d022c4d36b0fp-8,0x1.15dc9221c1a1p-13, -0x1.09c4342a2612p-18)

            P=evalpoly(x2,NA)
            Q=evalpoly(x2,DA)

            return fma(P / Q, x, x)
        end
    elseif (ia < 0x3ff40000)
    ## 0.84375 <= |x| < 1.25.  

        # Rational approximation on [0.84375, 1.25] 
        NB=( -0x1.359b8bef77538p-9, 0x1.a8d00ad92b34dp-2, -0x1.7d240fbb8c3f1p-2, 0x1.45fca805120e4p-2, -0x1.c63983d3e28ecp-4, 0x1.22a36599795ebp-5, -0x1.1bf380a96073fp-9 )
        DB=( 1, 0x1.b3e6618eee323p-4, 0x1.14af092eb6f33p-1, 0x1.2635cd99fe9a7p-4, 0x1.02660e763351fp-3, 0x1.bedc26b51dd1cp-7, 0x1.88b545735151dp-7 )

        C = 0x1.b0ac16p-1

        a = abs(x) - 1.0

        P=evalpoly(a,NB)
        Q=evalpoly(a,DB)


        if (sign)
            return -C - P / Q
        else
            return C + P / Q
        end
    elseif (ia < 0x40000000)
    ## 1.25 <= |x| < 2.0.  
        a = abs(x)
        a = a - 1.25
        
        # erfc polynomial approximation
        # Generated using Sollya::remez(f(c*x+d), deg, [(a-d)/c;(b-d)/c], 1, 1e-16), [|D ...|] with deg=15 a=1.25 b=2 c=1 d=1.25 
        PC=( 0x1.3bcd133aa0ffcp-4, -0x1.e4652fadcb702p-3, 0x1.2ebf3dcca0446p-2, -0x1.571d01c62d66p-3, 0x1.93a9a8f5b3413p-8, 0x1.8281cbcc2cd52p-5, -0x1.5cffd86b4de16p-6, -0x1.db4ccf595053ep-9, 0x1.757cbf8684edap-8, -0x1.ce7dfd2a9e56ap-11, -0x1.99ee3bc5a3263p-11, 0x1.3c57cf9213f5fp-12, 0x1.60692996bf254p-14, -0x1.6e44cb7c1fa2ap-14, 0x1.9d4484ac482b2p-16, -0x1.578c9e375d37p-19)

        # Obtains erfc of |x|
        r=evalpoly(a,PC)

        if (sign)
            return -1.0 + r
        else
            return 1.0 - r
        end
    elseif (ia < 0x400a0000)
    ## 2 <= |x| < 3.25.  
        a = abs(x)
        a = fma(0.5, a, -1.0)

        # Generated using Sollya::remez(f(c*x+d), deg, [(a-d)/c;(b-d)/c], 1, 1e-16), [|D ...|] with deg=17 a=2 b=3.25 c=2 d=2 
        PD=( 0x1.328f5ec350e5p-8, -0x1.529b9e8cf8e99p-5, 0x1.529b9e8cd9e71p-3, -0x1.8b0ae3a023bf2p-2, 0x1.1a2c592599d82p-1, -0x1.ace732477e494p-2, -0x1.e1a06a27920ffp-6, 0x1.bae92a6d27af6p-2, -0x1.a15470fcf5ce7p-2, 0x1.bafe45d18e213p-6, 0x1.0d950680d199ap-2, -0x1.8c9481e8f22e3p-3, -0x1.158450ed5c899p-4, 0x1.c01f2973b44p-3, -0x1.73ed2827546a7p-3, 0x1.47733687d1ff7p-4, -0x1.2dec70d00b8e1p-6, 0x1.a947ab83cd4fp-10 )


        # Obtains erfc of |x|
        r=evalpoly(a,PD)


        if (sign)
            return -1.0 + r
        else
            return 1.0 - r
        end
    elseif (ia < 0x40100000)
    ## 3.25 <= |x| < 4.0.  
        a = abs(x)
        a = a - 3.25

        # Generated using Sollya::remez(f(c*x+d), deg, [(a-d)/c;(b-d)/c], 1, 1e-16), [|D ...|] with deg=13 a=3.25 b=4 c=1 d=3.25 
        PE=( 0x1.20c13035539e4p-18, -0x1.e9b5e8d16df7ep-16, 0x1.8de3cd4733bf9p-14, -0x1.9aa48beb8382fp-13, 0x1.2c7d713370a9fp-12, -0x1.490b12110b9e2p-12, 0x1.1459c5d989d23p-12, -0x1.64b28e9f1269p-13, 0x1.57c76d9d05cf8p-14, -0x1.bf271d9951cf8p-16, 0x1.db7ea4d4535c9p-19, 0x1.91c2e102d5e49p-20, -0x1.e9f0826c2149ep-21, 0x1.60eebaea236e1p-23 )

        r=evalpoly(a,PE)


        if (sign)
            return -1.0 + r
        else
            return 1.0 - r
        end
    elseif (ia < 0x4017a000)
    ## 4 <= |x| < 5.90625.  
        a = abs(x)
        a = fma(0.5, a, -2.0)

        # Generated using Sollya::remez(f(c*x+d), deg, [(a-d)/c;(b-d)/c], 1, 1e-16), [|D ...|] with deg=16 a=4 b=5.90625 c=2 d=4 
        PF=( 0x1.08ddd130d1fa6p-26, -0x1.10b146f59ff06p-22, 0x1.10b135328b7b2p-19, -0x1.6039988e7575fp-17, 0x1.497d365e19367p-15, -0x1.da48d9afac83ep-14, 0x1.1024c9b1fbb48p-12, -0x1.fc962e7066272p-12, 0x1.87297282d4651p-11, -0x1.f057b255f8c59p-11, 0x1.0228d0eee063p-10, -0x1.b1b21b84ec41cp-11, 0x1.1ead8ae9e1253p-11, -0x1.1e708fba37fccp-12, 0x1.9559363991edap-14, -0x1.68c827b783d9cp-16, 0x1.2ec4adeccf4a2p-19 )

        r=evalpoly(a,PF)


        if (sign)
            return -1.0 + r
        else
            return 1.0 - r
        end
    else
        

        if (sign)
            return -1.0
        else
            return 1.0
        end

    end



end

function erfc(x::Float64)
    return 1-erf(x)
end

