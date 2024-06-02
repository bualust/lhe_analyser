ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      written by the UFO converter
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE COUP1()

      IMPLICIT NONE
      INCLUDE 'model_functions.inc'

      DOUBLE PRECISION PI, ZERO
      PARAMETER  (PI=3.141592653589793D0)
      PARAMETER  (ZERO=0D0)
      INCLUDE 'input.inc'
      INCLUDE 'coupl.inc'
      GC_1 = (MDL_EE*MDL_COMPLEXI)/3.000000D+00
      GC_2 = (-2.000000D+00*MDL_EE*MDL_COMPLEXI)/3.000000D+00
      GC_3 = -(MDL_EE*MDL_COMPLEXI)
      GC_135 = -(MDL_EE*MDL_COMPLEXI)/(2.000000D+00*MDL_CTH*MDL_STH)
      GC_195 = -(MDL_EE*MDL_COMPLEXI*MDL_STH)/(3.000000D+00*MDL_CTH)
      GC_196 = (2.000000D+00*MDL_EE*MDL_COMPLEXI*MDL_STH)/(3.000000D
     $ +00*MDL_CTH)
      GC_197 = -((MDL_EE*MDL_COMPLEXI*MDL_STH)/MDL_CTH)
      GC_223 = (4.000000D+00*MDL_COMPLEXI*MDL_GHAA)/MDL_VEVHAT
      GC_227 = (2.000000D+00*MDL_COMPLEXI*MDL_GHZA)/MDL_VEVHAT
      GC_284 = (MDL_EE__EXP__2*MDL_COMPLEXI*MDL_VEVHAT)/(2.000000D+00
     $ *MDL_CTH__EXP__2*MDL_STH__EXP__2)
      GC_451 = -((MDL_COMPLEXI*MDL_YB)/MDL_SQRT__2)
      GC_503 = -((MDL_CDHIM*MDL_VEVHAT__EXP__2*MDL_YB)
     $ /(MDL_LAMBDASMEFT__EXP__2*MDL_SQRT__2))
      GC_504 = (MDL_CDHRE*MDL_COMPLEXI*MDL_VEVHAT__EXP__2*MDL_YB)
     $ /(MDL_LAMBDASMEFT__EXP__2*MDL_SQRT__2)
      END
