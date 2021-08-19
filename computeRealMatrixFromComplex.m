function X = computeRealMatrixfromComplex(A)

    X = 1/2*(A + A.') + 1/(2*sqrt(-1))*(A - A.');