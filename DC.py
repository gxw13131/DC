'''
a script for computing Differential Coefficients.
'''
import numpy as np
from sympy import Integer, Rational, Matrix


def DC_sym(N, x_r):
    # x_r = Integer(x_r)
    postion_list = Matrix([i for i in range(x_r, x_r+N)])
    vandermonde = Matrix([[pos**i for i in range(N)]
                          for pos in postion_list])
    Coeff_M = Matrix([Rational(1, np.math.factorial(i)) for i in range(N)])
    C = Matrix([[vandermonde[i, j]*Coeff_M[j]
                 for j in range(N)] for i in range(N)])
    # print(C)
    inv_C = C.inv()
    # print(inv_C)
    return C, inv_C


def flux_sym(inv_C):
    N = inv_C.shape[0]
    Coeff = inv_C[1, :]
    b = [-Coeff[0]]
    for i in range(1, N):
        b.append(b[-1]-Coeff[i])
    return Matrix(b)


# def DC_np(n, left):
#     postion_list = [i for i in range(index_left, index_left+N)]
#     vandermonde = np.array([[pos**i for i in range(N)]
#                             for pos in postion_list])
#     Coeff_M = np.array([1/np.math.factorial(i) for i in range(N)])
#     C = np.matrix(vandermonde*Coeff_M)
#     inv_C = C.I
#     return C, inv_C


def print_invC(C_inv, x_begin):
    N = C_inv.shape[0]
    print('-'*80)
    for i in range(N):
        print('\n f{0:s} = '.format('\''*(i)+' '*(N-1-i)), end='')
        for j in range(0, N):
            if C_inv[i, j] != 0:
                print(
                    ' + {0:>6s}*f({1:2d})'.format(str(C_inv[i, j]), x_begin+j), end='')
    print('\n'+'-'*80)


def print_flux_R(flux_R, x_begin, x_base):
    N = flux_R.shape[0]
    print('\n f({0:d}+1/2) = '.format(x_base), end='')
    for j in range(0, N):
        if flux_R[j] != 0:
            print(
                ' + {0:>6s}*f({1:2d})'.format(str(flux_R[j]), x_begin+j), end='')
    print('\n'+'-'*80)


def read_params():
    N = input('Enter the order of scheme:')
    N = int(N)
    assert(N>0),'the order of scheme musut be larger than 1 !!'
    x_begin = input('Enter the begin index of stencil points:')
    x_begin = int(x_begin)
    x_base = input('Enter the base index of stencil(default value is 0):')
    x_base = 0 if len(x_base) == 0 else int(x_base)
    return N, x_begin, x_base


if __name__ == "__main__":
    while(1):
        try:
            N, x_begin, x_base = read_params()
            C, inv_C = DC_sym(N+1, x_begin-x_base)
            print_invC(inv_C, x_begin)
            flux_R = flux_sym(inv_C)
            print_flux_R(flux_R, x_begin, x_base)
            # eval_sym(inv_C, x_sample-x_base)
            print('press Ctrl+C to quit...')
        except KeyboardInterrupt:
            quit()
        except AssertionError:
            continue
