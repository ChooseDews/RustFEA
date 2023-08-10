

import sympy


def derive_eight_node():
    xi, eta, zeta = sympy.symbols('xi eta zeta') #local coordinates
    N1 = (1-xi)*(1-eta)*(1-zeta)/8
    N2 = (1+xi)*(1-eta)*(1-zeta)/8
    N3 = (1+xi)*(1+eta)*(1-zeta)/8
    N4 = (1-xi)*(1+eta)*(1-zeta)/8
    N5 = (1-xi)*(1-eta)*(1+zeta)/8
    N6 = (1+xi)*(1-eta)*(1+zeta)/8
    N7 = (1+xi)*(1+eta)*(1+zeta)/8
    N8 = (1-xi)*(1+eta)*(1+zeta)/8
    N = sympy.Matrix([N1, N2, N3, N4, N5, N6, N7, N8])
    dN = N.jacobian([xi, eta, zeta])
    return N, dN

def test_derive_eight_node():
    N,dN = derive_eight_node()
    print(N)
    print("Derivatives:")
    print(sympy.N(sympy.simplify(dN), 5))

    #

if __name__ == '__main__':
    test_derive_eight_node()
