#!/bin/bash

# Pressure measured across the spirometer is
#   P = Y = (1-H(t-x))A1e^(P1t) + H(t-x)A2e^(P2t)
# where x is the time that the decay lines cross

from sympy import *
from sympy.abc import s,t
import dill

def generate_equation(y, u):
    '''
    generate_equation(y, u)
    y = function(output(t))
    u = function(input(t))
    Generates time-domain equation relating
    input pressure to output pressure.
    '''
    A1 = Symbol("A1")
    A2 = Symbol("A2")
    P1 = Symbol("P1")
    P2 = Symbol("P2")

    # x is the time that the curves cross
    x = Symbol("x")

    # Laplace of pressure equation
    Y = A1/(s+P1) - A1*exp(-x*s)/(s+P1) + A2*exp(-x*s)/(s+P2)

    # Cross multiply the fractions
    Y = together(Y)

    # Separate the fraction
    numerator = fraction(Y)[0]
    denominator = fraction(Y)[1]

    # Add "fudge factor"
    # (easier substitution to get back to time-domain)
    numerator = numerator*s
    denominator = denominator*s

    # Expand all brackets
    numerator = numerator.expand()
    denominator = denominator.expand()

    # Collect powers of s
    numerator = collect(numerator, s)
    denominator = collect(denominator, s)

    y1 = diff(y, t)
    y2 = diff(y, t, t)
    y3 = diff(y, t, t, t)
    u1 = diff(u, t)
    u2 = diff(u, t, t)
    u3 = diff(u, t, t, t)

    def convert_to_time_domain(func, v, v1, v2, v3):
        # Assumes func was multiplied by s
        func_in_time_domain = func.subs(exp(-x*s), Heaviside(x))
        func_in_time_domain = func_in_time_domain.subs(s**4, v3)
        func_in_time_domain = func_in_time_domain.subs(s**3, v2)
        func_in_time_domain = func_in_time_domain.subs(s**2, v1)
        func_in_time_domain = func_in_time_domain.subs(s, v)
        return(func_in_time_domain)

    RHS = convert_to_time_domain(numerator, u, u1, u2, u3)
    LHS = convert_to_time_domain(denominator, y, y1, y2, y3)

    #print("Input")
    #pprint(simplify(RHS))
    #print("\nOutput")
    #pprint(simplify(LHS))

    return([LHS, RHS])

def solve_equation(LHS, RHS, u):
    '''
    Given time domain LHS (array of values) and RHS,
    solve for ODE(var) in RHS
    '''
    # Rearrange equation to equal zero
    eq = Eq(RHS - LHS)
    solved_equation = dsolve(eq, u)

    return(solved_equation)


def save_equations():
    y = Function("y")(t)
    u = Function("u")(t)
    x = Symbol("x")

    [LHS, RHS] = generate_equation(y, u)
    time_dependent_input = solve_equation(y, RHS*Heaviside(x), u)

    # Substitute into EoM_y here and get values
    y_of_t_minus_2= Symbol("y2t")
    y_of_t_minus_1 = Symbol("y1t")
    y_of_t = Symbol("yt")
    y_of_t_plus_1= Symbol("yt1")
    Ht = Symbol("Ht")
    u_of_t = Symbol("ut")
    dt = Symbol("dt")

    # Approximations to derivatives
    dy = (y_of_t - y_of_t_minus_1)/dt
    d2y = (y_of_t - 2*y_of_t_plus_1 + y_of_t_minus_2)/dt**2

    Approx_EoM_y = LHS.subs({Derivative(y,t,2):d2y})
    Approx_EoM_y = Approx_EoM_y.subs({Derivative(y,t):dy})
    Approx_EoM_y = Approx_EoM_y.subs({y:y_of_t})

    new_u = time_dependent_input.subs(u,u_of_t)
    new_u = new_u.subs(y,y_of_t)

    filename = 'symbolicEquations.pk'
    with open(filename, 'wb') as fileobject:
        dill.dump([Approx_EoM_y,new_u], fileobject)
