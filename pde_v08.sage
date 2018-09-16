###############################################################################
#   Solution of Partial Differential Equation by Oberst-Riquier Algorithm     #
#                                                                             #
#   Author: Ayan Sengupta, Debasattam Pal                                     #
#   Version: 1.0                                                              #
###############################################################################
reset()  # Resets all variables previously defined

import timeit as tm
import re

exp = 10


'''
Here Oberst-Riquier algorithm is used for solving the given set of PDE's. The
output of this algorithm is in a power series form. The program first asks for
the number of independent variables, after that it asks for the partial
differential equation(s). After that the initial conditions must be given as
specified. The final solution will be saved in a text file in the current
directory named 'output.txt'. To see the solution in the terminal type
'solution'. Once the number of independent variables are defined, variables
in d i.e. d1,d2,d3...and variables in x i.e. x1,x2,x3... are created as
required. Here di signifies partial derivative with respect to xi.
Now, the user needs to give the PDE's in polynomial form in terms of d1,d2,d3..

* added support for reading equations from a file.
* removed manual entry of pdes and number of variables in terminal.

'''
sys.path.append('~/Documents/Project/SAGE/pde_sem1_2016')

def Polynomial_ring(num, precision, mono_order='lex'):
    ''' This Function defines a polynomial ring with variables d1,d2,d3.... The number of variables defined will be given by the user'''
    RF = RealField(precision)  # having problems with fixed precision
    P = PolynomialRing(
        QQ, [
            'd%s' %
            p for p in range(
                1, num + 1)], order=mono_order)
    P.inject_variables()
    return P





def Groebner(eq):
    '''This function takes in the polynomials given by the user, generates the ideal and finally computes the reduced groebner basis'''
    global start_time
    start_time = tm.default_timer()
    I = ideal(eq)
    groebner = I.groebner_basis()
    return groebner


def generate_dv(exp):
    '''generates all the possible monomials with max degree of each variable defined by the user, d1 d2 with exp3 will generate 15 terms d1,d1^2,d1^3,d*1d2,d1^2*d2...d1^3*d1^3'''

    var_list = P.gens_dict().values()
    num = [exp + 1] * len(var_list)
    d_v = monomials(var_list, num)
    d_v.remove(1)
    #print("d_v generated...")
    return d_v


def standard_monomial():
    '''This function generates a set of standard monomials'''
    T = dict()
    T_values = []
    T_keys = []
    for i in range(len(d_v)):
        if d_v[i].reduce(groebner) != 0:
            # print d_v[i]
            T[d_v[i]] = d_v[i].reduce(groebner)

    T_keys = T.keys()
    T_values = T.values()
    T_values_coeff = [T_values[i].coefficients()[0]
                      for i in range(len(T_values))]

    std_monomial_list = [set(T_values[i].monomials())
                         for i in range(len(T_values))]
    std_monomial_set = set().union(*std_monomial_list)
    std_monomial = sorted(list(std_monomial_set))
    return T, T_values, T_keys, T_values_coeff, std_monomial


def standard_monomial_evaluated_over_initial_condition(
        initial_condition, std_monomial):

    std_mon_coeff_from_initial_condition = {}
    for i in range(len(std_monomial)):
        derivative_order = list(std_monomial[i].degrees())
        diff_parameter = [[x_var[k]] * derivative_order[k]
                          for k in range(len(derivative_order)) if derivative_order[k] != 0]
        flat_diff_parameter = [
            item for sublist in diff_parameter for item in sublist]
        coeff_temp = initial_condition.derivative(flat_diff_parameter)
        for p in range(len(x_var)):
            m = x_var[p]
            coeff_temp = coeff_temp.subs(m == 0)
        std_mon_coeff_from_initial_condition[std_monomial[i]] = coeff_temp

    return std_mon_coeff_from_initial_condition


def generate_var_x(num):
    '''Defines the variables in terms of x1, x2, x3, ... '''
    x_var = [var('x%d' % i) for i in range(1, num + 1)]
    return x_var


def generate_factorial(monomial):
    fac_temp = 1
    for p in range(len(monomial.degrees())):
        fac_temp = fac_temp * factorial(monomial.degrees()[p])
    return fac_temp


def coeff_temp(f, n):
    k = len(value[0].degrees())
    v = x_var
    for i in range(k):
        z = v[i]
        g = f.derivative(z, value[n].degrees()[i])
        f = g.subs(z == 0)

    return f


def calculate_a0(f):
    for i in range(len(value[0].degrees())):
        z = x_var[i]
        f = f.subs(z == 0)
    return f


def generate_solution(num,
        std_mon_coeff_from_initial_condition,
        initial_condition,
        prec,
        x_var):

    coeff_temp = []
    for i in range(len(value)):
        sub_monomial = value[i].monomials()
        sub_coefficient = value[i].coefficients()
        sub_monomial_value = [
            std_mon_coeff_from_initial_condition[
                sub_monomial[k]] for k in range(
                len(sub_monomial))]
        total_sub_monomial_sum = sum(
            [(a * b) for a, b in zip(sub_coefficient, sub_monomial_value)])
        coeff_temp.append(total_sub_monomial_sum / generate_factorial(key[i]))
    #solution_temp = sum([(a * b) for a, b in zip(coeff_temp, key)])
    #solution_temp = solution_temp + calculate_a0(initial_condition)
    solution_temp_approx = sum([(round(a, 3) * b)
                                for a, b in zip(coeff_temp, key)])
    solution_temp_approx = solution_temp_approx + \
        calculate_a0(initial_condition)
    solution_approx = solution_temp_approx.subs(
        {P.gens()[i]: x_var[i] for i in range(0, num)})

    solution_temp = [(a * b) for a, b in zip(coeff_temp, key)]
    solution_temp.append(calculate_a0(initial_condition))
    solution_sum = sum(solution_temp)
    solution = solution_sum.subs(
        {P.gens()[i]: x_var[i] for i in range(0, num)})

    return solution, solution_approx


def pde_solution_check(poly, solution):

    return nil


def print_initial_condition(std_monomial):

    value_for_infinity = 9

    monomial_deg = sum(std_monomial).degrees()
    monomial_var = sum(std_monomial).variables()
    #exp_initial_condition = [0]*len(monomial_var)
    for i in range(len(monomial_deg)):
        if monomial_deg[i] > value_for_infinity:
          #  exp_initial_condition[i] = value_for_infinity+1
            print 'P(a*{0})*exp(b*{1})'.format(x_var[i], x_var[i]),
        else:
            pass


def initial_condition_var(num,list_mono):
    l = max({(value[i].variables()) for i in range(len(value))})
    x_initial = [l[j].subs({P.gens()[i]:x_var[i]
                            for i in range(0, num)}) for j in range(len(l))]
    return x_initial


#############################################################
#################### MAIN PROGRAM ###########################
#############################################################

#num_var = input("Number of independent variables: ")
#o = raw_input("monomial ordering(lex, deglex, degrevlex ...): ")
fname = "example3.txt"
try:
    with open(fname) as f:
        content = f.readlines()
except IOError:
    print "can't find file {0}".format(fname)

content = [x.strip() for x in content]
f.close()
num = eval(content[0])
P = Polynomial_ring(num, 10)
content = content[1:]
equations = [eval(re.sub('\^', '**', p)) for p in content]
groebner = Groebner(equations)
time_groebner = tm.default_timer()
print('\ngroebner basis generated in {0}'.format(
    round(time_groebner - start_time, 3)))

d_v = generate_dv(exp)
time_allmonomials = tm.default_timer()
print('all monomials generated in {0}'.format(
    round(time_allmonomials - start_time, 3)))

x_var = generate_var_x(num)
dic, value, key, coeff, std_monomial = standard_monomial()

if len(dic) == 0:
    print('w = 0')
time_std_monomial = tm.default_timer()
print('standard monomial calculated in {0}\n'.format(
    round(time_std_monomial - start_time, 3)))


print('standard monomial set is {0}' .format(std_monomial))

#############################################################
##################### Initial condition #####################
#############################################################

print "\nGive initial condition of pde in terms of",
init = initial_condition_var(num, value)
for i in range(len(init)):
    print init[i],
print_initial_condition(std_monomial)
i = raw_input("\nEnter the initialcondition: ")

f(x) = i
print '\nevaluation complete'
monomial_coeff = standard_monomial_evaluated_over_initial_condition(
    f(x), std_monomial)

solution, solution_approx = generate_solution(num,monomial_coeff, f(x), 20, x_var)


##############################################################
#######################  OUTPUT  #############################
##############################################################
print("\nThe solution of the pde's:")
print(solution)

filepath = 'pde_output.txt'
f = open(filepath, 'w')
f.write(str(solution))
f.close()
