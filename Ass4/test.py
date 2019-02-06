# from sympy.solvers import solve
# from sympy import *
# from sympy.core.symbol import symbols
# from sympy.solvers.solveset import nonlinsolve
# from sympy import Symbol
import itertools

this = ["".join(seq) for seq in itertools.product("01", repeat=9)]
permutations = len(this)
print(permutations)
print(this)








# x1,x2 = symbols('x1 x2', real=True)
# f1 = 4*x1**2 - 42*x1 + 4*x1*x2 + 2*x2**2 - 14
# f2 = 2*x1**2 + 4*x1*x2 + 4*x2**3 - 26*x2 - 22
#
# x,y = symbols('x y', real=True)
# g = x**4 - 22*x**2 + x + 114
# g_solved = solve(g)
# print(g_solved)
#
#
#
# # print(f1, f2)
# #
# # solved = nonlinsolve([f1, f2], [x1,x2])
# # print(solved)
# #
# # f3 = f1-f2
# # print(f3)
# # solution = solve(f3, x1)[1]
# #
# # subbed = f3.subs(x1,solution)
# # print(subbed, type(subbed))
# #
# # def f(func, x2):
# #     H = func.subs(f)
# #
# # solved_fin = solve(subbed)
# # print(solved_fin)