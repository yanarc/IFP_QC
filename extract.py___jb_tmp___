from numpy import *
from sympy import *
import numpy as np
from dwave_qbsolv import QBSolv


def findmatch(exprr, pattern):
    return [ee.match(pattern) for ee in exprr.find(pattern)]


a = Wild("a", exclude=[Pow])
b = Wild("b", exclude=[1, Pow])
c = Wild("c", exclude=[1, Pow])
d = Wild("d", exclude=[1, Pow])
s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12 = symbols(
        's1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12')
s = [s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12]
f1 = -84*s1*s10 + 2*s1*s2 + 79.0*s1*s3 + 47.5*s1*s4 - 2*s1*s5 - 4*s1*s6 - 8*s1*s7 - 16*s1*s8 - 148*s1*s9 + 130.5*s1 + 6*s10*s2 - 84*s10*s4 + s10*s5 + 2*s10*s6 - 4*s10*s7 - 8*s10*s8 - 81*s10 - 124*s11*s2 - 124*s11*s4 + 2*s11*s5 + 4*s11*s6 - 8*s11*s7 - 16*s11*s8 + s11*s9 - 107*s11 - 84*s12*s2 - 84*s12*s3 + 6*s12*s4 + s12*s5 + 2*s12*s6 - 4*s12*s7 - 8*s12*s8 - 81*s12 + 47.5*s2*s3 + 71.0*s2*s4 - 8*s2*s5 - 16*s2*s6 + s2*s7 + 2*s2*s8 + 6*s2*s9 + 107.5*s2 + 2*s3*s4 - 2*s3*s5 - 4*s3*s6 - 8*s3*s7 - 16*s3*s8 - 148*s3*s9 + 130.5*s3 - 8*s4*s5 - 16*s4*s6 + s4*s7 + 2*s4*s8 + 6*s4*s9 + 107.5*s4 + 34*s5*s6 - 4*s5*s7 - 8*s5*s8 - 8*s5*s9 - 41.0*s5 - 8*s6*s7 - 16*s6*s8 - 16*s6*s9 - 82.0*s6 + 34*s7*s8 + 3.0*s7 + 6.0*s8 - 137*s9 + 808.0
# 以下的为替换的内容
dic = findmatch(f1, a * b * c)
h_list = np.zeros(12)
J_list = np.zeros((12, 12))
h = {}
J = {}
print(dic)
for i in range(len(dic)):
    if dic[i][b] in s:
        if s.index(dic[i][b]) < s.index(dic[i][c]):
            J_list[s.index(dic[i][b])][s.index(dic[i][c])] = dic[i][a]
        else:
            J_list[s.index(dic[i][c])][s.index(dic[i][b])] = dic[i][a]
    else:
        h_list[s.index(dic[i][c])] = dic[i][b]

for i in range(len(h_list)):
    h.setdefault(i, h_list[i])

for i in range(12):
    for j in range(i+1, 12):
        J.setdefault((i, j), J_list[i][j])
print(h)
print(J)