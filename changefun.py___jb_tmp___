from numpy import *
from sympy import *
import numpy as np
from dwave_qbsolv import QBSolv


def findmatch(exprr, pattern):
    return [ee.match(pattern) for ee in exprr.find(pattern)]


def changefunc(bb, cc):
    tt = 0
    if bb*cc == p1*q1:
        tt = t1
    elif bb*cc == p1*q2:
        tt = t2
    elif bb*cc == p2*q1:
        tt = t4
    elif bb*cc == p2*q2:
        tt = t3
    return tt


p1, p2, q1, q2, c1, c2, c3, c4, t1, t2, t3, t4 = symbols('p1, p2, q1, q2, c1, c2, c3, c4, t1, t2, t3, t4')
f1 = 68*c1*c2 - 8*c1*c3 - 16*c1*c4 - 16*c1*p1*q1 + 2*c1*p1*q2 - 4*c1*p1 + 2*c1*p2*q1 + 4*c1*p2*q2 - 16*c1*p2 - 4*c1*q1 - 16*c1*q2 + 43.0*c1 - 16*c2*c3 - 32*c2*c4 - 32*c2*p1*q1 + 4*c2*p1*q2 - 8*c2*p1 + 4*c2*p2*q1 + 8*c2*p2*q2 - 32*c2*p2 - 8*c2*q1 - 32*c2*q2 + 120.0*c2 + 68*c3*c4 - 8*c3*p1*q2 - 16*c3*p1 - 8*c3*p2*q1 - 16*c3*p2*q2 + 2*c3*p2 - 16*c3*q1 + 2*c3*q2 + 5.0*c3 - 16*c4*p1*q2 - 32*c4*p1 - 16*c4*p2*q1 - 32*c4*p2*q2 + 4*c4*p2 - 32*c4*q1 + 4*c4*q2 + 44.0*c4 + 2*p1*p2*q1*q2 + 12*p1*p2*q1 + 12*p1*p2*q2 + 4*p1*p2 + 12*p1*q1*q2 + 10.0*p1*q1 + 11.0*p1*q2 + 3.0*p1 + 12*p2*q1*q2 + 11.0*p2*q1 + 18.0*p2*q2 - 11.0*p2 + 4*q1*q2 + 3.0*q1 - 11.0*q2 + 14.0
f1 = f1.subs(p1 * p2 * q1 * q2, (t1 * t3 + 2 * (p2 * q2 - 2 * p2 * t3 - 2 * q2 * t3 + 3 * t3) + 2 * (
                p1 * q1 - 2 * p1 * t1 - 2 * q1 * t1 + 3 * t1)))

a = Wild("a", exclude=[Pow])
b = Wild("b", exclude=[1, Pow])
c = Wild("c", exclude=[1, Pow])
d = Wild("d", exclude=[1, Pow])

x = [p1, p2, q1, q2, c1, c2, c3, c4, t1, t2, t3, t4]
p_list = [p1, p2]
q_list = [q1, q2]
p_q_list = [p1, p2, q1, q2]
c_list = [c1, c2, c3, c4]
t_list = [t1, t2, t3, t4]
print(findmatch(f1, a * b * c * d))
change_1 = findmatch(f1, a * b * c * d)

ff = 0
# 消掉三项为cpq的
for i in range(len(change_1)):
    if change_1[i][b] in c_list:  # 判断是否为三次项
        f1 = f1 - change_1[i][a]*change_1[i][b]*change_1[i][c]*change_1[i][d]
        if change_1[i][a] > 0:   # 判断三次项前面常数正负
            ff = ff + change_1[i][a] * (changefunc(change_1[i][c], change_1[i][d]) * change_1[i][b] + 2 * (
                        change_1[i][c] * change_1[i][d] - 2 * change_1[i][c] * changefunc(change_1[i][c],
                                                                                          change_1[i][d]) -
                        2 * change_1[i][d] * changefunc(change_1[i][c], change_1[i][d]) + 3 * changefunc(
                         change_1[i][c], change_1[i][d])))
        elif change_1[i][a] < 0:  # 判断三次项前面常数正负
            ff = ff + change_1[i][a] * (changefunc(change_1[i][c], change_1[i][d]) * change_1[i][b] - 2 * (
                        change_1[i][c] * change_1[i][d] - 2 * change_1[i][c] * changefunc(change_1[i][c],
                                                                                          change_1[i][d]) -
                        2 * change_1[i][d] * changefunc(change_1[i][c], change_1[i][d]) + 3 * changefunc(
                         change_1[i][c], change_1[i][d])))


#  消掉3项为ppq,或者pqq的
for i in range(len(change_1)):
    if change_1[i][b] in p_q_list:
        print(change_1[i])
        f1 = f1 - change_1[i][a] * change_1[i][b] * change_1[i][c] * change_1[i][d]
        if change_1[i][c] in p_list: # 判断如果第二个未知数是p，那么b和d结合
            if change_1[i][a] > 0:
                ff = ff + change_1[i][a] * (changefunc(change_1[i][b], change_1[i][d]) * change_1[i][c] + 2 * (
                        change_1[i][b] * change_1[i][d] - 2 * change_1[i][b] * changefunc(change_1[i][b],
                                                                                          change_1[i][d]) -
                        2 * change_1[i][d] * changefunc(change_1[i][b], change_1[i][d]) + 3 * changefunc(
                    change_1[i][b], change_1[i][d])))
            elif change_1[i][a] < 0:
                ff = ff + change_1[i][a] * (changefunc(change_1[i][b], change_1[i][d]) * change_1[i][c] - 2 * (
                        change_1[i][b] * change_1[i][d] - 2 * change_1[i][b] * changefunc(change_1[i][b],
                                                                                          change_1[i][d]) -
                        2 * change_1[i][d] * changefunc(change_1[i][b], change_1[i][d]) + 3 * changefunc(
                    change_1[i][b], change_1[i][d])))
        if change_1[i][c] in q_list:  # 判断如果第二个未知数是p去，那么b和结合
            if change_1[i][a] > 0:
                ff = ff + change_1[i][a] * (changefunc(change_1[i][b], change_1[i][c]) * change_1[i][d] + 2 * (
                        change_1[i][b] * change_1[i][c] - 2 * change_1[i][b] * changefunc(change_1[i][b],
                                                                                          change_1[i][c]) -
                        2 * change_1[i][c] * changefunc(change_1[i][b], change_1[i][c]) + 3 * changefunc(
                    change_1[i][b], change_1[i][c])))
            elif change_1[i][a] < 0:
                ff = ff + change_1[i][a] * (changefunc(change_1[i][b], change_1[i][d]) * change_1[i][c] - 2 * (
                        change_1[i][b] * change_1[i][c] - 2 * change_1[i][b] * changefunc(change_1[i][b],
                                                                                          change_1[i][c]) -
                        2 * change_1[i][c] * changefunc(change_1[i][b], change_1[i][c]) + 3 * changefunc(
                    change_1[i][b], change_1[i][c])))


f1 = f1 + ff
print(f1 - (43*c1 + 120*c2 + 5*c3 + 44*c4 + 3*p1 - 11*p2 + 3*q1 - 11*q2 + 444*t1 + 252*t2 + 372*t3 + 252*t4 + 68*c1*c2 - 8*c1*c3
-16*c1*c4 - 16*c2*c3 - 32*c2*c4 + 68*c3*c4 - 4*c1*p1 - 16*c1*p2 - 8*c2*p1 - 32*c2*p2 - 16*c3*p1 + 2*c3*p2 - 32*c4*p1 + 4*c4*p2
-4*c1*q1 - 16*c1*q2 - 8*c2*q1 - 32*c2*q2 - 16*c3*q1 + 2*c3*q2 - 32*c4*q1 + 4*c4*q2 - 16*c1*t1 + 2*c1*t2 - 32*c2*t1 + 4*c1*t3
+4*c2*t2 + 2*c1*t4 + 8*c2*t3 - 8*c3*t2 + 4*c2*t4 - 16*c3*t3 - 16*c4*t2 - 8*c3*t4 - 32*c4*t3 - 16*c4*t4 + 4*p1*p2 + 158*p1*q1
+95*p1*q2 + 95*p2*q1 + 142*p2*q2 + 4*q1*q2 - 296*p1*t1 - 168*p1*t2 + 12*p2*t1 + 12*p2*t2 - 248*p2*t3 - 168*p2*t4 - 296*q1*t1
+12*q2*t1 - 168*q2*t2 - 168*q1*t4 - 248*q2*t3 + 12*q2*t4 + 2*t1*t3 + 14))



