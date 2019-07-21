from numpy import *
from sympy import *
import numpy as np
from dwave_qbsolv import QBSolv

init_printing(use_unicode=True)
print('要分解的值为143:')


# 第m块的项目和
def prodfunc(m):
    prodf = 0
    if m == 1:
        for r in range(2, 6):
            prodf = prodf + cols1[r, 1] * 2 ** 0
        for j in range(2, 6):
            prodf = prodf + cols1[j, 0] * 2 ** 1
    elif m == 2:
        for r in range(2, 6):
            prodf = prodf + cols2[r, 1] * 2 ** 0
        for j in range(2, 6):
            prodf = prodf + cols2[j, 0] * 2 ** 1
    elif m == 3:
        for r in range(2, 6):
            prodf = prodf + cols3[r, 2] * 2 ** 0
        for j in range(2, 6):
            prodf = prodf + cols3[j, 1] * 2 ** 1
    else:
        prodf = 0
    return prodf


# 第m-1块给第m块的进位
def cinvalfunc(m):
    if m == 2:
        cinval = cols2[6, 1] * 1 + cols2[6, 0] * 2
    elif m == 3:
        cinval = cols3[6, 2] * 1 + cols3[6, 1] * 2
    else:
        cinval = 0
    return cinval


# 第m块给第m+1块的进位
def max_carryfunc(m):
    if m == 1:
        max_carry = cols2[6, 1] * 4 + cols2[6, 0] * 8
    elif m == 2:
        max_carry = cols3[6, 2] * 4 + cols3[6, 1] * 8
    else:
        max_carry = 0
    return max_carry


# 第m块的目标值
def targetvaluefunc(m):
    if m == 0:
        targetvalue = 1
    elif m == 1:
        targetvalue = cols1[7, 1] * 1 + cols1[7, 0] * 2
    elif m == 2:
        targetvalue = cols2[7, 1] * 1 + cols2[7, 0] * 2
    elif m == 3:
        targetvalue = cols3[7, 2] * 1 + cols3[7, 1] * 2 + cols3[7, 0] * 4
    else:
        targetvalue = 0
    return targetvalue


# 定义匹配函数
def findmatch(exprr, pattern):
    return [ee.match(pattern) for ee in exprr.find(pattern)]


# 主函数
if __name__ == '__main__':
    # 对未知数定义为symbol对象
    p1, p2, q1, q2, c1, c2, c3, c4, t = symbols('p1, p2, q1, q2, c1, c2, c3, c4, t')
    n = list(bin(143))
    binary_n = n[2:]  # 输出二进制的n的数组形式

    # 对n的二进制字符串转为整型
    target_values = []
    for i in range(0, len(binary_n)):
        target_values = np.append(target_values, int(binary_n[i]))  # 字符串转为整型

    #  对乘法表每一行进行输入，这里还可以改进
    p = [0, 0, 0, 0, 1, p2, p1, 1] # 乘法表第1行
    q = [0, 0, 0, 0, 1, q2, q1, 1]  # 乘法表第2行
    product_terms1 = [0, 0, 0, 0, 1, p2, p1, 1]  # 乘法表第3行
    product_terms2 = [0, 0, 0, q1, p2 * q1, p1 * q1, q1, 0]  # 乘法表第4行
    product_terms3 = [0, 0, q2, p2 * q2, p1 * q2, q2, 0, 0]  # 乘法表第5行
    product_terms4 = [0, 1, p2, p1, 1, 0, 0, 0]  # 乘法表第6行
    carries = [0, c2, c1, 0, 0, 0, 0, 0]  # 乘法表第7行

    # 建立乘法表数组
    # AA = [p, q, product_terms1, product_terms2, product_terms3, product_terms4, carries, target_values]
    A1 = np.vstack((p, q, product_terms1, product_terms2, product_terms3, product_terms4, carries, target_values))
    print("二进制乘法表：")
    print(A1)
    # 对整个乘法表进行分块
    left0, cols0 = np.split(A1, [7], axis=1)  # 分出第一块
    left1, cols1 = np.split(left0, [5], axis=1)  # 分出第二块
    cols3, cols2 = np.split(left1, [3], axis=1)  # 分出第三块

    # 建立消耗函数
    f1 = 0
    for i in range(1, 4):
        f1 = f1 + (prodfunc(i) + cinvalfunc(i) - max_carryfunc(i) - targetvaluefunc(i)) ** 2

    # x的平方等于本身
    x = [p1, p2, q1, q2, c1, c2]
    f1 = expand(f1)
    for i in range(0, len(x)):
        f1 = f1.replace(x[i] ** 2, x[i])
    # print(f1)
    # 定义匹配模式对象
    a = Wild("a", exclude=[1, Pow])
    b = Wild("b", exclude=[1, Pow])
    c = Wild("c", exclude=[1, Pow])
    d = Wild("d", exclude=[1, Pow])
    # 多项转化为2项
    f1 = f1.subs(p1, 1 - q1)
    f1 = f1.subs(p2, 1 - q2)
    # print(f1)
    f1 = expand(f1)
    for i in range(0, len(x)):
        f1 = f1.replace(x[i] ** 2, x[i])
    # print(f1)
    # print(f1.find(a * b * c * d))
    f1 = f1.subs(16 * c1 * q1 * q2, 16 * (t * c1 + 2 * (q1 * q2 - 2 * q1 * t - 2 * q2 * t + 3 * t)))
    f1 = f1.subs(32 * c2 * q1 * q2, 32 * (t * c2 + 2 * (q1 * q2 - 2 * q1 * t - 2 * q2 * t + 3 * t)))
    for i in range(0, len(x)):
        f1 = f1.replace(x[i] ** 2, x[i])
    # print(f1)
    #  转化为粒子自旋形式
    x1 = [c1, c2, q1, q2, t]
    s1, s2, s3, s4, s5 = symbols(
        's1, s2, s3, s4, s5')  # 定义二进制变量
    s = [s1, s2, s3, s4, s5]
    for i in range(0, len(x1)):
        f1 = f1.subs(x1[i], (1 - s[i]) / 2)
    f1 = 2*expand(f1)
    print("ISing 模型多项式：")
    print(f1)

    #  提取出h
    print("提取出来的h")
    f2 = f1.replace(a*b*c, 0)
    h_list = np.zeros(5)
    f2 = findmatch(f2, a * b)
    for i in range(len(f2)):
        h_list[s.index(f2[i][b])] = f2[i][a]
    h = {}
    for i in range(len(h_list)):
        h.setdefault(i, h_list[i])
    print(h)

    # 提取出J
    print("提取出来的J")
    J_list = np.zeros((5, 5))
    f3 = findmatch(f1, a * b * c)
    for i in range(len(f3)):
        J_list[s.index(f3[i][b])][s.index(f3[i][c])] = f3[i][a]
    J = {}
    for i in range(5):
        for j in range(i+1, 5):
            J.setdefault((i, j), J_list[i][j])
    print(J)

    # 调用QBsolve第三方库
    response = QBSolv().sample_ising(h, J)
    print("粒子自旋态：")
    spin = list(response.samples())
    print(spin)
    print("粒子自旋态对应能量值：")
    energy = list(response.data_vectors['energy'])
    print(energy)

    # 反推q的值
    q1 = (1-spin[0][2])/2
    q2 = (1-spin[0][3])/2
    q = 1*2**3 + q2 * 2**2 + q1 * 2 + 1
    print('q的值为', end=" ")
    print(q, end=" ")
    q1 = (1 - spin[1][2]) / 2
    q2 = (1 - spin[1][3]) / 2
    q = 1 * 2 ** 3 + q2 * 2 ** 2 + q1 * 2 + 1
    print('或者', end=" ")
    print(q)

