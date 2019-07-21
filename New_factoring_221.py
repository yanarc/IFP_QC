from numpy import *
from sympy import *
import numpy as np
from dwave_qbsolv import QBSolv

init_printing(use_unicode=True)
encoding = 'utf-8-sig'
print('要分解的值为221:')


# 第m块的项目和
def prodfunc(m):
    prodf = 0
    for r in range(2, (cols[m].shape[0] - 2)):
        for u in range(cols[m].shape[1]):
            prodf = prodf + cols[m][r, u] * 2 ** (cols[m].shape[1]-1-u)
    return prodf


# 第m-1块给第m块的进位
def cinvalfunc(m):
    cinval = 0
    for r in range(cols[m].shape[1]):  # 第m块列数
        cinval = cinval + cols[m][(cols[m].shape[0]-2), r] * 2**(cols[m].shape[1]-1-r)
    return cinval


# 第m块给第m+1块的进位
def max_carryfunc(m):
    max_carry = 0
    if m < (len(cols)-1):
        for r in range(cols[m+1].shape[1]):
            max_carry = max_carry + cols[m+1][(cols[m].shape[0]-2), r] * 2**(cols[m+1].shape[1] + cols[m].shape[1]-1-r)
    return max_carry


# 第m块的目标值
def targetvaluefunc(m):
    targetvalue = 0
    for r in range(cols[m].shape[1]):
        targetvalue = targetvalue + cols[m][cols[m].shape[0]-1, r] * 2**(cols[m].shape[1]-1-r)
    return targetvalue


# 3项转为2项的函数
def changefunc(bb, cc):
    tt = 0
    if p3*p1 == bb*cc:
        tt = t5
    if p3*p2 == bb*cc:
        tt = t6
    return tt


# 定义匹配函数
def findmatch(exprr, pattern):
    return [ee.match(pattern) for ee in exprr.find(pattern)]


# 主函数
if __name__ == '__main__':
    # 对未知数定义为symbol对象
    p1, p2, p3, p4, p5, p6, q1, q2, q3, q4, q5, q6, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, = symbols(
        'p1, p2, p3, p4, p5, p6, q1, q2, q3, q4, q5, q6, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11')
    t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18 = symbols(
        't1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18')
    t19, t20, t21, t22, t23, t24, t25, t26, t27, t28, t29, t30, t31, t32, t33, t34, t35, t36 = symbols(
        't19, t20, t21, t22, t23, t24, t25, t26, t27, t28, t29, t30, t31, t32, t33, t34, t35, t36')

    x = [p3, q1, q2, c1, c2, c3, c4, t5, t6]
    p_list = [p3]
    q_list = [q1, q2]
    p_q_list = [p3, q1, q2]
    c_list = [c1, c2, c3, c4]
    t_list = [t5, t6]
    n = list(bin(221))
    binary_n = n[2:]  # 输出二进制的n的数组形式

    # 对n的二进制字符串转为整型
    target_values = []
    for i in range(0, len(binary_n)):
        target_values = np.append(target_values, int(binary_n[i]))  # 字符串转为整型

    #  对乘法表每一行进行输入，这里还可以改进
    p = [0, 0, 0, 1, p3, p2, p1, 1]
    q = [0, 0, 0, 0, 1, q2, q1, 1]
    product_terms1 = [0, 0, 0, 1, p3, p2, p1, 1]
    product_terms2 = [0, 0, q1, p3*q1, p2*q1, p1*q1, q1, 0]
    product_terms3 = [0, q2, p3*q2,	p2*q2, p1*q2, q2, 0, 0]
    product_terms4 = [1, p3, p2, p1, 1, 0, 0, 0]
    carries = [c4, c3, c2, c1, 0, 0, 0, 0]

    # 建立乘法表数组
    A1 = np.vstack((p, q, product_terms1, product_terms2, product_terms3, product_terms4, carries, target_values))
    print(A1)

    # 对整个乘法表进行分块
    left0, cols0 = np.split(A1, [7], axis=1)  # 分出第一块
    left1, cols1 = np.split(left0, [4], axis=1)  # 分出第二块
    cols3, cols2 = np.split(left1, [2], axis=1)   # 分出第三块
    cols = [cols0, cols1, cols2, cols3]

    # 建立消耗函数
    f1 = 0
    for i in range(1, len(cols)):
        f1 = f1 + (prodfunc(i) + cinvalfunc(i) - max_carryfunc(i) - targetvaluefunc(i)) ** 2

    # x的平方等于本身
    print('最开始的函数：')
    f1 = expand(f1)
    for i in range(0, len(x)):
        f1 = f1.replace(x[i] ** 2, x[i])
    print(f1)

    # 根据新算法减去某一些比特
    f1 = f1.subs(p1, (1-q1))
    f1 = f1.subs(p2, (1-q2))
    f1 = expand(f1)
    for i in range(0, len(x)):
        f1 = f1.replace(x[i] ** 2, x[i])
        f1 = f1.replace(x[i] ** 3, x[i])
    print(f1)

    # 定义匹配模式对象
    a = Wild("a", exclude=[Pow])
    b = Wild("b", exclude=[1, Pow])
    c = Wild("c", exclude=[1, Pow])
    d = Wild("d", exclude=[1, Pow])
    e = Wild("e", exclude=[1, Pow])
    y = Wild("y", exclude=[1, Pow])

    #  先单独消掉四次项
    change_1 = findmatch(f1, a * b * c * d * e)  # 将3次和4次项提取出来放入到dict中，包括a等于1的情况
    ff = 0
    for i in range(len(change_1)):
        if change_1[i][b] in p_list:  # 判断是否为四次项
            f1 = f1 - change_1[i][a] * change_1[i][b] * change_1[i][c] * change_1[i][d] * change_1[i][e]
            if change_1[i][a] > 0:
                ff = ff + change_1[i][a]*(changefunc(change_1[i][b], change_1[i][d])*change_1[i][c]*change_1[i][e]+2*(change_1[i][b]*
                          change_1[i][d]-2*change_1[i][b]*changefunc(change_1[i][b], change_1[i][d])-2*change_1[i][d]*changefunc(change_1[i][b], change_1[i][d])+3*changefunc(change_1[i][b],
                                                                                                    change_1[i][d])))
            if change_1[i][a] < 0:
                ff = ff + change_1[i][a] * (changefunc(change_1[i][b], change_1[i][d]) * change_1[i][c] * change_1[i][e] - 2 * (change_1[i][b] *
                          change_1[i][d] - 2 *change_1[i][b] * changefunc(change_1[i][b], change_1[i][d]) - 2 *
                                                            change_1[i][d] * changefunc(change_1[i][b], change_1[i][d]) + 3 *changefunc(change_1[i][b], change_1[i][d])))
    f1 = f1 + ff
    # 到这里的时候已经把四项全部转成3项了
    #  将3项转化为2项
    # 消掉三项为cpq的
    change_2 = findmatch(f1, a * b * c * d)
    ff = 0
    for i in range(len(change_2)):
        if change_2[i][b] in c_list:  # 判断是否为三次项
            f1 = f1 - change_2[i][a] * change_2[i][b] * change_2[i][c] * change_2[i][d]
            if change_2[i][a] > 0:  # 判断三次项前面常数正负
                ff = ff + change_2[i][a] * (changefunc(change_2[i][c], change_2[i][d]) * change_2[i][b] + 2 * (
                        change_2[i][c] * change_2[i][d] - 2 * change_2[i][c] * changefunc(change_2[i][c],
                                                                                          change_2[i][d]) -
                        2 * change_2[i][d] * changefunc(change_2[i][c], change_2[i][d]) + 3 * changefunc(
                         change_2[i][c], change_2[i][d])))
            elif change_2[i][a] < 0:  # 判断三次项前面常数正负
                ff = ff + change_2[i][a] * (changefunc(change_2[i][c], change_2[i][d]) * change_2[i][b] - 2 * (
                        change_2[i][c] * change_2[i][d] - 2 * change_2[i][c] * changefunc(change_2[i][c],
                                                                                          change_2[i][d]) -
                        2 * change_2[i][d] * changefunc(change_2[i][c], change_2[i][d]) + 3 * changefunc(
                         change_2[i][c], change_2[i][d])))

    #  消掉3项为ppq,或者pqq的
    for i in range(len(change_2)):
        if change_2[i][b] in p_q_list and change_2[i][d] in q_list:
            f1 = f1 - change_2[i][a] * change_2[i][b] * change_2[i][c] * change_2[i][d]
            if change_2[i][c] in p_list:  # 判断如果第二个未知数是p，那么b和d结合
                if change_2[i][a] > 0:
                    ff = ff + change_2[i][a] * (changefunc(change_2[i][b], change_2[i][d]) * change_2[i][c] + 2 * (
                            change_2[i][b] * change_2[i][d] - 2 * change_2[i][b] * changefunc(change_2[i][b],
                                                                                              change_2[i][d]) -
                            2 * change_2[i][d] * changefunc(change_2[i][b], change_2[i][d]) + 3 * changefunc(
                             change_2[i][b], change_2[i][d])))
                elif change_2[i][a] < 0:
                    ff = ff + change_2[i][a] * (changefunc(change_2[i][b], change_2[i][d]) * change_2[i][c] - 2 * (
                            change_2[i][b] * change_2[i][d] - 2 * change_2[i][b] * changefunc(change_2[i][b],
                                                                                              change_2[i][d]) -
                            2 * change_2[i][d] * changefunc(change_2[i][b], change_2[i][d]) + 3 * changefunc(
                             change_2[i][b], change_2[i][d])))
            if change_2[i][c] in q_list:  # 判断如果第二个未知数是p去，那么b和结合
                if change_2[i][a] > 0:
                    ff = ff + change_2[i][a] * (changefunc(change_2[i][b], change_2[i][c]) * change_2[i][d] + 2 * (
                            change_2[i][b] * change_2[i][c] - 2 * change_2[i][b] * changefunc(change_2[i][b],
                                                                                              change_2[i][c]) -
                            2 * change_2[i][c] * changefunc(change_2[i][b], change_2[i][c]) + 3 * changefunc(
                             change_2[i][b], change_2[i][c])))
                elif change_2[i][a] < 0:
                    ff = ff + change_2[i][a] * (changefunc(change_2[i][b], change_2[i][d]) * change_2[i][c] - 2 * (
                            change_2[i][b] * change_2[i][c] - 2 * change_2[i][b] * changefunc(change_2[i][b],
                                                                                              change_2[i][c]) -
                            2 * change_2[i][c] * changefunc(change_2[i][b], change_2[i][c]) + 3 * changefunc(
                             change_2[i][b], change_2[i][c])))

    # 消掉三项为pqt的
    for i in range(len(change_2)):
        if change_2[i][b] in p_list and change_2[i][d] in t_list:  # 判断是否为三次项pqt
            f1 = f1 - change_2[i][a] * change_2[i][b] * change_2[i][c] * change_2[i][d]
            if change_2[i][a] > 0:  # 判断三次项前面常数正负
                ff = ff + change_2[i][a] * (changefunc(change_2[i][b], change_2[i][c]) * change_2[i][d] + 2 * (
                            change_2[i][b] * change_2[i][c] - 2 * change_2[i][b] * changefunc(change_2[i][b],
                                                                                              change_2[i][c]) -
                            2 * change_2[i][c] * changefunc(change_2[i][b], change_2[i][c]) + 3 * changefunc(
                             change_2[i][b], change_2[i][c])))
            elif change_2[i][a] < 0:  # 判断三次项前面常数正负
                ff = ff + change_2[i][a] * (changefunc(change_2[i][b], change_2[i][c]) * change_2[i][d] - 2 * (
                        change_2[i][b] * change_2[i][c] - 2 * change_2[i][b] * changefunc(change_2[i][b],
                                                                                          change_2[i][c]) -
                        2 * change_2[i][c] * changefunc(change_2[i][b], change_2[i][c]) + 3 * changefunc(
                         change_2[i][b], change_2[i][c])))

    f1 = f1 + ff    # 转换后加上二次项
    print('转换成二次项后的函数：')
    print(f1)   # 到这里为止所有已经转成二次项

    # 转化为粒子自旋形式
    s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12 = symbols(
        's1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12')  # 定义二进制未知变量
    s13, s14, s15, s16, s17, s18, s19, s20, s21, s22, s23, s24 = symbols(
        's13, s14, s15, s16, s17, s18, s19, s20, s21, s22, s23, s24')  # 定义二进制未知变量
    s25, s26, s27, s28, s29, s30, s31, s32, s33, s34, s35, s36 = symbols(
        's25, s26, s27, s28, s29, s30, s31, s32, s33, s34, s35, s36')  # 定义二进制未知变量
    s37, s38, s39, s40, s41, s42, s43, s44, s45, s46, s47, s48 = symbols(
        's37, s38, s39, s40, s41, s42, s43, s44, s45, s46, s47, s48')  # 定义二进制未知变量
    s49, s50, s51, s52, s53, s54, s55, s56, s57, s58, s59 = symbols(
        's49, s50, s51, s52, s53, s54, s55, s56, s57, s58, s59')  # 定义二进制未知变量
    s = [s1, s2, s3, s4, s5, s6, s7, s8, s9]
    for i in range(len(x)):
        f1 = f1.subs(x[i], ((1 - s[i]) / 2))
    f1 = 2*f1
    f1 = expand(f1)
    print("ISing 模型多项式：")
    print(f1)

    #  提取h 和 J
    dic = findmatch(f1, a * b * c)
    h_list = np.zeros(len(x))
    J_list = np.zeros((len(x), len(x)))
    h = {}
    J = {}
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

    for i in range(len(x)):
        for j in range(i + 1, len(x)):
            J.setdefault((i, j), J_list[i][j])
    print('提取出h:')
    print(h)
    print('提取出J：')
    print(J)

    # 调用qbsolve第三方库
    count_all_first = 1  # 输入调用第三方库的次数
    correct_count = 0
    count = 0
    q_truth_1 = 13  # 输入质因子
    q_truth_2 = 17  # 输入质因子
    for i in range(count_all_first):
        response = QBSolv().sample_ising(h, J)
        print("粒子自旋态：")
        spin = list(response.samples())
        print(spin)
        print("粒子自旋态对应能量值：")
        energy = list(response.data_vectors['energy'])
        print(energy)

        # 反推回去q的值
        q = 2 ** (len(q_list) + 1) + 1
        for j in range(len(q_list)):
            q_list[j] = (1 - spin[0][len(p_list) + j]) / 2
        for j in range(len(q_list)):
            q = q + q_list[j] * 2 ** (j + 1)
        print('q的值为')
        print(q)
        if q == q_truth_1 or q == q_truth_2:
            correct_count = correct_count + 1
        else:
            # 如果能量值和第一个一样小，那么
            for k in range(1, len(energy)):
                if energy[k] <= energy[0]:
                    count = count + 1
                    q = 2 ** (len(q_list) + 1) + 1
                    for j in range(len(q_list)):  # 找到对应的q值
                        q_list[j] = (1 - spin[k][len(p_list) + j]) / 2
                    for j in range(len(q_list)):
                        q = q + q_list[j] * 2 ** (j + 1)
                    if q == q_truth_1 or q == q_truth_2:
                        correct_count = correct_count + 1
                    print('或者为', end='')
                    print(q)
    count_all = count_all_first + count
    print('h的最小值为：', end='')
    print(min(h.values()), end='')
    print(', 最大值为：', end='')
    print(max(h.values()))
    print('J的最小值为：', end='')
    print(min(J.values()), end='')
    print(', 最大值为：', end='')
    print(max(J.values()))
    print('总的最低能量获得次数：', end='')
    print(count_all)
    print('正确的次数为：', end='')
    print(correct_count)
    print('正确率为：', end='')
    print(correct_count / count_all)
    print('所用的比特数为：', end='')
    print(len(s))