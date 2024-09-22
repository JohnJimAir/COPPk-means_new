import numpy as np
import math
from fractions import Fraction as Frac
# import matplotlib.pyplot as plt


def custom_formatter(x):
    return "{:.9f}".format(x)
np.set_printoptions(formatter={'float': custom_formatter}, linewidth=1000)


def GetCoefficients_Powerof_MinusXof2(power):
    # n 最小为 0
    coefficients = np.zeros(2*power+1, dtype=object)
    # 先反着计算，因为 n+1 的情况只比 n 多了两个元素在前面，后面都是一样的
    for i in range(2*power+1):
        if i%4 == 0:
            coefficients[i] = 1
        elif i%4 == 2:
            coefficients[i] = -1
    # 最后再倒一下顺序
    coefficients = coefficients[::-1]
    return coefficients

def GetCoefficients_Powerof_OneMinusXof2(power):
    coefficients = GetCoefficients_Powerof_MinusXof2(power)
    for i in range(2*power+1):
        if i%2 == 0:
            coefficients[i] = math.comb(power, i//2) * coefficients[i]
    return coefficients

def GetCoefficients_SingleItem(i_asinpaper):
    coefficients = GetCoefficients_Powerof_OneMinusXof2(i_asinpaper)
    coefficients = np.append(coefficients, 0)
    coefficients = math.comb(2*i_asinpaper, i_asinpaper) * coefficients
    coefficients = Frac(1, 4**i_asinpaper) * coefficients
    return coefficients

def GetCoefficients_AllItems(n_asinpaper):
    coefficients = GetCoefficients_SingleItem(n_asinpaper)
    for i in range(n_asinpaper-1,-1,-1):
        coefficients_tmp = GetCoefficients_SingleItem(i)
        # 在前面填充 0，对齐长度
        coefficients_tmp = np.pad(coefficients_tmp, (2*(n_asinpaper-i), 0), 'constant', constant_values=0)
        coefficients = coefficients + coefficients_tmp
    return coefficients
    

def polynomial(x, coefficients):
    return np.polyval(coefficients, x)





coefficients1 = GetCoefficients_AllItems(31)
# coefficients2 = GetCoefficients_AllItems(1)



x_values = np.linspace(-1.0, 1.0, 20) 
y_values = polynomial(x_values, coefficients1)
# y_values = polynomial(y_values, coefficients2)
# print(y_values)



# level: 3 3,   n: 3 3
# [-1.0 -0.9999999999997824 -0.9999999917255403 -0.9999969741773597 -0.9998414783108927 -0.9972525571127747 -0.9772568735268367 -0.8908617940103182 -0.6588821472937999 -0.24786921765925937 0.2478692176592589 0.6588821472937999 0.8908617940103181 0.9772568735268367 0.9972525571127746 0.9998414783108927 0.9999969741773597 0.99999999172554 0.9999999999997825 1.0]
# level: 2 4,   n: 1 7
# [-1.0 -0.9999999999997903 -0.9999999911777542 -0.9999966269645878 -0.9998222527089993 -0.9969700524760128 -0.9755758972885566 -0.8859978203875248 -0.6518393163126185 -0.2442731316613606 0.24427313166136008 0.6518393163126185 0.8859978203875246 0.9755758972885569 0.9969700524760129 0.9998222527089993 0.9999966269645878 0.9999999911777542 0.9999999999997908 1.0]
# level: 4 2,   n: 7 1
# [-1.0 -0.999999999999559 -0.9999999856500602 -0.999995474022913 -0.9997913112386136 -0.9967295783816947 -0.9748334934653073 -0.8850067055780768 -0.6513535171656496 -0.244246367111853 0.2442463671118525 0.6513535171656496 0.8850067055780765 0.9748334934653072 0.9967295783816947 0.9997913112386136 0.999995474022913 0.9999999856500601 0.999999999999559 1.0]
# level: 5,     n: 7-15,    15
# [-1.0 -0.9999999999990552 -0.9999999713107409 -0.9999918722759652 -0.9996657647505807 -0.9953415404823583 -0.9679630706095451 -0.867344759020655 -0.6275278368756653 -0.23247819595750835 0.23247819595750788 0.6275278368756653 0.8673447590206549 0.9679630706095452 0.9953415404823586 0.9996657647505807 0.9999918722759642 0.9999999713107508 0.999999999999038 1.0]
# level: 6,     n: 17-31,   20
# [-1.0 -0.9999999999996889 -0.9999999998088157 -0.9999996954776397 -0.9999616242566645 -0.9988369359585872 -0.9861340319679678 -0.9156357386651809 -0.6939958563717128 -0.26561864235132643 0.2656186423513259 0.6939958563717128 0.9156357386651809 0.9861340319679678 0.9988369359585866 0.9999616242566645 0.9999996954776492 0.9999999998088532 1.0000000000002698 1.0]
# level: 6,     n: 17-31,   31
# [-1,0 -0.9999999997663357 -0.9999999999665047 -0.9999999997591822 -0.999999646143093 -0.9999411032293966 -0.9976632619238224 -0.9672260270112554 -0.7945563230824538 -0.3252971569147774 0.32529715691477673 0.7945563230824538 0.9672260270112553 0.9976632619238225 0.9999411032293951 0.999999646143093 0.9999999997590399 0.9999999999848097 0.9999999996376331 0.9999999928902401]

# 画出多项式的图像
# plt.plot(x_values, y_values, label='2x^2 - x + 3')
# plt.axhline(0, color='black',linewidth=0.5)
# plt.axvline(0, color='black',linewidth=0.5)
# plt.grid(color = 'gray', linestyle = '--', linewidth = 0.5)
# plt.title('Polynomial Plot')
# plt.xlabel('x')
# plt.ylabel('f(x)')
# plt.legend()
# plt.show()

