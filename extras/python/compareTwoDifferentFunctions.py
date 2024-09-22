import numpy as np
import matplotlib.pyplot as plt
import compute_coefficients as compcoeff
coefficients_3 = compcoeff.GetCoefficients_AllItems(3)
coefficients_15 = compcoeff.GetCoefficients_AllItems(15)

# 定义函数
def function_atom(x):
    return -1.5*x*(x-1)*(x+1)

def function_their(x, n):
    for i in range(n):
        x = function_atom(x)
    return x*1.73205080757

def function_my(x):
    # for i in range(2):
    #     x = np.polyval(coefficients_3, x)
    for i in range(2):
        x = np.polyval(coefficients_15, x)
    return x

tobetest = 0.001
print(function_their(tobetest, 5))
print(function_my(tobetest))



x_values = np.linspace(-1, 1, 1000)
y_values_their = function_their(x_values, 5)
y_values_my = function_my(x_values)

# 绘制图像
plt.plot(x_values, y_values_their, label='thier')
plt.plot(x_values, y_values_my, label='my')
plt.title('Function Plot')
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.grid(True)
plt.show()
