import numpy as np
import compute_coefficients as compcoeff


coefficients_1 = compcoeff.GetCoefficients_AllItems(1)
coefficients_3 = compcoeff.GetCoefficients_AllItems(3)
coefficients_4 = compcoeff.GetCoefficients_AllItems(4)
coefficients_15 = compcoeff.GetCoefficients_AllItems(15)


def computeDistance(dimension, points, center):
    sum = 0
    for i in range(dimension):
        sum += (points[i] - center[i]) **2
    return sum


def f1(x):
    return np.polyval(coefficients_1, x)

def f3(x):
    return np.polyval(coefficients_3, x)

def f4(x):
    return np.polyval(coefficients_4, x)

def f15(x):
    return np.polyval(coefficients_15, x)


def sign_4(x):
    for i in range(4):
        x = f4(x)
    return x

def sign_15(x):
    for i in range(3):
        x = f15(x)
    return x

def sign_customized(x):
    for i in range(2):
        x = f3(x)
    for i in range(2):
        x = f15(x)
    return x


def compare_4(x, y, scalar):  #哪个小，对应的bool就是1
    diff = x - y
    diff = diff / scalar
    bool1 = sign_4(diff)
    bool1 = (bool1 + 1.0) / 2.0
    bool0 = 1.0 - bool1
    return bool0, bool1

def compare_15(x, y, scalar):  #哪个小，对应的bool就是1
    diff = x - y
    diff = diff / scalar
    bool1 = sign_15(diff)
    bool1 = (bool1 + 1.0) / 2.0
    bool0 = 1.0 - bool1
    return bool0, bool1

def compare_customized(x, y, scalar):  #哪个小，对应的bool就是1
    diff = x - y
    diff = diff / scalar
    bool1 = sign_customized(diff)
    bool1 = (bool1 + 1.0) / 2.0
    bool0 = 1.0 - bool1
    return bool0, bool1


def compare_truly(x, y, scalar): # 这里其实不需要scalar，但是为了保持接口一致
    if x < y:
        return 1.0, 0.0
    elif x == y :
        return 0.5, 0.5
    else:
        return 0.0, 1.0


def compareDistance_chained(num_centers, distance_matrix, scalar):
    bool0, bool1 = compare_4(distance_matrix[0], distance_matrix[1], scalar)
    bool_matrix = np.array([bool0, bool1])  #密文情况下，新建一个切片
    newdis = bool0 * distance_matrix[0] + bool1 * distance_matrix[1] #用函数接口来实现
    
    for i in range(1, num_centers-1):    
        bool0, bool1 = compare_4(newdis, distance_matrix[i+1], scalar)
        bool_matrix = bool_matrix * bool0 #在密文情况下要写一个循环，可以写成一个函数接口，注意bootstrap
        bool_matrix = np.append(bool_matrix, [bool1], axis=0)  #这是往切片中新加元素
        newdis = bool0 * newdis + bool1 * distance_matrix[i+1]  #注意newdis的level
    return bool_matrix


def compareDistance_arg_onePoint(distances, num_centers, scalar):
    # 可以开一个矩阵来存放，对角线全置为1
    argMatrix = np.ones([num_centers, num_centers])
    for i in range(num_centers):
        for j in range(i+1, num_centers):
            bool0, bool1 = compare_customized(distances[i], distances[j], scalar)
            argMatrix[i,j] = bool0
            argMatrix[j,i] = bool1

    bool_matrix = np.zeros([num_centers])
    for i in range(num_centers):
        bool = argMatrix[i,0]
        for j in range(1, num_centers):
            bool = bool * argMatrix[i,j]
        bool_matrix[i] = bool

    return bool_matrix


def compareDistance_arg_manyPoints(distance_matrix, num_centers, num_points, scalar): #这种写法很耦合，不好。应该忽视num_points。如果要计算很多个点的话，在最外面套一个关于num_points的循环就好了。 
    # 可以开一个矩阵来存放，对角线全置为1
    argMatrix = np.ones([num_centers, num_centers, num_points])
    for i in range(num_centers):
        for j in range(i+1, num_centers):
            bool0, bool1 = compare_4(distance_matrix[i], distance_matrix[j], scalar)
            argMatrix[i,j,:] = bool0
            argMatrix[j,i,:] = bool1

    bool_matrix = np.zeros([num_centers, num_points])
    for i in range(num_centers):
        bool = argMatrix[i,0,:]
        for j in range(1, num_centers):
            bool = bool * argMatrix[i,j,:]
        bool_matrix[i,:] = bool

    return bool_matrix


def getNewCenter(points, center, bool, dimension):
    newCenter = np.zeros(dimension)
    for i in range(dimension):
        newCenter[i] = np.mean(points[i] * bool + center[i] * (1-bool))
    return newCenter