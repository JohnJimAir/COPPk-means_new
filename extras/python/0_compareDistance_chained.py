#  输入是距离矩阵，输出是bool标记矩阵 
import functions as func
import numpy as np
np.set_printoptions(suppress=True, precision=17)
np.set_printoptions(linewidth=1000)



dis0 = np.array([12.96, 10.24, 7.29, 0.16,  0.36,  1.44,  7.84, 14.44, 19.36, 26.01])
dis1 = np.array([31.36, 27.04, 22.09,  5.76,  1.96,  0.64,  0.64,  3.24,  5.76,  9.61])
dis2 = np.array([36.0,  31.36, 26.01,  7.84,  3.24,  1.44,  0.16,  1.96,  4.0,  7.29])
distance_matrix = np.array([dis0, dis1])
distance_matrix = np.append(distance_matrix, [dis2], axis=0)
scalar = 100.0


bool0, bool1 = func.compare_2(0.64, 23, scalar)
print(bool0, bool1)

bool_matrix = func.compareDistance_chained(3, distance_matrix, scalar)

print(bool_matrix[0])
print(bool_matrix[1])
print(bool_matrix[2])
 