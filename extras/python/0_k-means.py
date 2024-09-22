#  完整的k-means过程，多项式近似比较，stabilized，argmax
import functions as func
import numpy as np

def custom_formatter(x):
    return "{:.9f}".format(x)
np.set_printoptions(formatter={'float': custom_formatter}, linewidth=1000)


num_points = 4096
num_centers = 4
dimension = 16
start_interval = -1.0
end_interval = 1.0
length_interval = end_interval - start_interval
scalar = length_interval**2 * dimension

# 读文件
data_list = []
with open("/Users/johnjim/Applications/z_dataset/RandNumbers_" + str(num_points) + "_" + str(dimension) + "_" + str(start_interval) + "-" + str(end_interval) + "_.txt", "r") as file:
    for line in file:
        values_str = line.strip().split(" ")
        row_data = [float(value) for value in values_str]
        data_list.append(row_data)
numpy_array = np.array(data_list)
points = numpy_array

centers = np.zeros((num_centers, dimension))
for i in range(num_centers):
    centers[i,:] = points[:,i]


for k in range(2):
    print("第 %d 次迭代", k , end='')

    # 计算距离
    distance_matrix = np.zeros([num_centers, num_points])
    for i in range(num_centers):
        distance_matrix[i,:] = func.computeDistance(dimension, points, centers[i])
    print(distance_matrix[0][0])
    # 比较距离
    bool_matrix = np.zeros([num_centers, num_points])
    for i in range(num_points):
        bool_matrix[:,i] = func.compareDistance_arg_onePoint(distance_matrix[:,i], num_centers, scalar)
    # print(bool_matrix[:,0])
    # print(bool_matrix[:,1])
    # print(bool_matrix[:,2])
    # print(bool_matrix[:,3])
    # print(bool_matrix[:,4])
    # print(bool_matrix[:,5])


    # 更新中心点
    for i in range(num_centers):
        centers[i] = func.getNewCenter(points, centers[i], bool_matrix[i], dimension)
    print(centers[0])
    print(centers[1])
    print(centers[2])
    # print(centers[3])
    # print(centers[4])
    # print(centers[5])
    # print(centers[6])
    # print(centers[7])
