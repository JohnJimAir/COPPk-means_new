import random


start_interval = 0.0
end_interval = 1.0

num_points = 3276
dimension = 32


random_numbers_2d = [[random.uniform(start_interval, end_interval) for _ in range(num_points)] for _ in range(dimension)]
with open("/Users/johnjim/Applications/z_dataset/RandNumbers_" + str(num_points) + "_" + str(dimension) + "_" + str(start_interval) + "-" + str(end_interval) + "_.txt", "w") as file:
    for row in random_numbers_2d:
        row_str = " ".join(map(str, row))
        file.write(f"{row_str}\n")
print("随机数已保存到文件")