package main

import (
	"COPPk-means/plaintext"
	"COPPk-means/plaintext/check"
	"COPPk-means/utilities"
	"fmt"
	"strconv"
)

func main(){

	precision_print := 8

	num_points_whole := 60000
	num_points := 3158
	num_centers := 10
	dimension := 512

	file_path_label := "../../../z_dataset/mnist_train.csv"
	tmp := utilities.SelectColumn_GivenIndex(ReadMNIST(file_path_label), 0)
	labelSeq_benchmark := make([]int, num_points_whole)
	for i:=0;i<num_points_whole;i++ {
		labelSeq_benchmark[i] = int(tmp[i])+1
	}


	file_path_reduced := "../../../z_dataset/mnist_train_reduced_" + strconv.Itoa(dimension) + ".csv"
	points_whole, maximumDistance := utilities.ReadMNIST_Reduced(file_path_reduced)
	points_whole = utilities.Rescale(points_whole, maximumDistance)
	// index_random := utilities.GenerateRandomNumbers_GivenInterval(0, 60000, 2)
	index_random := []int{52078, 16909}
	points_allbatches := utilities.TransposeMatrix(points_whole)
	points_allbatches = append(points_allbatches, points_allbatches[index_random[0]])
	points_allbatches = append(points_allbatches, points_allbatches[index_random[1]])
	points_allbatches = utilities.TransposeMatrix(points_allbatches)

	num_batch := 19
	point_seperated := make([][][]float64, num_batch)
	for i:=0;i<num_batch;i++ {
		point_seperated[i] = utilities.SelectContinuousColumns(points_allbatches, i*num_points, (i+1)*num_points)
	}

	index_centers := make([]int, num_centers)
	for i := 0; i < num_centers; i++ {
		index_centers[i] = i
	}
	centers := plaintext.ExtractCentersFromPoints(utilities.SelectContinuousColumns(point_seperated[0], 0, num_points), index_centers)


	centers_notstabilized := centers
	centers_stabilized := centers
	var bool_matrix_notstabilized [][]float64
	var bool_matrix_stabilized [][]float64

	// begin computation
	
	iterations_updateCenters := 19*8
	for i := 0; i < iterations_updateCenters; i++{
		fmt.Printf("====== NOW AT THE %dth ITERATION ======\n", i)

		points := point_seperated[i%num_batch]
		// fmt.Println("== NOT STABILIZED ==")
		distances_notstabilized := plaintext.Compute_Distances(points, centers_notstabilized, dimension, num_points, num_centers)

		bool_matrix_notstabilized = plaintext.Compare_Distances(distances_notstabilized, num_points, num_centers, "Exactly")

		centers_notstabilized = plaintext.Update_Centers_NotStabilized(points, bool_matrix_notstabilized, num_points, num_centers, dimension)
		// fmt.Println(centers_notstabilized)
		
		// fmt.Println("== STABILIZED ==")
		distances_stabilized := plaintext.Compute_Distances(points, centers_stabilized, dimension, num_points, num_centers)

		bool_matrix_stabilized = plaintext.Compare_Distances(distances_stabilized, num_points, num_centers, "NotExactly")

		centers_stabilized = plaintext.Update_Centers_Stabilized(points, bool_matrix_stabilized, centers_stabilized, num_points, num_centers, dimension)
		
		// fmt.Println(utilities.FormatFloat(centers_notstabilized[0], 10))
		fmt.Println(utilities.FormatFloat(centers_stabilized[0], precision_print))
		// fmt.Println(utilities.FormatFloat(centers_notstabilized[1], 10))
		fmt.Println(utilities.FormatFloat(centers_stabilized[1], precision_print))
		// fmt.Println(utilities.FormatFloat(centers_notstabilized[63], 10))
		fmt.Println(utilities.FormatFloat(centers_stabilized[63], precision_print))	
	}

	// have obtained the centroids, now compute accuracy
	distances_notstabilized := plaintext.Compute_Distances(points_whole, centers_notstabilized, dimension, num_points_whole, num_centers)
	bool_matrix_notstabilized_whole := plaintext.Compare_Distances(distances_notstabilized, num_points_whole, num_centers, "Exactly")
	labelSeq_notstabilized := check.Transform_BoolMatrix_to_NumLabels(bool_matrix_notstabilized_whole)

	distances_stabilized := plaintext.Compute_Distances(points_whole, centers_stabilized, dimension, num_points_whole, num_centers)
	bool_matrix_stabilized_whole := plaintext.Compare_Distances(distances_stabilized, num_points_whole, num_centers, "NotExactly")
	labelSeq_stabilized := check.Transform_BoolMatrix_to_NumLabels(bool_matrix_stabilized_whole)

	fmt.Println(check.CalculateCorrelation_AllPermutations(labelSeq_benchmark, labelSeq_notstabilized, num_centers))
	fmt.Println(check.CalculateCorrelation_AllPermutations(labelSeq_benchmark, labelSeq_stabilized, num_centers))
	fmt.Println(check.CalculateCorrelation_AllPermutations(labelSeq_notstabilized, labelSeq_stabilized, num_centers))


	
}

// There is a first row of label in the input file.
func ReadMNIST(filePath string) ([][]float64) {

	matrix_string, err := utilities.ReadCSV(filePath)
	if err != nil {
		panic(err)
	}
	return utilities.ParseFloat_Matrix(matrix_string)
}