package utilities

import (
	"bufio"
	"encoding/csv"
	"fmt"
	"os"
	"strconv"
	"strings"
)

func ReadPEGASUS(file_directory string, num_points int, dimension int, interval_start float64, interval_end float64, num_centers int) (points [][]float64, maximumDistance float64) {

	filepath := file_directory + "/RandNumbers_" + strconv.Itoa(num_points) + "_" + strconv.Itoa(dimension) + "_" + strconv.FormatFloat(interval_start, 'f', 1, 64) + "-" + strconv.FormatFloat(interval_end, 'f', 1, 64) + "_.txt"
	points, err := ReadFloatArrayFromFile(filepath)
	if err != nil {
		fmt.Println("Error reading file:", err)
		return
	}

	maximumDistance = CalculateMaximumDistance_UpperBound(points)

	return points, maximumDistance
}

func ReadFCPS(filepath string) (points [][]float64, dimension int, num_points int, num_centers int, maximumDistance float64, labelSeq_benchmark []int) {

	records, err := ReadCSV(filepath)
	if err != nil {
		panic(err)
	}
	data_includeLabel := ParseFloat_Matrix(records[1:])
	data_includeLabel = TransposeMatrix(data_includeLabel)
	
	points = data_includeLabel[1:]
	dimension = len(points)
	num_points = len(points[0])
	maximumDistance = CalculateMaximumDistance_UpperBound(points)

	labelSeq_benchmark = make([]int, num_points)
	for i, value := range data_includeLabel[0] {
		labelSeq_benchmark[i] = int(value)
	}

	num_centers = CountUniqueValues(labelSeq_benchmark)

	return points, dimension, num_points, num_centers, maximumDistance, labelSeq_benchmark
}

func ReadG2(file_directory string, dimension int) (points [][]float64, num_points int, num_centers int, maximumDistance float64, labelSeq_benchmark []int) {

	filepath := file_directory + "/g2-" + strconv.Itoa(dimension) + "-20.txt"
	points, err := ReadFloatArrayFromFile(filepath)
	if err != nil {
		fmt.Println("Error reading file:", err)
		return
	}
	points = TransposeMatrix(points)

	maximumDistance = CalculateMaximumDistance_UpperBound(points)

	num_points = 2048
	num_centers = 2
	labelSeq_benchmark = make([]int, 2048)
	for i := 0; i < 1024; i++ {
		labelSeq_benchmark[i] = 1
		labelSeq_benchmark[i+1024] = 2
	}
	return points, num_points, num_centers, maximumDistance, labelSeq_benchmark
}

// There is not the first row of label in the input file.
func ReadMNIST_Reduced(filePath string) (points [][]float64, maximumDistance float64) {
	
	matrix_string, err := ReadCSV(filePath)
	if err != nil {
		panic(err)
	}
	points = ParseFloat_Matrix(matrix_string)

	points = TransposeMatrix(points)
	maximumDistance = CalculateMaximumDistance_UpperBound(points)

	return points, maximumDistance
}

func ReadFloatArrayFromFile(filename string) ([][]float64, error) {
	var result [][]float64

	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()

	// Set the size of buffer of Scanner explicitly
	scanner := bufio.NewScanner(file)
	const maxScanTokenSize = 64 * 20480 // Set the size of buffer so as to accommodate the largest line
	buf := make([]byte, maxScanTokenSize)
	scanner.Buffer(buf, maxScanTokenSize)

	for scanner.Scan() {
		line := scanner.Text()
		values := strings.Fields(line)
		row := make([]float64, len(values))

		for i, value := range values {
			num, err := strconv.ParseFloat(value, 64)
			if err != nil {
				return nil, err
			}
			row[i] = num
		}
		result = append(result, row)
	}
	if err := scanner.Err(); err != nil {
		return nil, err
	}
	return result, nil
}

func ReadCSV(filePath string) ([][]string, error) {
	
	file, err := os.Open(filePath)
	if err != nil {
		fmt.Println("Error opening file:", err)
		return nil, err
	}
	defer file.Close()

	reader := csv.NewReader(file)
	records, err := reader.ReadAll()
	if err != nil {
		fmt.Println("Error reading file:", err)
		return nil, err
	}

	return records, nil
}

// Parse a string matrix to a float64 matrix.
func ParseFloat_Matrix(matrix_string [][]string) (matrix_float64 [][]float64) {
	
	matrix_float64 = make([][]float64, len(matrix_string))

	for i := 0; i < len(matrix_string); i++ {
		matrix_float64[i] = make([]float64, len(matrix_string[i]))
		for j := 0; j < len(matrix_string[i]); j++ {
			value, err := strconv.ParseFloat(matrix_string[i][j], 64)
			if err != nil {
				return nil
			}
			matrix_float64[i][j] = value
		}
	}

	return matrix_float64
}