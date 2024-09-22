package utilities

import (
	"fmt"
	"math"
	"math/rand"
)

func TransposeMatrix(matrix [][]float64) [][]float64 {

	rows := len(matrix)
	cols := len(matrix[0])

	transposed := make([][]float64, cols)
	for i := range transposed {
		transposed[i] = make([]float64, rows)
	}

	for i := 0; i < rows; i++ {
		for j := 0; j < cols; j++ {
			transposed[j][i] = matrix[i][j]
		}
	}

	return transposed
}

// half-open interval [interval_start, interval_end)
func GenerateRandomNumbers_GivenInterval(interval_start int, interval_end int, num_need int) (randnum []int) {

	for len(randnum) < num_need {
		number := rand.Intn(interval_end - interval_start) + interval_start

		found := false
		for _, num := range randnum {
			if num == number {
				found = true
				break
			}
		}

		if !found {
			randnum = append(randnum, number)
		}
	}
	return randnum
}	

func GetSizeofRange(values []float64) (size_band float64) {
	
	min := values[0]
	max := values[0]
	for i:=0;i<len(values);i++ {
		if min > values[i] {
			min = values[i]
		}
		if max < values[i] {
			max = values[i]
		}
	}
	return max - min
}

func CountUniqueValues(values []int) (count int) {
	
	values_unique := make([]int, 0)
	values_unique = append(values_unique, values[0])
	for i:=1;i<len(values);i++ {
		found := false
		for j:=0;j<len(values_unique);j++ {
			if values[i] == values_unique[j] {
				found = true
				break
			}
		}
		if !found {
			values_unique = append(values_unique, values[i])
		}
	}
	return len(values_unique)
}

// Divide by the specified scale.
func Rescale(matrix [][]float64, scale float64) ([][]float64) {
	
	for i:=0;i<len(matrix);i++ {
		for j:=0;j<len(matrix[0]);j++ {
			matrix[i][j] = matrix[i][j] / scale
		}
	}
	return matrix
}

// Control the printing precision of a float64 slice by first converting it into a string slice.
func FormatFloat(slice_float64 []float64, precision int) (slice_string []string) {
	
	slice_string = make([]string, len(slice_float64))
	for i:=0;i<len(slice_float64);i++ {
		slice_string[i] = fmt.Sprintf("%.*f", precision, slice_float64[i])
	}
    return slice_string
}

// half-open interval [index_start, index_end)
func SelectContinuousColumns(matrix [][]float64, index_start int, index_end int) [][]float64 {

	matrix = TransposeMatrix(matrix)
	return TransposeMatrix(matrix[index_start : index_end])
}

func SelectColumn_GivenIndex(matrix [][]float64, index int) []float64 {

	matrix = TransposeMatrix(matrix)
	return matrix[index]
}

// Calculate an upper bound of the maximum distance between all the samples in a matrix, each coloum representing a sample, each row representing a dimension of all samples.
func CalculateMaximumDistance_UpperBound(matrix [][]float64) (maximumDistance float64) {
	
	maximumDistance = 0.0
	for i := 0; i < len(matrix); i++ {
		maximumDistance += GetSizeofRange(matrix[i]) * GetSizeofRange(matrix[i])
	}
	return math.Sqrt(maximumDistance)
}