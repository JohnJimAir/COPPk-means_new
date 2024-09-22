package main

import (
	"COPPk-means/utilities"
	"encoding/csv"
	"fmt"
	"os"
	"strconv"

	"github.com/sjwhitworth/golearn/pca"
	"gonum.org/v1/gonum/mat"
)

func main() {

    df := ReadMNIST("../../data/MNIST_CSV/mnist_train.csv")
    cx := RemoveFirstColumn(df)
    cv := FlattenMatrix(cx)

    X := mat.NewDense(60000, 784, cv)
    dimension_reduced := 128
    pca := pca.NewPCA(dimension_reduced)
    so := pca.FitTransform(X).RawMatrix().Data

    lk := SliceToMatrix(so, dimension_reduced)

    write_path := "../../data/dataset/mnist_train_reduced_" + strconv.Itoa(dimension_reduced) + ".csv"
    
    if err := WriteCSV(write_path, lk); err != nil {
        panic(err)
    }

    fmt.Println("Done.")
}

// There is a first row of label in the input file.
func ReadMNIST(filePath string) ([][]float64) {

	matrix_string, err := utilities.ReadCSV(filePath)
	if err != nil {
		panic(err)
	}
	return utilities.ParseFloat_Matrix(matrix_string)
}

func RemoveFirstColumn(matrix [][]float64) [][]float64 {
    
	numRows := len(matrix)
    if numRows == 0 {
        return nil
    }
    numCols := len(matrix[0])
    
    newMatrix := make([][]float64, numRows)
    for i := range newMatrix {
        newMatrix[i] = make([]float64, numCols-1)
        copy(newMatrix[i], matrix[i][1:])
    }
    
    return newMatrix
}

// Transform a long slice into a matrix, according the specified number of elements per row.
func SliceToMatrix(slice []float64, length_row int) [][]float64 {
    
	numRows := (len(slice) + length_row - 1) / length_row
    matrix := make([][]float64, numRows)
    
    for i := range matrix {
        start := i * length_row
        end := start + length_row
        if end > len(slice) {
            end = len(slice)
        }
        matrix[i] = slice[start:end]
    }
    
    return matrix
}

// Flatten a matrix into a long slice in a row-major manner.
func FlattenMatrix(matrix [][]float64) []float64 {
    
	numRows := len(matrix)
    if numRows == 0 {
        return nil
    }
    numCols := len(matrix[0])
    
    flattened := make([]float64, numRows*numCols)
    idx := 0
    for _, row := range matrix {
        for _, val := range row {
            flattened[idx] = val
            idx++
        }
    }
    
    return flattened
}

func WriteCSV(filename string, data [][]float64) error {
    
	file, err := os.Create(filename)
    if err != nil {
        return err
    }
    defer file.Close()

    writer := csv.NewWriter(file)
    defer writer.Flush()

    for _, row := range data {
        stringRow := make([]string, len(row))
        for i, val := range row {
            stringRow[i] = strconv.FormatFloat(val, 'f', -1, 64)
        }
        if err := writer.Write(stringRow); err != nil {
            return err
        }
    }

    return nil
}