package main

import "fmt"

func main() {
	array := []float64{24.809902119, 24.518310939}
	sum := Sum(array)
	fmt.Println(sum)

}
func Sum(array []float64) (result float64) {
	result = 0.0
	for i:=0;i<len(array);i++ {
		result += array[i]
	}
	return result
}