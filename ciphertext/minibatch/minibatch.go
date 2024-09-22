package minibatch

import (
	"goisbest/ciphertext/elaborate"

	"github.com/tuneinsight/lattigo/v5/core/rlwe"
	"github.com/tuneinsight/lattigo/v5/he/hefloat"
)

func EncryptPoints_OneBatch(params hefloat.Parameters, encoder *hefloat.Encoder, encryptor *rlwe.Encryptor, points [][]float64, dimension int, num_points int) (points_ct []*rlwe.Ciphertext){

	capacity, quantity := elaborate.CalCapacityAndQuantity(params.MaxSlots(), num_points, dimension)
	points_ct = make([]*rlwe.Ciphertext, quantity)

	values := make([][]float64, quantity)
	for i:=0;i<dimension;i++ {	
		indexSplit := i/capacity
		values[indexSplit] = append(values[indexSplit], points[i]...)
	}
	for i:=0;i<quantity;i++ { 
		values[i] = append(values[i], make([]float64, params.MaxSlots()-len(values[i]))...) // make the length of values[i] equal to params.MaxSlots()
	}

	for i:=0;i<quantity;i++ {
		points_pt := hefloat.NewPlaintext(params, params.MaxLevel())
		if err := encoder.Encode(values[i], points_pt); err != nil {
			panic(err)
		}
		var err error
		if points_ct[i], err = encryptor.EncryptNew(points_pt); err != nil {
			panic(err)
		}
	}

	return points_ct
}

// make sure that num_points*num_centers is smaller than the number of slots in a ciphertext
func ExtractPoints_OneBatch(eval *hefloat.Evaluator, points_ct []*rlwe.Ciphertext, num_points int, dimension int, num_centers int) (points_formated_ct []*rlwe.Ciphertext) {
	// another method, which use rotation directly instead of calling the map function, could be faster?
	points_formated_ct = make([]*rlwe.Ciphertext, dimension)

	var map_location []elaborate.Pair
	for i:=0;i<dimension;i++ {
		for j:=0;j<num_centers;j++ {
			map_location = append(map_location, elaborate.Pair{First: i, Second: i*num_centers + j})
		}
	}

	points_formated_ct = elaborate.Map_MultiIn_MultiOut(eval, points_ct, map_location, num_points)

	return points_formated_ct
}