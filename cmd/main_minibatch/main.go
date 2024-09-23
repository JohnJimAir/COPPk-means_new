package main

import (
	"COPPk-means/ciphertext/elaborate"
	"COPPk-means/ciphertext/minibatch"
	"COPPk-means/ciphertext/printers"
	"COPPk-means/plaintext"
	"COPPk-means/utilities"
	"flag"
	"fmt"
	"time"

	"github.com/tuneinsight/lattigo/v5/core/rlwe"
	"github.com/tuneinsight/lattigo/v5/he/hefloat"
	"github.com/tuneinsight/lattigo/v5/he/hefloat/bootstrapping"
	"github.com/tuneinsight/lattigo/v5/ring"
	"github.com/tuneinsight/lattigo/v5/utils"
)

var flagShort = flag.Bool("short", false, "run the example with a smaller and insecure ring degree.")

func main() {

	precision_print := 8
	num_points := 3158
	num_centers := 10
	dimension := 64
	 
	flag.Parse()
	LogN := 16
	if *flagShort {
		LogN -= 3
	}

	startTime := time.Now()
	var params hefloat.Parameters
	var err error
	if params, err = hefloat.NewParametersFromLiteral(
		hefloat.ParametersLiteral{
			LogN:            LogN,                                              // Log2 of the ring degree
			LogQ:            []int{55, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40}, // Log2 of the ciphertext prime moduli
			LogP:            []int{61, 61, 61},                                 // Log2 of the key-switch auxiliary prime moduli
			LogDefaultScale: 40,                                                // Log2 of the scale
			Xs:              ring.Ternary{H: 192},
			RingType:        0,
		}); err != nil {
		panic(err)
	}
	kgen := rlwe.NewKeyGenerator(params)
	sk, pk := kgen.GenKeyPairNew()
	rlk := kgen.GenRelinearizationKeyNew(sk)
	evk := rlwe.NewMemEvaluationKeySet(rlk)

	encoder := hefloat.NewEncoder(params)
	encryptor := rlwe.NewEncryptor(params, pk)
	decryptor := rlwe.NewDecryptor(params, sk)
	eval := hefloat.NewEvaluator(params, evk)

	// LogSlots := params.LogMaxSlots()
	// Slots := 1 << LogSlots

	fmt.Println("== GENERATE ROTATION EVALUATOR ==")
	galEls := []uint64{
		params.GaloisElement(num_points), params.GaloisElement(-num_points), params.GaloisElement(1), params.GaloisElement(-1), params.GaloisElement(num_centers),
	}
	for i := elaborate.GetPowerof2_SmallerOrEqual(num_points); i >= 2; i /= 2 {
		galEls = append(galEls, params.GaloisElement(-i))
	}
	eval = eval.WithKey(rlwe.NewMemEvaluationKeySet(rlk, kgen.GenGaloisKeysNew(galEls, sk)...))

	fmt.Println("== GENERATE BOOTSTRAP EVALUATOR ==")
	btpParametersLit := bootstrapping.ParametersLiteral{
		LogN: utils.Pointy(LogN),
	}
	btpParams, err := bootstrapping.NewParametersFromLiteral(params, btpParametersLit)
	if err != nil {
		panic(err)
	}
	fmt.Println("Generating bootstrapping evaluation keys...")
	evk_boot, _, err := btpParams.GenEvaluationKeys(sk)
	if err != nil {
		panic(err)
	}
	fmt.Println("Done")
	var eval_boot *bootstrapping.Evaluator
	if eval_boot, err = bootstrapping.NewEvaluator(btpParams, evk_boot); err != nil {
		panic(err)
	}
	fmt.Printf("Pointer address: %p\n", eval_boot)

	fmt.Println("== GENERATE INNRESUM EVALUATOR ==")
	batch := 1
	n := num_points
	eval_Innersum := eval.WithKey(rlwe.NewMemEvaluationKeySet(rlk, kgen.GenGaloisKeysNew(params.GaloisElementsForInnerSum(batch, n), sk)...))
	endTime := time.Now()
	fmt.Printf("PARAMETERS GENERATION TIME: %s\n", endTime.Sub(startTime))


	fmt.Printf("== ENCRYPT POINTS AND CENTERS ==\n")

	file_path := "../../data/dataset/mnist_train_reduced_64.csv"
	points_allbatches, maximumDistance := utilities.ReadMNIST_Reduced(file_path)
	points_allbatches = utilities.Rescale(points_allbatches, maximumDistance)
	// index_random := utilities.GenerateRandomNumbers_GivenInterval(0, 60000, 2)
	index_random := []int{52078, 16909}
	points_allbatches = utilities.TransposeMatrix(points_allbatches)
	points_allbatches = append(points_allbatches, points_allbatches[index_random[0]])
	points_allbatches = append(points_allbatches, points_allbatches[index_random[1]])
	points_allbatches = utilities.TransposeMatrix(points_allbatches)

	num_batch := 19
	points_allbatches_ct := make([][]*rlwe.Ciphertext, num_batch)
	var centers_ct []*rlwe.Ciphertext
	for i := 0; i < num_batch; i++ {
		points := utilities.SelectContinuousColumns(points_allbatches, i*num_points, (i+1)*num_points)

		if i == 0 {
			// index_centers := utilities.GenerateRandomNumbers_GivenInterval(0, num_points, num_centers)
			index_centers := make([]int, num_centers)
			for i := 0; i < num_centers; i++ {
				index_centers[i] = i
			}
			centers := plaintext.ExtractCentersFromPoints(points, index_centers)
			centers_ct = elaborate.EncryptCenters_NoNeedExtract_2(params, encoder, encryptor, centers, dimension, num_points, num_centers)
		}

		points_ct := minibatch.EncryptPoints_OneBatch(params, encoder, encryptor, points, dimension, num_points)
		points_allbatches_ct[i] = points_ct	
	}
	fmt.Printf("== DONE. ENCRYPT POINTS AND CENTERS ==\n")


	var bool_ct *rlwe.Ciphertext

	start_time2 := time.Now()
	iterations_updateCenters := 19
	for k := 0; k <= iterations_updateCenters; k++ {

		fmt.Printf("=============================== NOW AT THE %dth ITERATION ===============================\n", k)
		
		fmt.Printf("== EXTRACT NEW DATA POINTS ==\n")
		points_formated_ct := minibatch.ExtractPoints_OneBatch(eval, points_allbatches_ct[k % num_batch], num_points, dimension, num_centers)
		printers.PrintFirstOfFragment(points_formated_ct[0], encoder, decryptor, num_points, precision_print)

		fmt.Printf("== COMPUTE DISTANCE ==\n")
		distances_ct := elaborate.ComputeDistance_2(eval, points_formated_ct, centers_ct, dimension)
		printers.PrintFirstOfFragment(distances_ct, encoder, decryptor, num_points, precision_print)

		fmt.Printf("== DISTANCE COMPARISON ==\n")
		bool_ct = elaborate.CompareDistance_2(params, eval, distances_ct, dimension, num_points, num_centers, eval_boot)

		if k == iterations_updateCenters {
			break
		}
		fmt.Printf("== UPDATE CENTERS ==\n")
		centers_ct = elaborate.UpdateAllCenters_2(params, eval, eval_Innersum, points_formated_ct, bool_ct, centers_ct, dimension, num_points, num_centers, eval_boot)
		printers.PrintFirstOfFragment(centers_ct[0], encoder, decryptor, num_points, precision_print)
		printers.PrintFirstOfFragment(centers_ct[1], encoder, decryptor, num_points, precision_print)
		printers.PrintFirstOfFragment(centers_ct[63], encoder, decryptor, num_points, precision_print)
	}
	end_time2 := time.Now()
	fmt.Printf("PROGRAM EXECUTION TIME: %s\n", end_time2.Sub(start_time2))

	fmt.Printf("=============================== PRINT CHECK RESULT ===============================\n")



}
