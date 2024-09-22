//============================================================================
// Name        : Logic_Heaan.cpp
// Author      : Dongwoo Kim
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================


#include <iostream>
#include <iomanip>
#include <HEAAN.h>
#include <complex>
#include <math.h>
#include <algorithm>
#include <NTL/BasicThreadPool.h>

using namespace std;

int precision = 16;


// order: from c_0 to c_n
double fcoeff[4][5] = { { 3.0 / 2, -1.0 / 2, 0, 0, 0 },
{ 15.0 / 8, -10.0 / 8, 3.0 / 8, 0, 0 },
{ 35.0 / 16, -35.0 / 16, 21.0 / 16, -5.0 / 16, 0 },
{ 315.0 / 128, -420.0 / 128, 378.0 / 128, -180.0 / 128, 35.0 / 128 } };

// order: from c_0 to c_n
double encfcoeff[5] = {315.0, -420.0, 378.0, -180.0, 35.0};					// div by 128 = 2^7
unsigned long intCoeffNoSign[5] = {315, 420, 378, 180, 35};								// div by 128 = 2^7
double encfcoeff3[4] = { 35.0, -35.0, 21.0, -5.0};							// div by 16 = 2^4
double encfcoeff1[2] = { 3.0, -1.0};										// div by 2 = 2^1

//double encgcoeff[5] = {5.7129, -34.1540, 94.7418, -110.8312, 45.5305};
double encgcoeff[5] = {5.713801350009095, -34.145770080711120, 94.690589987263750, -110.7468429981996, 45.488221741637894 };
double encgcoeff3[4] = {4.4820, -16.1884, 25.0141, -12.5577};
double encgcoeff1[2] = {2.0766, -1.3266};


// order: from c_0 to c_n
double gcoeff[4][5] = { { 2.0766, -1.3266, 0, 0, 0 },
{ 3.2565, -5.9644, 3.7079, 0, 0 },
{ 4.4820, -16.1884, 25.0141, -12.5577, 0 },
{ 5850.0/1024, -34974.0/1024, 97015.0/1024, -113492.0/1024, 46623.0/1024 } };

// order: from c_0 to c_n
unsigned long intCoeffNoSignG[5] = {5850, 34974, 97015, 113492, 46623};								// div by 2^10


// compute f_n(x)
double evalf(double x, int n) {
	double xsqr = x*x;
	double xtmp = x;
	double res = 0;
	for (int i = 0; i < 5; ++i) {
		res += xtmp * fcoeff[n - 1][i];
		xtmp *= xsqr;
	}

	return res;
}

// compute g_n(x)
double evalg(double x, int n) {
	double xsqr = x*x;
	double xtmp = x;
	double res = 0;
	for (int i = 0; i < 5; ++i) {
		res += xtmp * gcoeff[n - 1][i];
		xtmp *= xsqr;
	}

	return res;
}


/// Encrypted Algorithm

void pow2(Ciphertext& res, Ciphertext input, long logn, long logp, Scheme scheme){
	res.copy(input);
	for (long i = 0; i < logn; ++i){
		scheme.multAndEqual(res, res);
		scheme.reScaleByAndEqual(res, logp);
	}
}


void evalInv(Ciphertext& res, Ciphertext input, long iter, long logp, Scheme scheme){
	Ciphertext _input;
	scheme.negate(_input, input);
	Ciphertext a;
	Ciphertext b;

	scheme.addConst(a, _input, 2.0, logp);		// a_0
	scheme.modDownByAndEqual(a, logp);
	scheme.addConst(b, _input, 1.0, logp);		// b_0

	Ciphertext tmp;
	for (long i = 0; i < iter; ++i){
		scheme.multAndEqual(b, b);				// b' = b^2
		scheme.reScaleByAndEqual(b, logp);
		scheme.addConst(tmp, b, 1.0, logp);		// b' + 1
		scheme.multAndEqual(a, tmp);			// a' = a * (b' + 1)
		scheme.reScaleByAndEqual(a, logp);
		cout << "inv a.logq: " << a.logq << endl;
	}
	
	res.copy(a);
	cout << "inv end" <<endl;
}

	
void evalSqr(Ciphertext& res, Ciphertext input, long iter, long logp, Scheme scheme){
	Ciphertext a = Ciphertext(input);
	Ciphertext b = Ciphertext(input);
	scheme.addConstAndEqual(b, -1.0, logp);

	Ciphertext tmp;
	for (long i = 0; i < iter; ++i){
		scheme.divByPo2(tmp, b, 1);
		scheme.negateAndEqual(tmp);
		scheme.addConstAndEqual(tmp, 1.0, logp);	// 1 - b/2
		scheme.modDownToAndEqual(a, tmp.logq);
		scheme.multAndEqual(a, tmp);
		scheme.reScaleByAndEqual(a, logp);			// a' = a ( 1 - b/2)

		scheme.addConst(tmp, b, -3.0, logp);
		scheme.divByPo2AndEqual(tmp, 2);			// (b-3)/4

		scheme.squareAndEqual(b);
		scheme.reScaleByAndEqual(b, logp);			// b' = b^2
		scheme.modDownToAndEqual(tmp, b.logq);
		scheme.multAndEqual(b, tmp);
		scheme.reScaleByAndEqual(b, logp);			// b' = b^2 * (b-3)/4

		cout << "a logq: " << a.logq << endl;
	}
	res.copy(a);
}


void evalMinMax(Ciphertext& min, Ciphertext& max, Ciphertext in_a, Ciphertext in_b, long iter, long logp, Scheme scheme){
	Ciphertext x;
	scheme.add(x, in_a, in_b);
	scheme.divByPo2AndEqual(x, 1);

	Ciphertext y;
	scheme.sub(y, in_a, in_b);
	cout << "y logq: " << y.logq <<endl;
	scheme.divByPo2AndEqual(y, 1);
	scheme.squareAndEqual(y);
	scheme.reScaleByAndEqual(y, logp);

	Ciphertext z;
	evalSqr(z, y, iter, logp, scheme);
	cout << "z logq: " << z.logq <<endl;

	scheme.modDownToAndEqual(x, z.logq);
	scheme.sub(min, x, z);
	scheme.add(max, x, z);
}


void evalComp(Ciphertext& resa, Ciphertext& resb, Ciphertext in_a, Ciphertext in_b, long iter, long* iter_inv, long logm, long logp, Scheme scheme){
	Ciphertext a = Ciphertext(in_a);
	Ciphertext b = Ciphertext(in_b);

	Ciphertext avg;
	scheme.add(avg, a, b);
	scheme.divByPo2AndEqual(avg, 1);			// avg = (a+b)/2
	
	Ciphertext inv;
	evalInv(inv, avg, iter_inv[0], logp, scheme);

	scheme.divByPo2AndEqual(a, 1);
	scheme.modDownToAndEqual(a, inv.logq);
	scheme.multAndEqual(a, inv);
	scheme.reScaleByAndEqual(a, logp);			// a_0 = a/2 * inv (avg, iter_inv)

	scheme.negate(b, a);
	scheme.addConstAndEqual(b, 1.0, logp); 		// b_0 = 1 - a_0

	Ciphertext am;
	Ciphertext bm;
	for (long i = 0; i < iter; ++i){
		pow2(am, a, logm, logp, scheme);
		pow2(bm, b, logm, logp, scheme);
		Ciphertext tmp;
		scheme.add(tmp, am, bm);
		evalInv(inv, tmp, iter_inv[1], logp, scheme);		// approx. of  1 / (a^m + b^m)

		if (i%2 != 0){
			scheme.modDownToAndEqual(am, inv.logq);
			scheme.mult(a, am, inv);
			scheme.reScaleByAndEqual(a, logp);					// a' = a^m * approx. of ''

			scheme.negate(b, a);
			scheme.addConstAndEqual(b, 1.0, logp);				// b' = 1 - a'
		}
		else {
			scheme.modDownToAndEqual(bm, inv.logq);
			scheme.mult(b, bm, inv);
			scheme.reScaleByAndEqual(b, logp);					// b' = b^m * approx. of ''

			scheme.negate(a, b);
			scheme.addConstAndEqual(a, 1.0, logp);				// a' = 1 - b'
		}
		cout << "a, b logq: " << a.logq << endl;
		cout << i << "iter" <<endl;
	}

	resa.copy(a);
	resb.copy(b);

	cout << "comp end" <<endl;
}


// For fdeg = 1
void evalF1(Ciphertext& res, Ciphertext in, long logp, Scheme scheme) {
	Ciphertext sqrx;
	Ciphertext x = Ciphertext(in);
	scheme.square(sqrx, in);								// t = x^2
	scheme.reScaleByAndEqual(sqrx, logp);

	// Compute a1t + a0

	scheme.multByConst(res, sqrx, encfcoeff1[1], logp);		// a1t
	scheme.reScaleByAndEqual(res, logp);
	scheme.addConstAndEqual(res, encfcoeff1[0], logp);		// a1t + a0

	scheme.modDownToAndEqual(x, res.logq);
	scheme.multAndEqual(res, x);							// res = (a1t + a0)x
	scheme.reScaleByAndEqual(res, logp);
	scheme.divByPo2AndEqual(res, 1);
}

// For fdeg = 3
void evalF3(Ciphertext& res, Ciphertext in, long logp, Scheme scheme) {
	Ciphertext sqrx, sqrx2;
	Ciphertext x = Ciphertext(in);
	scheme.square(sqrx, in);								// t = x^2
	scheme.reScaleByAndEqual(sqrx, logp);
	scheme.square(sqrx2, sqrx);								// t^2
	scheme.reScaleByAndEqual(sqrx2, logp);
	scheme.modDownToAndEqual(sqrx, sqrx2.logq);

	// Compute (a3t)t^2 + a2t^2 + a1t + a0

	Ciphertext tmp1, tmp2, tmp3;

	scheme.multByConst(tmp3, sqrx, encfcoeff3[3], logp);		// a3t
	scheme.multByConst(tmp2, sqrx2, encfcoeff3[2], logp);		// a2t^2
	scheme.multByConst(tmp1, sqrx, encfcoeff3[1], logp);		// a1t

	scheme.reScaleByAndEqual(tmp3, logp);
	scheme.modDownToAndEqual(sqrx2, tmp3.logq);

	scheme.multAndEqual(tmp3, sqrx2);						// (a3t)t^2
	scheme.reScaleByAndEqual(tmp3, logp);
	scheme.addAndEqual(tmp2, tmp1);							// a2t^2 + a1t
	scheme.addConstAndEqual(tmp2, encfcoeff3[0], 2 * logp);	// a2t^2 + a1t + a0
	scheme.reScaleByAndEqual(tmp2, logp);
	
	scheme.modDownToAndEqual(tmp2, tmp3.logq);
	scheme.add(res, tmp2, tmp3);							// res = (a3t)t^2 + a2t^2 + a1t + a0
	scheme.modDownToAndEqual(x, res.logq);
	scheme.multAndEqual(res, x);							// res = ((a4t^2 + a3t)t^2 + a2t^2 + a1t + a0)x
	scheme.reScaleByAndEqual(res, logp);
	scheme.divByPo2AndEqual(res, 4);
}

//// For fdeg = 4 , simple addition
//void evalF(Ciphertext& res, Ciphertext in, long logp, Scheme scheme, SecretKey secretKey) {
//	Ciphertext sqrx, sqrx2;
//	Ciphertext x = Ciphertext(in);
//	scheme.square(sqrx, x);								// t = x^2
//	scheme.reScaleByAndEqual(sqrx, logp);
//	scheme.square(sqrx2, sqrx);								// t^2
//	scheme.reScaleByAndEqual(sqrx2, logp);
//	scheme.modDownToAndEqual(sqrx, sqrx2.logq);
//
//	// Compute (a4t^2 + a3t)t^2 + a2t^2 + a1t + a0
//
//
//	Ciphertext tmp[5];
//
//	tmp[2].copy(sqrx2);
//	tmp[4].copy(sqrx2);
//	tmp[1].copy(sqrx);
//	tmp[3].copy(sqrx);
//
//
//	for (int j = 1; j < 5; ++j) {
//		tmp[0].copy(tmp[j]);
//		for (int i = 0; i < intCoeffNoSign[j]; ++i)
//			scheme.addAndEqual(tmp[j], tmp[0]);
//	}
//
//	complex<double>* check = scheme.decrypt(secretKey, tmp[1]);
//	cout << "tmp1: " << fixed << setprecision(precision) << real(check[0]) << endl;
//	check = scheme.decrypt(secretKey, tmp[2]);
//	cout << "tmp2: " << fixed << setprecision(precision) << real(check[0]) << endl;
//	check = scheme.decrypt(secretKey, tmp[3]);
//	cout << "tmp3: " << fixed << setprecision(precision) << real(check[0]) << endl;
//	check = scheme.decrypt(secretKey, tmp[4]);
//	cout << "tmp4: " << fixed << setprecision(precision) << real(check[0]) << endl;
//
//	scheme.subAndEqual(tmp[4], tmp[3]);										// a4t^2 + a3t
//	scheme.subAndEqual(tmp[2], tmp[1]);										// a2t^2 + a1t
//	scheme.addConstAndEqual(tmp[2], (double)intCoeffNoSign[0], logp);		// a2t^2 + a1t + a0
//
//	scheme.divByPo2AndEqual(tmp[4],7);
//	scheme.divByPo2AndEqual(tmp[2],7);
//
//	check = scheme.decrypt(secretKey, tmp[4]);
//	cout << "add final tmp4: " << fixed << setprecision(precision) << real(check[0]) << endl;
//
//	check = scheme.decrypt(secretKey, sqrx2);
//	cout << "sqrx2: " << fixed << setprecision(precision) << real(check[0]) << endl;
//
//	//cout <<"tmp4 : " <<tmp[4].logq <<", " <<tmp[4].logp <<" sqrx2: "<<sqrx2.logq <<", " <<sqrx2.logp<<endl;
//
//	scheme.modDownToAndEqual(sqrx2, tmp[4].logq);
//	scheme.multAndEqual(tmp[4], sqrx2);										// (a4t^2 + a3t)t^2
//	scheme.reScaleByAndEqual(tmp[4], logp);
//	scheme.modDownToAndEqual(tmp[2], tmp[4].logq);
//
//	check = scheme.decrypt(secretKey, tmp[2]);
//	cout << "final tmp2: " << fixed << setprecision(precision) << real(check[0]) << endl;
//	check = scheme.decrypt(secretKey, tmp[4]);
//	cout << "final tmp4: " << fixed << setprecision(precision) << real(check[0]) << endl;
//
//
//	scheme.add(res, tmp[2], tmp[4]);								// res = (a4t^2 + a3t)t^2 + a2t^2 + a1t + a0
//
//	//cout << "preres, logp: " << res.logp << " , " << res.logq << endl;
//	check = scheme.decrypt(secretKey, res);
//	cout << "preres: " << fixed << setprecision(precision) << real(check[0]) << endl;
//
//	scheme.modDownToAndEqual(x, res.logq);
//	scheme.multAndEqual(res, x);									// res = ((a4t^2 + a3t)t^2 + a2t^2 + a1t + a0)x
//	scheme.reScaleByAndEqual(res, logp);
//	//scheme.divByPo2AndEqual(res, 7);
//
//	cout << "res, logp: " << res.logp << " , " << res.logq << endl;
//	check = scheme.decrypt(secretKey, res);
//	cout << "res: " << fixed << setprecision(precision) << real(check[0]) << endl;
//}

// For fdeg = 4
void evalF(Ciphertext& res, Ciphertext in, long logp, Scheme scheme, SecretKey secretKey) {
	Ciphertext sqrx, sqrx2;
	Ciphertext x = Ciphertext(in);
	scheme.square(sqrx, x);								// t = x^2
	scheme.reScaleByAndEqual(sqrx, logp);
	scheme.square(sqrx2, sqrx);								// t^2
	scheme.reScaleByAndEqual(sqrx2, logp);
	scheme.modDownToAndEqual(sqrx, sqrx2.logq);

	// Compute (a4t^2 + a3t)t^2 + a2t^2 + a1t + a0

	Ciphertext tmp[5];
	Ciphertext ctmp[5];

	tmp[2].copy(sqrx2);
	tmp[4].copy(sqrx2);
	tmp[1].copy(sqrx);
	tmp[3].copy(sqrx);

	unsigned int index[5];
	bool init[5];
	for (int i = 1; i < 5; ++i) {
		index[i] = intCoeffNoSign[i];
		init[i] = true;
	}

	// coeffs (bounded by 9 bit)
	for (int i = 0; i < 9; ++i) {
		for (int j = 1; j < 5; ++j){

			// compute a_j t^j
			if (index[j] & 1) {
				if (init[j]) {
					ctmp[j].copy(tmp[j]);
					init[j] = false;
				}
				else {
					scheme.addAndEqual(ctmp[j], tmp[j]);
				}
			}
			scheme.addAndEqual(tmp[j],tmp[j]);
			index[j] >>= 1;
		}
	}

	//complex<double>* check;
	//check = scheme.decrypt(secretKey, ctmp[1]);
	//cout << "tmp1: " << fixed << setprecision(precision) << real(check[0]) << endl;
	//check = scheme.decrypt(secretKey, ctmp[2]);
	//cout << "tmp2: " << fixed << setprecision(precision) << real(check[0]) << endl;
	//check = scheme.decrypt(secretKey, ctmp[3]);
	//cout << "tmp3: " << fixed << setprecision(precision) << real(check[0]) << endl;
	//check = scheme.decrypt(secretKey, ctmp[4]);
	//cout << "tmp4: " << fixed << setprecision(precision) << real(check[0]) << endl;

	scheme.subAndEqual(ctmp[4], ctmp[3]);									// a4t^2 + a3t
	scheme.subAndEqual(ctmp[2], ctmp[1]);									// a2t^2 + a1t
	scheme.addConstAndEqual(ctmp[2], (double) intCoeffNoSign[0], logp);		// a2t^2 + a1t + a0
	scheme.divByPo2AndEqual(ctmp[4], 7);
	scheme.divByPo2AndEqual(ctmp[2], 7);

	//check = scheme.decrypt(secretKey, ctmp[4]);
	//cout << "second" << endl;
	//cout << "tmp4: " << fixed << setprecision(precision) << real(check[0]) << endl;
	//check = scheme.decrypt(secretKey, ctmp[2]);
	//cout << "tmp2: " << fixed << setprecision(precision) << real(check[0]) << endl;

	scheme.modDownToAndEqual(sqrx2, ctmp[4].logq);
	scheme.multAndEqual(ctmp[4], sqrx2);										// (a4t^2 + a3t)t^2
	scheme.reScaleByAndEqual(ctmp[4], logp);
	scheme.modDownToAndEqual(ctmp[2], ctmp[4].logq);

	scheme.add(res, ctmp[2], ctmp[4]);						// res = (a4t^2 + a3t)t^2 + a2t^2 + a1t + a0

	scheme.modDownToAndEqual(x, res.logq);
	scheme.multAndEqual(res, x);							// res = ((a4t^2 + a3t)t^2 + a2t^2 + a1t + a0)x
	scheme.reScaleByAndEqual(res, logp);

	//cout << "res, logp: " << res.logp << " , " << res.logq << endl;
	//check = scheme.decrypt(secretKey, res);
	//cout << "res: " << fixed << setprecision(precision) << real(check[0]) << endl;
}
 
//// original fdeg = 4
//void evalF(Ciphertext& res, Ciphertext in, long logp, Scheme scheme, SecretKey secretKey) {
//	Ciphertext sqrx, sqrx2;
//	Ciphertext x = Ciphertext(in);
//	scheme.square(sqrx, x);								// t = x^2
//	scheme.reScaleByAndEqual(sqrx, logp);
//	scheme.square(sqrx2, sqrx);								// t^2
//	scheme.reScaleByAndEqual(sqrx2, logp);
//	scheme.modDownToAndEqual(sqrx, sqrx2.logq);
//
//	 // Compute (a4t^2 + a3t)t^2 + a2t^2 + a1t + a0
//
//	Ciphertext tmp1, tmp2, tmp3, tmp4;
//
//	scheme.multByConst(tmp4, sqrx2, encfcoeff[4], logp);	// a4t^2
//	scheme.multByConst(tmp2, sqrx2, encfcoeff[2], logp);	// a2t^2
//
//	scheme.multByConst(tmp3, sqrx, encfcoeff[3], logp);		// a3t
//	scheme.multByConst(tmp1, sqrx, encfcoeff[1], logp);		// a1t
//
//	scheme.addAndEqual(tmp4, tmp3);							// a4t^2 + a3t
//	scheme.addAndEqual(tmp2, tmp1);							// a2t^2 + a1t
//	scheme.addConstAndEqual(tmp2, encfcoeff[0], 2 * logp);	// a2t^2 + a1t + a0
//	scheme.reScaleByAndEqual(tmp4, logp);
//	scheme.reScaleByAndEqual(tmp2, logp);
//
//	scheme.modDownToAndEqual(sqrx2, tmp4.logq);
//
//	scheme.multAndEqual(tmp4, sqrx2);						// (a4t^2 + a3t)t^2
//	scheme.reScaleByAndEqual(tmp4, logp);
//	scheme.modDownToAndEqual(tmp2, tmp4.logq);
//
//	scheme.add(res, tmp2, tmp4);							// res = (a4t^2 + a3t)t^2 + a2t^2 + a1t + a0
//	scheme.modDownToAndEqual(x, res.logq);
//	scheme.multAndEqual(res, x);							// res = ((a4t^2 + a3t)t^2 + a2t^2 + a1t + a0)x
//	scheme.reScaleByAndEqual(res, logp);
//	scheme.divByPo2AndEqual(res, 7);
//}


// For gdeg = 1
void evalG1(Ciphertext& res, Ciphertext in, long logp, Scheme scheme) {
	Ciphertext sqrx;
	Ciphertext x = Ciphertext(in);
	scheme.square(sqrx, in);								// t = x^2
	scheme.reScaleByAndEqual(sqrx, logp);

	// Compute a1t + a0

	scheme.multByConst(res, sqrx, encgcoeff1[1], logp);		// a1t
	scheme.reScaleByAndEqual(res, logp);
	scheme.addConstAndEqual(res, encgcoeff1[0], logp);		// a1t + a0

	scheme.modDownToAndEqual(x, res.logq);
	scheme.multAndEqual(res, x);							// res = (a1t + a0)x
	scheme.reScaleByAndEqual(res, logp);
}


// For gdeg = 3
void evalG3(Ciphertext& res, Ciphertext in, long logp, Scheme scheme) {
	Ciphertext sqrx, sqrx2;
	Ciphertext x = Ciphertext(in);
	scheme.square(sqrx, in);								// t = x^2
	scheme.reScaleByAndEqual(sqrx, logp);
	scheme.square(sqrx2, sqrx);								// t^2
	scheme.reScaleByAndEqual(sqrx2, logp);
	scheme.modDownToAndEqual(sqrx, sqrx2.logq);

	// Compute (a3t)t^2 + a2t^2 + a1t + a0

	Ciphertext tmp1, tmp2, tmp3;

	scheme.multByConst(tmp3, sqrx, encgcoeff3[3], logp);		// a3t
	scheme.multByConst(tmp2, sqrx2, encgcoeff3[2], logp);	// a2t^2
	scheme.multByConst(tmp1, sqrx, encgcoeff3[1], logp);		// a1t

	scheme.reScaleByAndEqual(tmp3, logp);
	scheme.modDownToAndEqual(sqrx2, tmp3.logq);

	scheme.multAndEqual(tmp3, sqrx2);						// (a3t)t^2
	scheme.reScaleByAndEqual(tmp3, logp);
	scheme.addAndEqual(tmp2, tmp1);							// a2t^2 + a1t
	scheme.addConstAndEqual(tmp2, encgcoeff3[0], 2 * logp);	// a2t^2 + a1t + a0
	scheme.reScaleByAndEqual(tmp2, logp);

	scheme.modDownToAndEqual(tmp2, tmp3.logq);
	scheme.add(res, tmp2, tmp3);							// res = (a3t)t^2 + a2t^2 + a1t + a0
	scheme.modDownToAndEqual(x, res.logq);
	scheme.multAndEqual(res, x);							// res = ((a4t^2 + a3t)t^2 + a2t^2 + a1t + a0)x
	scheme.reScaleByAndEqual(res, logp);
}

//// For gdeg = 4		// need to be changed // original G
//void evalG(Ciphertext& res, Ciphertext in, int gdeg, long logp, Scheme scheme) {
//	Ciphertext sqrx, sqrx2;
//	Ciphertext x = Ciphertext(in);
//	scheme.square(sqrx, in);								// t = x^2
//	scheme.reScaleByAndEqual(sqrx, logp);
//	scheme.square(sqrx2, sqrx);								// t^2
//	scheme.reScaleByAndEqual(sqrx2, logp);
//	scheme.modDownToAndEqual(sqrx, sqrx2.logq);
//
//	// Compute (a4t^2 + a3t)t^2 + a2t^2 + a1t + a0
//
//	Ciphertext tmp1, tmp2, tmp3, tmp4;
//
//	scheme.multByConst(tmp4, sqrx2, encgcoeff[4], logp);	// a4t^2
//	scheme.multByConst(tmp2, sqrx2, encgcoeff[2], logp);	// a2t^2
//
//	scheme.multByConst(tmp3, sqrx, encgcoeff[3], logp);		// a3t
//	scheme.multByConst(tmp1, sqrx, encgcoeff[1], logp);		// a1t
//
//	scheme.addAndEqual(tmp4, tmp3);							// a4t^2 + a3t
//	scheme.addAndEqual(tmp2, tmp1);							// a2t^2 + a1t
//	scheme.addConstAndEqual(tmp2, encgcoeff[0], 2 * logp);	// a2t^2 + a1t + a0
//	scheme.reScaleByAndEqual(tmp4, logp);
//	scheme.reScaleByAndEqual(tmp2, logp);
//	scheme.modDownToAndEqual(sqrx2, tmp4.logq);
//
//	scheme.multAndEqual(tmp4, sqrx2);						// (a4t^2 + a3t)t^2
//	scheme.reScaleByAndEqual(tmp4, logp);
//	scheme.modDownToAndEqual(tmp2, tmp4.logq);
//
//	scheme.add(res, tmp2, tmp4);							// res = (a4t^2 + a3t)t^2 + a2t^2 + a1t + a0
//	scheme.modDownToAndEqual(x, res.logq);
//	scheme.multAndEqual(res, x);							// res = ((a4t^2 + a3t)t^2 + a2t^2 + a1t + a0)x
//	scheme.reScaleByAndEqual(res, logp);
//}

// For gdeg = 4
void evalG(Ciphertext& res, Ciphertext in, long logp, Scheme scheme, SecretKey secretKey) {
	Ciphertext sqrx, sqrx2;
	Ciphertext x = Ciphertext(in);
	scheme.square(sqrx, x);								// t = x^2
	scheme.reScaleByAndEqual(sqrx, logp);
	scheme.square(sqrx2, sqrx);								// t^2
	scheme.reScaleByAndEqual(sqrx2, logp);
	scheme.modDownToAndEqual(sqrx, sqrx2.logq);

	// Compute (a4t^2 + a3t)t^2 + a2t^2 + a1t + a0

	Ciphertext tmp[5];
	Ciphertext ctmp[5];

	tmp[2].copy(sqrx2);
	tmp[4].copy(sqrx2);
	tmp[1].copy(sqrx);
	tmp[3].copy(sqrx);
	
	unsigned long index[5];
	bool init[5];
	for (int i = 1; i < 5; ++i) {
		index[i] = intCoeffNoSignG[i];
		init[i] = true;
	}

	// coeffs (bounded by 17 bit)
	for (int i = 0; i < 17; ++i) {
		for (int j = 1; j < 5; ++j) {

			// compute a_j t^j
			if (index[j] & 1) {
				if (init[j]) {
					ctmp[j].copy(tmp[j]);
					init[j] = false;
				}
				else {
					scheme.addAndEqual(ctmp[j], tmp[j]);
				}
			}
			scheme.addAndEqual(tmp[j], tmp[j]);
			index[j] >>= 1;
		}
	}

	//complex<double>* check;
	//check = scheme.decrypt(secretKey, ctmp[1]);
	//cout << "tmp1: " << fixed << setprecision(precision) << real(check[0]) << endl;
	//check = scheme.decrypt(secretKey, ctmp[2]);
	//cout << "tmp2: " << fixed << setprecision(precision) << real(check[0]) << endl;
	//check = scheme.decrypt(secretKey, ctmp[3]);
	//cout << "tmp3: " << fixed << setprecision(precision) << real(check[0]) << endl;
	//check = scheme.decrypt(secretKey, ctmp[4]);
	//cout << "tmp4: " << fixed << setprecision(precision) << real(check[0]) << endl;


	scheme.subAndEqual(ctmp[4], ctmp[3]);									// a4t^2 + a3t
	scheme.subAndEqual(ctmp[2], ctmp[1]);									// a2t^2 + a1t
	scheme.addConstAndEqual(ctmp[2], (double)intCoeffNoSignG[0], logp);		// a2t^2 + a1t + a0
	scheme.divByPo2AndEqual(ctmp[4], 10);
	scheme.divByPo2AndEqual(ctmp[2], 10);

	//check = scheme.decrypt(secretKey, ctmp[4]);
	//cout << "second" << endl;
	//cout << "tmp4: " << fixed << setprecision(precision) << real(check[0]) << endl;
	//check = scheme.decrypt(secretKey, ctmp[2]);
	//cout << "tmp2: " << fixed << setprecision(precision) << real(check[0]) << endl;


	scheme.modDownToAndEqual(sqrx2, ctmp[4].logq);
	scheme.multAndEqual(ctmp[4], sqrx2);										// (a4t^2 + a3t)t^2
	scheme.reScaleByAndEqual(ctmp[4], logp);
	scheme.modDownToAndEqual(ctmp[2], ctmp[4].logq);

	scheme.add(res, ctmp[2], ctmp[4]);						// res = (a4t^2 + a3t)t^2 + a2t^2 + a1t + a0

	scheme.modDownToAndEqual(x, res.logq);
	scheme.multAndEqual(res, x);							// res = ((a4t^2 + a3t)t^2 + a2t^2 + a1t + a0)x
	scheme.reScaleByAndEqual(res, logp);
}




void evalNewComp(Ciphertext& resa, Ciphertext in_a, Ciphertext in_b, long iter, long logp, Scheme scheme, SecretKey secretKey) {

	Ciphertext x = Ciphertext(in_a);
	scheme.subAndEqual(x, in_b);

	cout << "initial x logq: " << x.logq << endl;
	for (int i = 0; i < iter; ++i) {
		evalF(x, x, logp, scheme, secretKey);
		cout << "after " << i+1 << "-th iter" << endl;
		cout << "x logq: " << x.logq << endl;		
	}
	
	scheme.addConst(resa, x, 1.0, logp);
	scheme.divByPo2AndEqual(resa, 1);
}


void evalNewCompG(Ciphertext& resa, Ciphertext in_a, Ciphertext in_b, long iter_f, long iter_g, long logp, Scheme scheme, SecretKey secretKey) {
	
	Ciphertext x = Ciphertext(in_a);
	scheme.subAndEqual(x, in_b);

	cout << "initial x logq: " << x.logq << endl;
	for (int i = 0; i < iter_g; ++i) {
		evalG(x, x, logp, scheme, secretKey);
		cout << "after " << i + 1 << "-th iter" << endl;
		cout << "x logq: " << x.logq << endl;
	}
	cout << "end g" << endl;
	cout << "initial x logq: " << x.logq << endl;
	for (int i = 0; i < iter_f; ++i) {
		evalF(x, x, logp, scheme, secretKey);
		cout << "after " << i + 1 << "-th iter" << endl;
		cout << "x logq: " << x.logq << endl;
	}

	scheme.addConst(resa, x, 1.0, logp);
	scheme.divByPo2AndEqual(resa, 1);
}


void evalNewMax(Ciphertext& resa, Ciphertext in_a, Ciphertext in_b, long iter, long logp, Scheme scheme, SecretKey secretKey) {

	Ciphertext x = Ciphertext(in_a);
	scheme.subAndEqual(x, in_b);
	Ciphertext tmp = Ciphertext(x);
	scheme.add(resa, in_a, in_b);
	scheme.divByPo2AndEqual(resa, 1);
	
	cout << "initial x logq: " << x.logq << endl;
	for (int i = 0; i < iter; ++i) {
		evalF(x, x, logp, scheme, secretKey);
		cout << "after " << i + 1 << "-th iter" << endl;
		cout << "x logq: " << x.logq << endl;
	}

	scheme.modDownToAndEqual(tmp, x.logq);
	scheme.multAndEqual(x, tmp);
	scheme.reScaleByAndEqual(x, logp);
	
	scheme.divByPo2AndEqual(x, 1);
	scheme.modDownToAndEqual(resa, x.logq);
	scheme.addAndEqual(resa, x);	
}



//// plain algorithm


// compute x^{2^d}
double repeatsq(double x, int d) {
	double sqr = x;
	for (int i = 0; i < d; i++) {
		sqr = sqr * sqr;
	}
	return sqr;
}

// compute inv = inv(x; d) with d iteration
void inverse(double& inv, double x, long d) {
	double a = 2 - x;
	double b = 1 - x;
	for (long i = 0; i < d; i++) {
		b = b * b;
		a = a * (1 + b);
	}
	inv = a;
}

// compute root = sqr(x; d) with d iteration, c: bound of relative error
void sqroot(double& root, double& c, double x, int d) {
	root = x;
	c = x - 1;
	for (int i = 0; i < d; i++) {
		root = root * (1 - c / 2);
		c = c * c * ((c - 3) / 4);
	}
}

// compute maxv = max(x, y) with d iteration, alpha: bound of error (precision)
void sqrmax(double& maxv, double x, double y, int d, long alpha) {
	double tmp1 = (x - y) * (x - y), tmp2, error;
	sqroot(tmp2, error, tmp1, d);
	maxv = (x + y) / 2 + tmp2 / 2;

	if (abs(maxv - max(x, y)) >= pow(2, -alpha))
		cout << "error";
}


// compute res = comp(in_a, in_b) with iter iteration, iter_inv for inverse, logm
void comp(double* res, double in_a, double in_b, long iter, long* iter_inv, long logm) {
	double inv;
	double a;
	double b;
	double avg = (in_a + in_b) / 2;

	inverse(inv, avg, iter_inv[0]);

	a = in_a * inv * 0.5;
	b = in_b * inv * 0.5;

	for (long i = 0; i < iter; ++i) {

		double am = repeatsq(a, logm);
		double bm = repeatsq(b, logm);
		inverse(inv, am + bm, iter_inv[1]);

		if (i % 2 != 0) {
			a = am * inv;
			b = 1 - a;
		}
		else {
			b = bm * inv;
			a = 1 - b;
		}
	}
	res[0] = a;
	res[1] = b;
}


void printVecSqrmax(complex<double>* input, long slots, long alpha, long worstC) {
	for (long i = 0; i < slots; ++i) {
		cout << fixed << setprecision(precision) << real(input[i]) << ", ";

		if ((abs(real(input[i]) - 1) >= pow(2, -alpha)) && (i == 0))
			cout << "error";
		if ((abs(real(input[i]) - worstC*pow(2, -alpha)) >= pow(2, -alpha)) && (i == 1))
			cout << "erorr";
		if ((abs(real(input[i]) - 1) >= pow(2, -alpha)) && (i == 2))
			cout << "erorr";
		if ((abs(real(input[i]) - pow(2, -alpha)) >= pow(2, -alpha)) && (i == 3))
			cout << "erorr";
	}
	cout << "." << endl;
}


void printVecComp(complex<double>* input, long slots, long alpha) {
	for (long i = 0; i < slots; ++i) {
		cout << fixed << setprecision(precision) << real(input[i]) << ", ";

		if ((abs((real(input[i]))) > pow(2, -alpha)) && (i == 0))
			cout << "error ";
		if ((abs(real(input[i]) - 1) > pow(2, -alpha)) && (i == 1))
			cout << "error ";
	}
	cout << "." << endl;
}


// compute res = NewMax(in_a, in_b) with iter and f_n
void NewMax(double& res, double in_a, double in_b, long iter, int n) {
	double x = in_a - in_b;
	res = (in_a + in_b) / 2;
	for (int i = 0; i < iter; ++i)
		x = evalf(x, n);

	res += (in_a - in_b)*x / 2;
}


// compute res = NewComp(in_a, in_b) with iter and f_n
void NewComp(double& res, double in_a, double in_b, long iter, int n) {
	double x = in_a - in_b;
	for (int i = 0; i < iter; ++i)
		x = evalf(x, n);

	res = (x + 1.0) / 2;
}

// compute res = NewComp(in_a, in_b) with iter and f_n and g_n
void NewCompG(double& res, double in_a, double in_b, long iter_f, long iter_g, int n) {
	double x = in_a - in_b;
	for (int i = 0; i < iter_g; ++i)
		x = evalg(x, n);
	for (int i = 0; i < iter_f; ++i)
		x = evalf(x, n);

	res = (x + 1.0) / 2;
}



int main() {

	// 4.1 NewMax (Renewed!)
	/*
	/// Parameters for HEAAN
	//// Modify corresponding parameters (logN, logQ) in Params.h file also !!, Then rebuild HEAAN!

	long logq = 1574; 			// the highest modulus of encryption (only for estimation of logQ)
	long logp = 40;				// precisionBits

	/// Params for Algorithm
	long iter = 4;				// Needs to be the same as iter_out estimated by program!
	long alpha = 8;
	long deg = 4;					// degree of f

	double prec = pow(2, -alpha);
	long num_inputs = (1 << alpha);				// total number of inputs = 2^alpha
	long pack_num = (1LL << 16);				// max number of slots in one ciphertext
	long num_slots;								// number of slots in test ciphertext
	long ctx_num = num_inputs / pack_num;		// number of ciphertexts
	if (ctx_num == 0) {
		num_slots = num_inputs;
		ctx_num = 1;
	}
	else
		num_slots = pack_num;
	

	SetNumThreads(8);
	srand(time(NULL));


	double* Avec = new double[num_inputs];		// values : 1*2^(-alpha), 2*2^(-alpha), 3*2^(-alpha), ... , (2^alpha - 2)*2^(-alpha), (2^alpha - 1)*2^(-alpha), 1
	double* Bvec = new double[num_inputs];		// values : all 0

	double* Rvec = new double[num_inputs];		// stores real max results (A > B)
	double* Pvec = new double[num_inputs];		// stores plain max results NewComp(A,B)
	complex<double>* Evec;
	// double* Evec = new double[num_inputs];		// stores enc max results NewComp(A,B) over HEAAN

	double* Perr = new double[num_inputs];		// stores errors (log2) of plain comparison results NewMax(A,B)
	double* Eerr = new double[num_inputs]; 		// stores errors (log2) of enc comparison results NewMax(A,B) over HEAAN

	double MaxErr = -100;						// stores max error
	long argMax = 0;							// stores the position of max error

	double worstC = pow(2, alpha + 1) / sqrt(1 + pow(2, iter + 2));						// worst value for NewMax
	
	long iter_exp = 0;
	long err_cnt;
	

	for (long i = 0; i < num_inputs; ++i) {
		Avec[i] = prec*(i + 1);
		Bvec[i] = 0;
		Rvec[i] = max(Avec[i], Bvec[i]);
	}

	//cout << "New Max (with Plain Algorithm) Expected iterations: ";

	//double maxout = 0;
	//for (; abs(maxout - worstC) > prec; ++iter_exp)
	//	NewMax(maxout, worstC, 0, iter_exp, deg);
	//cout << iter_exp << endl;
	
	cout << "New Max (with Plain Algorithm) Error Test: ";
	err_cnt = 0;
	for (long i = 0; i < num_inputs; ++i) {
		//sqrmax(Pvec[i], Avec[i], Bvec[i], iter, alpha);			// oldold max
		NewMax(Pvec[i], Avec[i], Bvec[i], iter, deg);				// new max (deg = 1 or 4) 
		Perr[i] = abs(Rvec[i] - Pvec[i]);
		if (Perr[i] > prec) {
			//cout << "error!" << endl;
			err_cnt++;
		}
		Perr[i] = log2(Perr[i]);
		if (Perr[i] > MaxErr) {
			MaxErr = Perr[i];
			argMax = i;
		}
	}
	if (err_cnt != 0) {
		cout << "Error occured: " << err_cnt++; cout << " at " << argMax << endl;
		cout << endl;
	}
	else
	cout << "Done." << endl;



	//// Start Enc Algorithm
	// Scheme & RotKey Gen
	TimeUtils timeutils;
	timeutils.start("Scheme Generation");
	Ring ring;
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);
	timeutils.stop("Scheme Generation");
	// scheme.addLeftRotKeys(secretKey);

	cout << "logN: " << logN << endl;
	cout << "logQ: " << logQ << endl;
	cout << "logp: " << logp << endl;

	cout << "-----------------------" << endl;

	long tot_err = 0;

	for (long test_iter = 0; test_iter < ctx_num; ++test_iter) {

		timeutils.start("Encrypt");
		Ciphertext ctxa, ctxb;
		scheme.encrypt(ctxa, &Avec[test_iter*num_slots], num_slots, logp, logq);
		scheme.encrypt(ctxb, &Bvec[test_iter*num_slots], num_slots, logp, logq);
		timeutils.stop("Encrypt");

		Ciphertext resa;

		timeutils.start("Evaluation");
		evalNewMax(resa, ctxa, ctxb, iter, logp, scheme, secretKey);
		timeutils.stop("Evaluation");

		Evec = scheme.decrypt(secretKey, resa);

		MaxErr = -100;
		cout << "New Comp G (with Enc Algorithm) Error Test: ";
		err_cnt = 0;
		for (long i = 0; i < num_slots; ++i) {
			//cout << Pvec[i] <<", "<< real(Evec[i]) << endl;
			Eerr[i] = abs(Rvec[test_iter*num_slots + i] - real(Evec[i]));
			if (Eerr[i] > prec) {
				//cout << "error!" << endl;
				err_cnt++;
				tot_err++;
			}
			Eerr[i] = log2(Eerr[i]);
			if (Eerr[i] > MaxErr) {
				MaxErr = Perr[i];
				argMax = i;
			}
		}
		if (err_cnt != 0) {
			cout << "Error occured: " << err_cnt++; cout << " at " << argMax << endl;
			cout << endl;
		}
		else
			cout << test_iter + 1 << "-th test / " << ctx_num << " Done." << endl;
	}

	cout << "total error: " << tot_err << endl;
	*/



	// 5.1 NewComp
	/// Parameters for HEAAN
	//// Modify corresponding parameters (logN, logQ) in Params.h file also !!, Then rebuild HEAAN!

	long logq = 2084; 			// the highest modulus of encryption (only for estimation of logQ)
	long logp = 40;				// precisionBits

/// Params for Algorithm
	long iter = 15;				// Needs to be the same as iter_out estimated by program!
	long iter_g = 10;
	long iter_f = 2;
	long alpha = 24;
	long deg = 4;					// degree of f, g

	double prec = pow(2, -alpha);
	long num_inputs = (1 << alpha);				// total number of inputs = 2^alpha
	long pack_num = (1LL << 16);				// max number of slots in one ciphertext
	long num_slots;								// number of slots in test ciphertext
	long ctx_num = num_inputs / pack_num;		// number of ciphertexts
	if (ctx_num == 0) {
		num_slots = num_inputs;
		ctx_num = 1;
	}
	else
		num_slots = pack_num;
		

	
	SetNumThreads(8);
	srand(time(NULL));


	double* Avec = new double[num_inputs];		// values : 1*2^(-alpha), 2*2^(-alpha), 3*2^(-alpha), ... , (2^alpha - 2)*2^(-alpha), (2^alpha - 1)*2^(-alpha), 1
	double* Bvec = new double[num_inputs];		// values : all 0

	double* Rvec = new double[num_inputs];		// stores real comparison results (A > B)
	double* Pvec = new double[num_inputs];		// stores plain comparison results NewComp(A,B)
	complex<double>* Evec;
	// double* Evec = new double[num_inputs];		// stores enc comparison results NewComp(A,B) over HEAAN

	double* Perr = new double[num_inputs];		// stores errors (log2) of plain comparison results NewComp(A,B)
	double* Eerr = new double[num_inputs]; 		// stores errors (log2) of enc comparison results NewComp(A,B) over HEAAN

	double MaxErr = -100;						// stores max error
	long argMax = 0;							// stores the position of max error

	double worstF = prec;						// worst value for NewComp
	double worstG[2];							// worst value for NewCompG, C = 0.75
	worstG[0] = prec;
	worstG[1] = 0.75;

	long iter_exp = 0;
	long err_cnt;


	for (long i = 0; i < num_inputs; ++i) {
		Avec[i] = prec*(i + 1);
		Bvec[i] = 0;
		Rvec[i] = (Avec[i] > Bvec[i]);
	}

	cout << "New Comp (with Plain Algorithm) Expected iterations: ";

	double test_iter = worstF;
	for (; abs(1 - test_iter) > prec; ++iter_exp)
		test_iter = evalf(test_iter, deg);
	cout << iter_exp << endl;

	cout << "New Comp G (with Plain Algorithm) Expected iterations: ";
	iter_exp = 0;
	test_iter = worstG[0];
	for (; test_iter < 0.75; ++iter_exp)
		test_iter = evalg(test_iter, deg);
	cout << "g: " << iter_exp;

	iter_exp = 0;
	test_iter = worstG[1];
	for (; abs(1 - test_iter) > prec; ++iter_exp)
		test_iter = evalf(test_iter, deg);
	cout << ", f: " << iter_exp << endl;


	cout << "New Comp (with Plain Algorithm) Error Test: ";
	err_cnt = 0;
	for (long i = 0; i < num_inputs; ++i) {
		NewComp(Pvec[i], Avec[i], Bvec[i], iter, deg);
		Perr[i] = abs(Rvec[i] - Pvec[i]);
		if (Perr[i] > prec) {
			//cout << "error!" << endl;
			err_cnt++;
		}
		Perr[i] = log2(Perr[i]);
		if (Perr[i] > MaxErr) {
			MaxErr = Perr[i];
			argMax = i;
		}
	}
	if (err_cnt != 0) {
		cout << "Error occured: " << err_cnt++; cout << " at " << argMax << endl;
		cout << endl;
	}
	else
		cout << "Done." << endl;


	cout << "New Comp G (with Plain Algorithm) Error Test: ";
	err_cnt = 0;
	for (long i = 0; i < num_inputs; ++i) {
		NewCompG(Pvec[i], Avec[i], Bvec[i], iter_f, iter_g, deg);
		Perr[i] = abs(Rvec[i] - Pvec[i]);
		if (Perr[i] > prec) {
			//cout << "error!" << endl;
			err_cnt++;
		}
		Perr[i] = log2(Perr[i]);
		if (Perr[i] > MaxErr) {
			MaxErr = Perr[i];
			argMax = i;
		}
	}
	if (err_cnt != 0) {
		cout << "Error occured: " << err_cnt++; cout << " at " << argMax << endl;
		cout << endl;
	}
	else
		cout << "Done." << endl;


	//// Uncomment two of the following lines for comparison of choice (exact, c = 1.01, c = 1.05 ; for description of c, see Paper p.15 or p.27 )
	//
	//	//double avec[2] = {1, 1};			// for exact comparison
	//	double avec[2] = {pow(2, -alpha), 1};			// for exact comparison
	//	double bvec[2] = {0, 0};
	////	double avec[2] = {1.5 - 0.01485, 1.5};			// for c = 1.01
	////	double bvec[2] = {1.5, 1.5 - 0.01485};
	////	double avec[2] = {1.5 - 0.07143, 1.5};			// for c = 1.05
	////	double bvec[2] = {1.5, 1.5 - 0.07143};
	//
	//	double tmpa = avec[0];
	//	long iter_out = 0;
	//	long iter_g_out = 0;
	//	long iter_f_out = 0;
	//
	//	for (; tmpa < 1 - pow(2, -alpha); ++iter_out) {
	//		tmpa = evalf(tmpa, deg);
	//	}
	//
	//	cout << "For NewComp " << endl;
	//	cout << "iter: " << iter_out << endl;
	//	cout << "tmpa: " << tmpa << endl;
	//
	//	tmpa = avec[0];
	//
	//	for (; tmpa < 0.75; ++iter_g_out) {
	//		tmpa = evalg(tmpa, deg);
	//	}
	//	double tmpaa = 0.75;
	//	for (; tmpaa < 1-pow(2, -alpha); ++iter_f_out) {
	//		tmpaa = evalf(tmpaa, deg);
	//	}
	//	cout << "For NewComp with g," << endl;
	//	cout << "tmpa: " << tmpa << endl;
	//	cout << "tmpaa: " << tmpaa << endl;
	//	cout << "iter_g: " << iter_g_out << endl;
	//	cout << "iter_f: " << iter_f_out << endl;
	//
	//	
	//	cout << "comp (real): \t \t \t";
	//	double realres[2];
	//	for (long i = 0; i < (sizeof(avec) / sizeof(double)); ++i) {
	//		if (avec[i] > bvec[i])
	//			realres[i] = 1;
	//		else
	//			realres[i] = 0;
	//
	//		cout << fixed << setprecision(precision) << realres[i] << ", ";
	//	}
	//	cout << endl;
	//
	//	cout << "New Comp (with Plain Algorithm): ";
	//	double res[2];
	//	for (long i = 0; i < 2; ++i) {		
	//		//NewComp(res[i], avec[i], bvec[i], iter, deg);
	//		NewCompG(res[i], avec[i], bvec[i], iter_f_out, iter_g_out, deg);
	//		cout << fixed << setprecision(precision) << res[i] << ", ";
	//		if ((abs(res[i] - realres[i]) > pow(2, -alpha)))
	//			cout << "error ";
	//	}
	//	cout << endl;
	//
	//	cout << "error (log): \t \t \t";
	//	for (long i = 0; i < (sizeof(avec) / sizeof(double)); ++i) 
	//		cout << fixed << setprecision(precision) << log(abs(res[i] - realres[i])) / log(2) << ", ";
	//	cout << endl;
	//
	//	cout << "\t \t \t alpha: \t " << alpha << endl;
	//	
	//	cout << "New Comp with g (with Plain Algorithm): ";
	//	for (long i = 0; i < 2; ++i) {
	//		NewCompG(res[i], avec[i], bvec[i], iter_f, iter_g, deg);
	//		cout << fixed << setprecision(precision) << res[i] << ", ";
	//		if ((abs(res[i] - realres[i]) > pow(2, -alpha)))
	//			cout << "error ";
	//	}
	//	cout << endl;
	//
	//	cout << "error (log): \t \t \t";
	//	for (long i = 0; i < (sizeof(avec) / sizeof(double)); ++i)
	//		cout << fixed << setprecision(precision) << log(abs(res[i] - realres[i])) / log(2) << ", ";
	//	cout << endl;
	//
	//	cout << "\t \t \t alpha: \t " << alpha << endl;
	//

//// Start Enc Algorithm
	 // Scheme & RotKey Gen
	TimeUtils timeutils;
	timeutils.start("Scheme Generation");
	Ring ring;
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);
	timeutils.stop("Scheme Generation");
	// scheme.addLeftRotKeys(secretKey);

	cout << "logN: " << logN << endl;
	cout << "logQ: " << logQ << endl;
	cout << "logp: " << logp << endl;

	cout << "-----------------------" << endl;

	long tot_err = 0;

	for (long test_iter = 0; test_iter < ctx_num; ++test_iter) {

		timeutils.start("Encrypt");
		Ciphertext ctxa, ctxb;
		scheme.encrypt(ctxa, &Avec[test_iter*num_slots], num_slots, logp, logq);
		scheme.encrypt(ctxb, &Bvec[test_iter*num_slots], num_slots, logp, logq);
		timeutils.stop("Encrypt");

		Ciphertext resa;

		timeutils.start("Evaluation");
		//evalNewComp(resa, ctxa, ctxb, iter, logp, scheme, secretKey);
		evalNewCompG(resa, ctxa, ctxb, iter_f, iter_g, logp, scheme, secretKey);
		timeutils.stop("Evaluation");

		Evec = scheme.decrypt(secretKey, resa);

		MaxErr = -100;
		cout << "New Comp G (with Enc Algorithm) Error Test: ";
		err_cnt = 0;
		for (long i = 0; i < num_slots; ++i) {
			//cout << Pvec[i] <<", "<< real(Evec[i]) << endl;
			Eerr[i] = abs(Rvec[test_iter*num_slots+i] - real(Evec[i]));
			if (Eerr[i] > prec) {
				//cout << "error!" << endl;
				err_cnt++;
				tot_err++;
			}
			Eerr[i] = log2(Eerr[i]);
			if (Eerr[i] > MaxErr) {
				MaxErr = Perr[i];
				argMax = i;
			}
		}
		if (err_cnt != 0) {
			cout << "Error occured: " << err_cnt++; cout << " at " << argMax << endl;
			cout << endl;
		}
		else
			cout << test_iter+1 <<"-th test / " << ctx_num <<" Done." << endl;
	}
	
	cout << "total error: " << tot_err << endl;


	return 0;
}
