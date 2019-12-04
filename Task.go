package main

import (
	"fmt"
	"math"

	"gonum.org/v1/gonum/mat"
)

func main() {
	// Physics params
	var (
		u1             = func(t float64) float64 { return math.Sin(0.01*t) + 273.15 }
		u0     float64 = 0. + 273.15
		lambda float64 = 42.1
		cc     float64 = 445.
		ro     float64 = 7800.
		gamma  float64 = 140.
		R      float64 = 0.5
		T      float64 = 60 * 1500
	)
	v1 := func(t float64) float64 { return (u1(t) - u0) / u0 }

	// Physics trans
	var (
		T1     = lambda / (cc * ro * R * R) * T
		gamma1 = R / lambda * gamma
	)

	// Numeric params
	var (
		N, M int     = 10, 4000
		sig  float64 = 0.5
	)

	var tao, h float64 = T1 / float64(M), 1 / float64(N)

	var x []float64
	for i := 0; i <= N; i++ {
		x = append(x, float64(i)*h)
	}

	var dx = make([]float64, N+1)
	dx[0] = (x[1]*x[1] - x[0]*x[0]) / (2 * h)
	dx[N] = (x[N]*x[N] - x[N-1]*x[N-1]) / (2 * h)

	for i := 1; i < N; i++ {
		dx[i] = (x[i+1]*x[i+1] - x[i-1]*x[i-1]) / (4 * h)
	}

	var dp = make([]float64, N+1)
	for i := 1; i < N; i++ {
		dp[i] = x[i] - h/2
	}

	var y = make([]float64, N+1)

	for i := 0; i <= N; i++ {
		y[i] = 0
	}

	printArr(y)
	fmt.Println()

	for j := 0; j < M; j++ {
		t := float64(j+1) * T / float64(M)
		b := make([]float64, N+1)
		c := make([]float64, N+1)
		d := make([]float64, N+1)
		phi := make([]float64, N+1)

		A := mat.NewDense(N+1, N+1, nil)
		yy := mat.NewDense(N+1, 1, nil)
		xx := mat.NewDense(N+1, 1, nil)

		b[0] = sig * tao / (h * h) * dp[1]
		c[0] = -dx[0]*0.5 - b[0]
		phi[0] = -0.5*dx[0]*y[0] - (1-sig)*tao/(h*h)*dp[1]*(y[1]-y[0])

		A.Set(0, 1, b[0])
		A.Set(0, 0, c[0])
		yy.Set(0, 0, phi[0])

		for i := 1; i < N; i++ {
			d[i] = -sig * tao / (h * h) * (x[i] - h/2)
			b[i] = -sig * tao / (h * h) * (x[i+1] - h/2)
			c[i] = x[i] - b[i] - d[i]
			phi[i] = x[i]*y[i] + tao*(1-sig)/(h*h)*(dp[i+1]*(y[i+1]-y[i])-dp[i]*(y[i]-y[i-1]))

			A.Set(i, i-1, d[i])
			A.Set(i, i, c[i])
			A.Set(i, i+1, b[i])
			yy.Set(i, 0, phi[i])
		}
		d[N] = sig * tao / (h * h) * dp[N]
		c[N] = -sig*tao/h*gamma1*x[N] - 0.5*dx[N] - d[N]
		phi[N] = (1-sig)*tao/h*gamma1*x[N]*y[N] - tao/h*x[N]*gamma1*v1(t) - 0.5*dx[N]*y[N] - (1-sig)*tao/(h*h)*dp[N]*(y[N]-y[N-1])

		A.Set(N, N, c[N])
		A.Set(N, N-1, d[N])
		yy.Set(N, 0, phi[N])

		xx.Solve(A, yy)

		// u^(j+1) <- u^j
		for i := range y {
			y[i] = xx.At(i, 0)
		}

		printArr(y)
		fmt.Println()
	}
}

func printArr(arr []float64) {
	var u0 float64 = 273.15
	for _, x := range arr {
		fmt.Printf("%14.9f", ((x+1)*u0 - 273.15))
		//if i != len(arr)-1 {
		//	fmt.Print("&")
		//} else {
		//	fmt.Print(" \\\\")
		//}
	}
}
