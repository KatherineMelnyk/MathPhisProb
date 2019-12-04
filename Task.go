package main

import (
	"fmt"
	//"math"

	"gonum.org/v1/gonum/mat"
)

func main() {
	// Physics params
	var (
		//u1             = func(t float64) float64 { return 0*math.Sin(0.01*t) + 273.15 }
		u1             = func(t float64) float64 { return 800 + 273.15 }
		u0     float64 = 0. + 273.15
		lambda float64 = 42.1
		cc     float64 = 445.
		ro     float64 = 7800.
		gamma  float64 = 140.
		R      float64 = 0.1
		T      float64 = 60 * 15
	)
	v1 := func(t float64) float64 { return (u1(t) - u0) / u0 } // обезрозмірена u

	// Physics trans
	var (
		T1     = lambda / (cc * ro * R * R) * T
		gamma1 = R / lambda * gamma
	)

	// Numeric params
	var (
		N, M int     = 50, 600
		sig  float64 = 1
	)

	var tao, h float64 = T1 / float64(M), 1 / float64(N)
	fmt.Printf("tao = %11.6f", tao)
	fmt.Print("\n")
	fmt.Printf("h = %11.6f", h)
	fmt.Print("\n")

	var x []float64
	for i := 0; i <= N; i++ {
		x = append(x, float64(i)*h)
	}

	var dx = make([]float64, N+1) // х хвилька m
	dx[0] = (x[1]*x[1]*x[1] - x[0]*x[0]*x[0]) / (3 * h)
	dx[N] = (x[N]*x[N]*x[N] - x[N-1]*x[N-1]*x[N-1]) / (3 * h)

	for i := 1; i < N; i++ {
		dx[i] = (x[i+1]*x[i+1]*x[i+1] - x[i-1]*x[i-1]*x[i-1]) / (6 * h)
	}

	var dp = make([]float64, N+1) // p хвилька
	for i := 1; i <= N; i++ {
		dp[i] = (x[i] - h/2) * (x[i] - h/2)
	}

	var y = make([]float64, N+1)

	for i := 0; i <= N; i++ {
		y[i] = 0
	}

	printArr(y)
	fmt.Println()

	for j := 0; j < M; j++ {// iterations by time
		t := float64(j+1) * T / float64(M) //moment t for current iteration
		b := make([]float64, N+1)
		c := make([]float64, N+1)
		d := make([]float64, N+1)
		phi := make([]float64, N+1)

		A := mat.NewDense(N+1, N+1, nil)
		yy := mat.NewDense(N+1, 1, nil)
		xx := mat.NewDense(N+1, 1, nil)

		b[0] = sig * tao / (h * h) * dp[1]
		c[0] = -dx[1]*0.5 - b[0]
		phi[0] = -0.5*dx[1]*y[0] - (1-sig)*tao/(h*h)*dp[1]*(y[1]-y[0])

		A.Set(0, 1, b[0])
		A.Set(0, 0, c[0])
		yy.Set(0, 0, phi[0])

		for i := 1; i < N; i++ {
			d[i] = -sig * tao / (h * h) * dp[i]
			b[i] = -sig * tao / (h * h) * dp[i+1]
			c[i] = dx[i] - b[i] - d[i]
			phi[i] = dx[i]*y[i] + tao*(1-sig)/(h*h)*(dp[i+1]*(y[i+1]-y[i])-dp[i]*(y[i]-y[i-1]))

			A.Set(i, i-1, d[i])
			A.Set(i, i, c[i])
			A.Set(i, i+1, b[i])
			yy.Set(i, 0, phi[i])
		}
		d[N] = sig * tao / (h * h) * dp[N]
		c[N] = -sig*tao/h*gamma1*x[N]*x[N] - 0.5*dx[N] - d[N]
		phi[N] = (1-sig)*tao/h*gamma1*x[N]*x[N]*y[N] - tao/h*x[N]*x[N]*gamma1*v1(t) - 0.5*dx[N]*y[N] - (1-sig)*tao/(h*h)*dp[N]*(y[N]-y[N-1])

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
		fmt.Printf("%11.6f", ((x+1)*u0 - 273.15))
	}
}
