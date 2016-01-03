package linalg

import (
	"flag"
	"fmt"
	"os"
	"testing"
)

var v VectorStructure
var m MatrixStructure

func TestMain(main *testing.M) {
	flag.Parse()
	v = NewVector([]float64{1, 1, 1, 1})
	m = NewMatrix([][]float64{
		{1, 2, 3, 4},
		{5, 6, 7, 8},
		{9, 10, 11, 12},
	})
	os.Exit(main.Run())
}

func TestVectorSize(t *testing.T) {
	t.Log("Calculating vector size of [1, 1, 1, 1]... (expected size: 4)")
	if size := v.Size(); size != 4 {
		t.Errorf("Expected vector size of 4, but it was %0.2f...", size)
	}
}

func TestVectorNorm(t *testing.T) {
	t.Log("Calculating vector norm of [1, 1, 1, 1]... (expected norm: 2)")
	if norm := v.Norm(2); norm != 2 {
		t.Errorf("Expected vector norm of 2, but it was %0.2f...", norm)
	}

}

func TestDotProduct(t *testing.T) {
	t.Log("Calculating dot product between [1, 1, 1, 1] and [0, 1, 0, 1]... (expected value: 2)")
	v2 := NewVector([]float64{0, 1, 0, 1})
	if dot := v.Dot(v2); dot != 2 {
		t.Errorf("Expected dot product of 2, but it was %0.2f...", dot)
	}
}

func TestVectorScalarMul(t *testing.T) {
	t.Log("Calculating the element-wise scalar multiplication of [1, 1, 1, 1] by 3... (expected vector: [3, 3, 3, 3])")
	v2 := v.ScalarMul(3)
	if result := NewVector([]float64{3, 3, 3, 3}); !result.Equal(v2) {
		t.Errorf("Expected vector [3, 3, 3, 3], but it was %f...", v2)
	}
}

func TestVectorEqual(t *testing.T) {
	t.Log("Testing for vector equality...")
	if v2 := NewVector([]float64{1, 1, 1, 1}); !v.Equal(v2) {
		t.Errorf("Expected true for equality, but it was false...")
	}

}

func TestVectorSum(t *testing.T) {
	t.Log("Calculating the sum between [1, 1, 1, 1] and [4, 4, 4, 4]... (expected vector: [5, 5, 5, 5])")
	v2 := NewVector([]float64{4, 4, 4, 4})
	result := NewVector([]float64{5, 5, 5, 5})
	if sum := v.Sum(v2); !result.Equal(sum) {
		t.Errorf("Expected vector [5, 5, 5, 5], but it was %f...", sum)
	}
}

func TestVectorMinus(t *testing.T) {
	t.Log("Calculating the difference between [1, 1, 1, 1] and [4, 4, 4, 4]... (expected vector: [-3, -3, -3, -3])")
	v2 := NewVector([]float64{4, 4, 4, 4})
	result := NewVector([]float64{-3, -3, -3, -3})
	if difference := v.Minus(v2); !result.Equal(difference) {
		t.Errorf("Expected vector [-3, -3, -3, -3], but it was %f...", difference)
	}
}

func TestVectorInsert(t *testing.T) {
	t.Log("Inserting a single 0 in the 2nd-index position of v...")
	v2 := v.Insert(2, 0)
	if v2[2] != 0 {
		t.Errorf("Expected 2nd-index of v2 to be 0, but it was %0.2f...", v2[2])
	}
	if size := v2.Size(); size != 5 {
		t.Errorf("Expected a new vector with size equal to 5, but it was %d...", size)
	}
}

func TestVectorRemove(t *testing.T) {
	t.Log("Removing first element of v... (expected new size: 3)")
	v2 := v.Remove(0)
	if size := v2.Size(); size != 3 {
		t.Errorf("Expected a new vector with size equal to 3, but it was %d...", size)
	}
}

func TestMatrixShape(t *testing.T) {
	t.Log("Calculating matrix shape of m... (expected shape: 3x4)")
	if shape := m.Shape(); shape[0] != 3 || shape[1] != 4 {
		t.Errorf("Expected shape of m to be 3x4, but it was %d...", shape)
	}
}

// TODO: Test transposed matrix in an element-wise way.
func TestMatrixTranspose(t *testing.T) {
	t.Log("Transposing matrix m...")
	mT := m.Transpose()
	if shape := mT.Shape(); shape[0] != 4 || shape[1] != 3 {
		t.Errorf("Expected shape of m.Transpose() to be 4x3, but it was %d...", shape)
	}
}

func TestMatrixMultiplication(t *testing.T) {
	t.Log("Multiplicating matrix m...")
	m2 := NewMatrix([][]float64{
		{1, 0},
		{0, 1},
		{0, 0},
		{0, 0},
	})
	result := m.Mul(m2)
	if shape := result.Shape(); shape[0] != 3 || shape[1] != 2 {
		t.Errorf("Expected shape of multiplied matrix to be 3x2, but it was %d...", shape)
	}
}

func TestMatrixMap(t *testing.T) {

}

func TestMatrixLUDecomposition(t *testing.T) {
	// L, U := m.Mul(m.Transpose()).LUDecomposition()
	// fmt.Println(m.Mul(m.Transpose()), L, U, L.Mul(U))
}

func TestSolveDiagonalSystem(t *testing.T) {
	// testM := Eye(3)
	// testM[0][1] = 3
	// testV := NewVector([]float64{1, 2, 4})
	// fmt.Println(testM.SolveDiagonalSystem(testV))
}

func TestSolveSystem(t *testing.T) {
	testM := NewMatrix([][]float64{
		{1, 1, 1},
		{0, 2, 5},
		{2, 5, -1},
	})
	testV := NewVector([]float64{6, -4, 27})
	fmt.Println(testM.SolveSystem(testV))
}

func TestMatrixIsUpperDiagonal(t *testing.T) {
	// temp := NewMatrix([][]float64{
	// {1, 0, 0, 0, 0},
	// {3, 4, 0, 0, 0},
	// {0, 2, 3, 0, 0},
	// {3, 2, 1, 4, 0},
	// })
	// fmt.Println(temp.IsUpperDiagonal())
}

func TestMatrixEye(t *testing.T) {
	// fmt.Println(Eye(5))
}

func TestMatrixRowOperation(t *testing.T) {
	// temp := Eye(5)
	// fmt.Println(temp.RowOperation(2, 1, 4))
}
