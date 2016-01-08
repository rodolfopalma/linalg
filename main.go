package linalg

import (
	"fmt"
	"math"
)

const (
	NEARLY_EQUAL_TOLERANCE = 0.00001
)

type VectorStructure []float64

func NewVector(values []float64) VectorStructure {
	vector := make(VectorStructure, 0)
	return append(vector, values...)
}

func (v VectorStructure) Size() int {
	return len(v)
}

func (v VectorStructure) Norm(order float64) float64 {
	norm := 0.0
	switch {
	case order > 0:
		for _, value := range v {
			norm = norm + math.Pow(math.Abs(value), order)
		}
		norm = math.Pow(norm, 1/order)
	case order == 0:
		currentMax := math.Inf(-1)
		for _, value := range v {
			if value > currentMax {
				currentMax = value
			}
		}
		norm = currentMax
	}
	return norm
}

func (v VectorStructure) Dot(v2 VectorStructure) float64 {
	if v.Size() != v2.Size() {
		panic("Vectors should've the same size.")
	}
	sum := 0.0
	for i, value := range v {
		sum = sum + value*v2[i]
	}
	return sum
}

func (v VectorStructure) ScalarMul(lambda float64) VectorStructure {
	newValues := make([]float64, 0)
	for _, value := range v {
		newValues = append(newValues, lambda*value)
	}
	return NewVector(newValues)
}

func (v VectorStructure) Equal(v2 VectorStructure) bool {
	if v.Size() != v2.Size() {
		return false
	}
	for i, value := range v {
		if math.Abs(value-v2[i]) > NEARLY_EQUAL_TOLERANCE {
			return false
		}
	}
	return true
}

func (v VectorStructure) Sum(v2 VectorStructure) VectorStructure {
	if v.Size() != v2.Size() {
		panic("Vectors should've the same size.")
	}
	newValues := make([]float64, 0)
	for i, value := range v {
		newValues = append(newValues, value+v2[i])
	}
	return NewVector(newValues)
}

func (v VectorStructure) Minus(v2 VectorStructure) VectorStructure {
	return v.Sum(v2.ScalarMul(-1))
}

func (v VectorStructure) Insert(i int, value float64) VectorStructure {
	return NewVector(append(v[:i], append([]float64{value}, v[i:]...)...))
}

func (v VectorStructure) Remove(i int) VectorStructure {
	return NewVector(append(v[:i], v[i+1:]...))
}

type MatrixStructure [][]float64

func NewMatrix(values [][]float64) MatrixStructure {
	firstRowLength := len(values[0])
	matrix := make(MatrixStructure, 0)
	for _, row := range values {
		if firstRowLength == len(row) {
			matrix = append(matrix, row)
		} else {
			panic("Matrix shape is not valid.")
		}

	}
	return matrix
}

func Eye(n int) MatrixStructure {
	matrix := make(MatrixStructure, 0)
	for i := 0; i < n; i++ {
		row := make([]float64, n)
		row[i] = 1
		matrix = append(matrix, row)
	}
	return matrix
}

func (m MatrixStructure) Copy() MatrixStructure {
	matrix := make(MatrixStructure, 0)
	shape := m.Shape()
	for i := 0; i < shape[0]; i++ {
		row := make([]float64, 0)
		row = append(row, m[i]...)
		matrix = append(matrix, row)
	}
	return matrix
}

func (m MatrixStructure) Shape() []int {
	return []int{len(m), len(m[0])}
}

func (m MatrixStructure) Transpose() MatrixStructure {
	matrixValues := make([][]float64, 0)
	origShape := m.Shape()
	for i := 0; i < origShape[1]; i++ {
		rowValues := make([]float64, 0)
		for j := 0; j < origShape[0]; j++ {
			rowValues = append(rowValues, m[j][i])
		}
		matrixValues = append(matrixValues, rowValues)
	}
	// Currently doesn't make sense to call NewMatrix, but probably in the
	// short term the MatrixStructure will get more complex.
	return NewMatrix(matrixValues)
}

func (m MatrixStructure) Mul(m2 MatrixStructure) MatrixStructure {
	leftShape := m.Shape()
	rightShape := m2.Shape()
	if leftShape[1] != rightShape[0] {
		panic("Matrices shapes are not valid.")
	}
	matrixValues := make([][]float64, 0)
	for i := 0; i < leftShape[0]; i++ {
		rowValues := make([]float64, 0)
		for j := 0; j < rightShape[1]; j++ {
			sumValue := 0.0
			for k := 0; k < leftShape[1]; k++ {
				sumValue = sumValue + m[i][k]*m2[k][j]
			}
			rowValues = append(rowValues, sumValue)
		}
		matrixValues = append(matrixValues, rowValues)
	}
	return NewMatrix(matrixValues)
}

func (m MatrixStructure) Map(f func(float64) float64) MatrixStructure {
	matrixValues := make([][]float64, 0)
	origShape := m.Shape()
	for i := 0; i < origShape[0]; i++ {
		rowValues := make([]float64, 0)
		for j := 0; j < origShape[1]; j++ {
			rowValues = append(rowValues, f(m[i][j]))
		}
		matrixValues = append(matrixValues, rowValues)
	}
	return NewMatrix(matrixValues)
}

func (m MatrixStructure) LUDecomposition() (MatrixStructure, MatrixStructure) {
	shape := m.Shape()
	L := Eye(shape[0])
	U := m.Copy()
	for !U.IsUpperDiagonal() {
		for col := 0; col < shape[1]-1; col++ {
			for row := col + 1; row < shape[0]; row++ {
				value := U[row][col]
				if value != 0 {
					factor := value / U[col][col]
					L[row][col] = factor
					U = U.RowOperation(row, col, -factor)
				}
			}
		}
	}
	return L, U
}

func (m MatrixStructure) IsLowerDiagonal() bool {
	shape := m.Shape()
	for i := 0; i < shape[0]; i++ {
		for j := i + 1; j < shape[1]; j++ {
			if m[i][j] != 0 {
				return false
			}
		}
	}
	return true
}

func (m MatrixStructure) IsUpperDiagonal() bool {
	return m.Transpose().IsLowerDiagonal()
}

func (m MatrixStructure) IsDiagonal() bool {
	return m.IsLowerDiagonal() || m.IsUpperDiagonal()
}

func (m MatrixStructure) SolveDiagonalSystem(v VectorStructure) VectorStructure {
	// Upper or lower diagonal?
	// TODO: assuming lower diagonal matrix, make this work in a flexible way.
	matrixShape := m.Shape()
	vectorSize := v.Size()
	if matrixShape[0] != vectorSize {
		panic("Vector size and matrix shape aren't consistent.")
	}
	// Create slice of the values of the solution vector, initializing it with
	// the first diagonal element of matrix.
	// TODO: DRY!!
	solutionValues := make([]float64, matrixShape[0])
	switch {
	case m.IsLowerDiagonal():
		solutionValues[0] = v[0] / m[0][0]
		for row := 1; row < matrixShape[0]; row++ {
			value := v[row]
			for col := 0; col < row; col++ {
				value = value - m[row][col]*solutionValues[col]
			}
			solutionValues[row] = value / m[row][row]
		}
	case m.IsUpperDiagonal():
		solutionValues[matrixShape[0]-1] = v[matrixShape[0]-1] / m[matrixShape[0]-1][matrixShape[0]-1]
		for row := matrixShape[0] - 2; row >= 0; row-- {
			value := v[row]
			for col := row + 1; col < matrixShape[1]; col++ {
				value = value - m[row][col]*solutionValues[col]
			}
			solutionValues[row] = value / m[row][row]
		}
	default:
		panic("Not a diagonal system.")
	}
	return NewVector(solutionValues)
}

func (m MatrixStructure) SolveSystem(v VectorStructure) VectorStructure {
	if m.IsDiagonal() {
		return m.SolveDiagonalSystem(v)
	} else {
		L, U := m.LUDecomposition()
		y := L.SolveDiagonalSystem(v)
		return U.SolveDiagonalSystem(y)
	}
}

// Performs the following row operation:
// m[row] = m[row] + factor * m[otherRow]
func (m MatrixStructure) RowOperation(row int, otherRow int, factor float64) MatrixStructure {
	E := Eye(m.Shape()[0])
	E[row][otherRow] = factor
	return E.Mul(m)
}

func (m MatrixStructure) Equal(m2 MatrixStructure) bool {
	shape_1 := m.Shape()
	shape_2 := m2.Shape()
	if shape_1[0] != shape_2[0] || shape_1[1] != shape_2[1] {
		return false
	}
	for i, row := range m {
		for j, value := range row {
			if math.Abs(value-m2[i][j]) > NEARLY_EQUAL_TOLERANCE {
				return false
			}
		}
	}
	return true
}

func (m MatrixStructure) String() string {
	str := "-\n"
	for i, row := range m {
		if i > 0 {
			str += "\n"
		}
		str += "|"
		for j, value := range row {
			str += fmt.Sprintf("%0.2f", value)
			if j != len(row)-1 {
				str += " "
			}
		}
	}
	return str + "\n"
}
