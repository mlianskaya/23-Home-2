#include <iostream>
#include <fstream>

using namespace std;

template<typename T>
class Matrix
{
private:
	T **A;//матрица
	int m; //количество строк
	int n;//количество столбцов
public:
	//конструктор
	Matrix()
	{
		m = 0;
		n = 0;
		A = nullptr;
	}
	Matrix(int x, int y) : m(x), n(y)
	{
		A = new T*[m];
		for (int i = 0; i < m; i++) {
			A[i] = new T[n];
		}
	}
	Matrix(const Matrix& matr)
	{
		m = matr.m;
		n = matr.n;
		A = new T*[m];
		for (int i = 0; i < m; i++) {
			A[i] = new T[n];
			for (int j = 0; j < n; j++) {
				A[i][j] = matr.A[i][j];
			}
		}
	}
	Matrix(const string& filename) {
		ifstream file(filename);
		if (file.is_open()) {
			file >> m >> n;
			A = new T*[m];
			for (int i = 0; i < m; i++) {
				A[i] = new T[n];
				for (int j = 0; j < n; j++) {
					file >> A[i][j];
				}
			}
			file.close();
		}
		else {
			cerr << "Unable to open file\n ";
		}
	}


	friend ostream& operator<<(ostream& os, const Matrix& matrix) {
		for (int i = 0; i < matrix.m; ++i) {
			for (int j = 0; j < matrix.n; ++j)
				os << matrix.A[i][j] << " ";
			os << endl;
		}
		return os;
	}

	friend istream& operator>>(istream& is, Matrix& matrix) {
		cout << "enter number of rows and columns of the matrix:\n";
		is >> matrix.m >> matrix.n;
		matrix.A = new T*[matrix.m];
		for (int i = 0; i < matrix.m; ++i) {
			matrix.A[i] = new T[matrix.n];
			for (int j = 0; j < matrix.n; ++j)
			{
				cout << "Enter the number: row - " << i + 1 << "; column - " << j + 1 << ": ";
				is >> matrix.A[i][j];
			}
		}
		return is;
	}

	friend ifstream& operator>>(ifstream& filename, Matrix& matr) {
		filename >> matr.m;
		filename >> matr.n;
		matr.matrix = new T*[matr.m];
		for (int i = 0; i < matr.m; i++) {
			matr.matrix[i] = new T[matr.n];
			for (int j = 0; j < matr.n; j++) {
				filename >> matr.A[i][j];
			}
		}
		return filename;
	}

	friend ofstream& operator<<(ofstream& filename, const Matrix& matr) {
		for (int i = 0; i < matr.m; i++) {
			for (int j = 0; j < matr.n; j++) {
				filename << matr.A[i][j] << " ";
			}
			filename << endl;
		}
		return filename;
	}

	static Matrix zero(int numM, int numN) {
		Matrix zeroMatrix(numM, numN);
		for (int i = 0; i < numM; ++i) {
			for (int j = 0; j < numN; ++j) {
				zeroMatrix.A[i][j] = 0;
			}
		}
		return zeroMatrix;
	}

	static Matrix identity(int size) {
		Matrix identityMatrix(size, size);
		for (int i = 0; i < size; ++i) {
			for (int j = 0; j < size; ++j) {
				if (i == j) {
					identityMatrix.A[i][j] = 1;
				}
				else {
					identityMatrix.A[i][j] = 0;
				}
			}
		}
		return identityMatrix;
	}

	// перегрузка оператора * 
	Matrix operator * (const Matrix &matr) const
	{
		if (n != matr.m)
		{
			cerr << "impossible to multiply!\n";
			return Matrix();
		}
		Matrix result(m, matr.n);
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < matr.n; j++) {
				result.A[i][j] = 0;
				for (int k = 0; k < n; k++) {
					result.A[i][j] += A[i][k] * matr.A[k][j];
				}
			}
		}
		return result;
	}
	Matrix operator * (const T& scalar) const
	{
		Matrix result(m, n);
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				result.A[i][j] = A[i][j] * scalar;
			}
		}
		return result;
	}
	//перегрузка оператора +
	Matrix operator + (const Matrix &matr) const
	{
		if (m != matr.m || n != matr.n)
		{
			cerr << "impossible to add!\n";
			return Matrix();
		}
		Matrix result(m, n);
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				result.A[i][j] = A[i][j] + matr.A[i][j];
			}
		}
		return result;
	}
	//перегрузка оператора -
	Matrix operator - (const Matrix &matr) const
	{
		if (m != matr.m || n != matr.n)
		{
			cerr << "impossible to substract!\n";
			return Matrix();
		}
		Matrix result(m, n);
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				result.A[i][j] = A[i][j] - matr.A[i][j];
			}
		}
		return result;
	}

	Matrix& operator=(const Matrix& matr) {
		if (this == &matr) {
			return *this;
		}
		if (A != nullptr) {
			for (int i = 0; i < m; ++i) {
				delete[] A[i];
			}
			delete[] A;
		}
		m = matr.m;
		n = matr.n;
		A = new T*[m];
		for (int i = 0; i < m; ++i) {
			A[i] = new T[n];
			for (int j = 0; j < n; ++j) {
				A[i][j] = matr.A[i][j];
			}
		}
		return *this;
	}


	double determinant() {
		if (m != n) {
			throw "Unable to calculate determinant";
			return 0;
		}
		if (m == 1) {
			return A[0][0];
		}

		if (m == 2) {
			return A[0][0] * A[1][1] - A[0][1] * A[1][0];
		}
		double det = 0;
		Matrix minor(m - 1, n - 1);
		for (int j = 0; j < n; j++) {
			for (int i = 1; i < m; i++) {
				int k = 0;
				for (int l = 0; l < n; l++) {
					if (l != j) {
						minor.A[i - 1][k] = A[i][l];
						k++;
					}
				}
			}
			det += pow(-1, j)*A[0][j] * minor.determinant();
		}
		return det;
	}

	double cofactor(int row, int column) {
		Matrix minor(m - 1, n - 1);
		int m1, n1, sign;
		m1 = 0;
		for (int i = 0; i < m - 1; i++) {
			if (i == row - 1) {
				m1 = 1;
			}
			n1 = 0;
			for (int j = 0; j < m - 1; j++) {
				if (j == column - 1) {
					n1 = 1;
				}
				minor.A[i][j] = A[i + m1][j + n1];
			}
		}
		if ((row + column) % 2 == 0) {
			sign = 1;
		}
		else {
			sign = -1;
		}
		return sign * minor.determinant();
	}

	Matrix cofactorMatrix() {
		Matrix copy(m, n);
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				copy.A[i][j] = A[i][j];
			}
		}
		Matrix cofactorM(m, n);
		if (m != n) {
			throw "Unable to find the cofactor matrix!";
		}
		else {
			for (int i = 0; i < m; i++) {
				for (int j = 0; j < n; j++) {
					cofactorM.A[i][j] = copy.cofactor(i + 1, j + 1);
				}
			}
		}
		return cofactorM;
	}

	Matrix transposeMatrix() {
		Matrix transpose(n, m);

		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < m; ++j) {
				transpose.A[i][j] = A[j][i];
			}
		}

		return transpose;
	}

	Matrix operator !()
	{
		Matrix copy(m, n);
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				copy.A[i][j] = A[i][j];
			}
		}
		double det = copy.determinant();
		if ((m != n) || (det == 0)) {
			throw "Unable to find the inverse matrix";
		}
		else {
			Matrix cofactorMat = copy.cofactorMatrix();
			Matrix adjointMat = cofactorMat.transposeMatrix();
			Matrix inverse = adjointMat * (1 / det);
			return inverse;
		}
	}


	~Matrix() {
		if (A != nullptr)
		{
			for (int i = 0; i < m; i++) {
				delete[] A[i];
			}
			delete[] A;
		}
	}
};

int main()
{
	Matrix<int> A;
	cin >> A;
	cout << "matrix 1:\n";
	cout << A;
	cout << endl;

	Matrix<int> B("input.txt");
	cout << "matrix 2:\n";
	cout << B;
	cout << endl;

	cout << "matrix addition result:\n";
	Matrix<int> C = A + B;
	cout << C;
	cout << endl;

	cout << "matrix substraction result:\n";
	Matrix<int> D = A - B;
	cout << D;
	cout << endl;

	cout << "matrix multiplication result:\n";
	Matrix<int> E = A * B;
	cout << E;
	cout << endl;

	int scal;
	cout << "enter the number to multiply the matrix by the scalar:\n";
	cin >> scal;
	Matrix<int> F = A * scal;
	cout << "the result of multiplying a matrix by a number:\n";
	cout << F;
	cout << endl;

	cout << "the result of assigning a matrix of certain values:\n";
	Matrix<int> G = F;
	cout << G;
	cout << endl;

	cout << "zero matrix\n";
	int i, j;
	cout << "enter number of rows:\n";
	cin >> i;
	cout << "enter number of columns:\n";
	cin >> j;
	Matrix<int> zeroMatrix = Matrix<int>::zero(i, j);
	cout << zeroMatrix;
	cout << endl;

	cout << "identity matrix\n";
	int k;
	cout << "enter number of rows and columns:\n";
	cin >> k;
	Matrix<double> identityMatrix = Matrix<double>::identity(k);
	cout << identityMatrix;
	double scal2;
	cout << "enter the number to multiply the matrix by the scalar:\n";
	cin >> scal2;
	Matrix<double> P = identityMatrix * scal2;
	ofstream file("output.txt");
	file << P;
	cout << "the result is written to the file\n";
	cout << endl;

	Matrix<double> N;
	cin >> N;
	try {
		Matrix<double> I = !N;
		cout << "inverse matrix:\n";
		cout << I;
	}
	catch (const char* error_message) {
		cout << error_message << endl;
	}
	return 0;
};