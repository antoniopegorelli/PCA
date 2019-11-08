#include <iostream>

using namespace std;

// Matrix class declaration
class matrix
{
private:
	int matrixRows;
	int matrixCols;
	double matrixData[100][100];
public:
	matrix(int, int);
	~matrix();
	void addData(double, int, int);
	double getData(int, int);
	int getRows(void);
	int getCols(void);
	void print(void);

};

// Matrix class contructor
matrix::matrix(int rows, int cols)
{
	matrixRows = rows;
	matrixCols = cols;
	matrixData[0][0] = 0;
}

// Matrix class destructor
matrix::~matrix()
{

}

// Matrix class to add data
void matrix::addData(double data, int row, int col)
{
	matrixData[row][col] = data;
}

// Matrix class to get data
double matrix::getData(int row, int col)
{
	return matrixData[row][col];
}

// Matrix class to get row size
int matrix::getRows(void)
{
	return matrixRows;
}

// Matrix class to get column size
int matrix::getCols(void)
{
	return matrixCols;
}

// Matrix class to print it on screen
void matrix::print(void)
{
	cout << "Matrix " << matrixRows << " x " << matrixCols << ":" << endl << endl;
	
	for (int i = 0; i < this->getRows(); i++)
	{
		for (int j = 0; j < this->getCols(); j++)
		{
			cout << this->getData(i, j) << "    ";
		}
		cout << endl << endl;
	}
	cout << endl;
}

// Alps Water data
double alps[17][2] = { {194.5, 20.79},
{194.3, 20.79},
{197.9, 22.4},
{198.4, 22.67},
{199.4, 23.15},
{199.9, 23.35},
{200.9, 23.89},
{201.1, 23.99},
{201.4, 24.02},
{201.3, 24.01},
{203.6, 25.14},
{204.6, 26.57},
{209.5, 28.49},
{208.6, 27.76},
{210.7, 29.04},
{211.9, 29.88},
{212.2, 30.06} };

// Books x Grades data
double books[40][3] = { {0,	9,	45},
{1,	15,	57},
{0,	10,	45},
{2,	16,	51},
{4,	10,	65},
{4,	20,	88},
{1,	11,	44},
{4,	20,	87},
{3,	15,	89},
{0,	15,	59},
{2,	8,	66},
{1,	13,	65},
{4,	18,	56},
{1,	10,	47},
{0,	8,	66},
{1,	10,	41},
{3,	16,	56},
{0,	11,	37},
{1,	19,	45},
{4,	12,	58},
{4,	11,	47},
{0,	19,	64},
{2,	15,	97},
{3,	15,	55},
{1,	20,	51},
{0,	6,	61},
{3,	15,	69},
{3,	19,	79},
{2,	14,	71},
{2,	13,	62},
{3,	17,	87},
{2,	20,	54},
{2,	11,	43},
{3,	20,	92},
{4,	20,	83},
{4,	20,	94},
{3,	9,	60},
{1,	8,	56},
{2,	16,	88},
{0,	10,	62} };

// US Census Dataset data
double census[11][2] = { {1900,	75.9950},
{1910,	91.9720},
{1920,	105.7110},
{1930,	123.2030},
{1940,	131.6690},
{1950,	150.6970},
{1960,	179.3230},
{1970,	203.2120},
{1980,	226.5050},
{1990,	249.6330},
{2000,	281.4220} };

// Matrix transposer
matrix* matTrans(matrix* mat)
{
	matrix* result = new matrix(mat->getCols(), mat->getRows());

	for (int i = 0; i < mat->getRows(); i++)
	{
		for (int j = 0; j < mat->getCols(); j++)
		{
			result->addData(mat->getData(i, j), j, i);
		}
	}
	return result;
}

// Matrix multiplier
matrix* matMult(matrix* matA, matrix* matB)
{

	if (matA->getCols() != matB->getRows())
	{
		cout << "Matrix multiplication incompatibility" << endl;
	}

	matrix* result = new matrix(matA->getRows(), matB->getCols());
	double total;

	for (int i = 0; i < matA->getRows(); i++)
	{
		for (int j = 0; j < matB->getCols(); j++)
		{
			total = 0;
			for (int a = 0; a < matA->getCols(); a++)
			{
				total += matA->getData(i, a) * matB->getData(a, j);
			}
			result->addData(total, i, j);
		}
	}
	return result;
}

// Matrix minor build
matrix* matMinor(matrix* mat, int row, int col)
{
	matrix* calc = new matrix(mat->getRows() - 1, mat->getCols() - 1);

	for (int i = 0; i < mat->getRows(); i++)
	{
		for (int j = 0; j < mat->getCols(); j++)
		{
			if (i != row && j != col)
			{
				int a = i;
				int b = j;
				if (i > row)
				{
					a--;
				}
				if (j > col)
				{
					b--;
				}
				calc->addData(mat->getData(i, j), a, b);
			}
		}
	}
	return calc;
}

// Matrix determinant calculator
double matDet(matrix* mat)
{
	double det = 0;

	//for (int i = 0; i < mat->getRows(); i++)
	//{
	//	for (int j = 0; j < mat->getCols(); j++)
	//	{
	//		cout << mat->getData(i, j) << "    ";
	//	}
	//	cout << endl << endl;
	//}

	if (mat->getRows() == 2)
	{
		det = (mat->getData(0, 0) * mat->getData(1, 1)) - (mat->getData(0, 1) * mat->getData(1, 0));
	}
	else
	{
		for (int a = 0; a < mat->getCols(); a++)
		{
			matrix* calc = matMinor(mat, 0, a);

			if (a % 2 == 0)
			{
				det += mat->getData(0, a) * (matDet(calc));
			}
			else
			{
				det -= mat->getData(0, a) * (matDet(calc));
			}
		}
	}
	//cout << "determinant = " << det << endl << endl;
	return det;
}

// Matrix inverter
matrix* matInv(matrix* mat)
{
	double det = matDet(mat);
	matrix* calc = new matrix(mat->getRows(), mat->getCols());

	if (mat->getRows() == 2)
	{
		calc->addData(mat->getData(1, 1) / det, 0, 0);
		calc->addData(-mat->getData(0, 1) / det, 0, 1);
		calc->addData(-mat->getData(1, 0) / det, 1, 0);
		calc->addData(mat->getData(0, 0) / det, 1, 1);
	}
	else
	{
		matrix* cofactors = new matrix(mat->getRows(), mat->getCols());
		matrix* adjugate = new matrix(mat->getRows(), mat->getCols());

		for (int i = 0; i < mat->getRows(); i++)
		{
			for (int j = 0; j < mat->getCols(); j++)
			{
				if (((i % 2 == 0) && (j % 2 == 0)) || ((i % 2 != 0) && (j % 2 != 0)))
				{
					cofactors->addData(matDet(matMinor(mat, i, j)), i, j);
				}
				else
				{
					cofactors->addData(-matDet(matMinor(mat, i, j)), i, j);
				}
			}
		}

		for (int i = 0; i < cofactors->getRows(); i++)
		{
			for (int j = i; j < cofactors->getCols(); j++)
			{
				if (i == j)
				{
					adjugate->addData(cofactors->getData(i, j), i, j);
				}
				else
				{
					adjugate->addData(cofactors->getData(i, j), j, i);
					adjugate->addData(cofactors->getData(j, i), i, j);
				}
			}
		}

		for (int i = 0; i < adjugate->getRows(); i++)
		{
			for (int j = 0; j < adjugate->getCols(); j++)
			{
				calc->addData(adjugate->getData(i, j) / det, i, j);
			}
		}
	}
	return calc;
}

// Linear regression
matrix* linear(matrix* x, matrix* y)
{
	cout << "Linear:" << endl << endl;

	matrix* xt = matTrans(x);
	matrix* xtx = matMult(xt, x);
	matrix* xty = matMult(xt, y);
	matrix* xtx1 = matInv(xtx);
	matrix* beta = matMult(xtx1, xty);

	beta->print();

	return beta;
}

// Quadratic regression
matrix* quad(matrix* x, matrix* y)
{
	cout << "Quadratic:" << endl << endl;

	matrix* x2 = new matrix(x->getRows(), (x->getCols() * 2) - 1);

	for (int i = 0; i < x->getRows(); i++)
	{
		x2->addData(1, i, 0);
		for (int j = 1; j < ((x->getCols() * 2) - 1); j++)
		{
			if (j < x->getCols())
			{
				x2->addData(x->getData(i, j), i, j);
			}
			else
			{
				x2->addData(pow(x->getData(i, j - x->getCols() + 1), 2), i, j);
			}
		}
	}

	matrix* xt = matTrans(x2);
	matrix* xtx = matMult(xt, x2);
	matrix* xty = matMult(xt, y);
	matrix* xtx1 = matInv(xtx);
	matrix* beta = matMult(xtx1, xty);

	beta->print();

	return beta;
}

// Robust regression
matrix* robust(matrix* x, matrix* y, matrix* ref)
{
	cout << "Robust:" << endl << endl;

	matrix* w = new matrix(x->getRows(), 1);
	for (int i = 0; i < x->getRows(); i++)
	{
		w->addData(abs(1 / (y->getData(i, 0) - (ref->getData(0, 0) + x->getData(i, 1) * ref->getData(1, 0)))), i, 0);
	}
	//w->print();

	matrix* wx = new matrix(x->getRows(), x->getCols());
	for (int i = 0; i < x->getRows(); i++)
	{
		for (int j = 0; j < x->getCols(); j++)
		{
			wx->addData(x->getData(i, j) * w->getData(i, 0), i, j);
		}

	}

	matrix* wy = new matrix(y->getRows(), y->getCols());
	for (int i = 0; i < y->getRows(); i++)
	{
		wy->addData(y->getData(i, 0) * w->getData(i, 0), i, 0);
	}

	matrix* xt = matTrans(x);
	matrix* xtx = matMult(xt, wx);
	matrix* xty = matMult(xt, wy);
	matrix* xtx1 = matInv(xtx);
	matrix* beta = matMult(xtx1, xty);

	beta->print();

	return beta;
}

// Covariance calculation
matrix* covariance(matrix* ref)
{
	matrix* cov = new matrix(ref->getCols(), ref->getCols());
	//ref->print();
	// Averages calculation
	matrix* avg = new matrix(ref->getCols(), 1);
	for (int i = 0; i < avg->getRows(); i++)
	{
		double sum = 0;
		for (int j = 0; j < ref->getRows(); j++)
		{
			sum += ref->getData(j, i);
		}
		double mean = sum / ref->getRows();
		avg->addData(mean, i, 0);
	}

	avg->print();

	// Covariance calculation
	for (int i = 0; i < ref->getRows(); i++)
	{
		for (int j = i; j < ref->getCols(); j++)
		{
			double sum = 0;
			for (int a = 0; a < ref->getRows(); a++)
			{
				sum += (ref->getData(a, i) - avg->getData(i, 0)) * (ref->getData(a, j) - avg->getData(j, 0));
			}
			cov->addData(sum / double(ref->getRows() - (double)1), i, j);
			if (i != j) { cov->addData(cov->getData(i, j), j, i); }
		}
	}
	cout << "Covariance: " << endl;
	cov->print();
	return cov;
}

// Eigenvalue calculation
matrix* eigenvalues(matrix* ref)
{
	matrix* eigenval;
	
	if (ref->getRows() == 2)
	{
		// Based on 2x2 determinant (y as lambda): y^2 - (a + d)*y + (a*d - c*b)
		// For matrix:
		// |a  b|
		// |c  d|
		double a = 1;
		double b = -(ref->getData(0, 0) + ref->getData(1, 1));
		double c = (ref->getData(0, 0) * ref->getData(1, 1)) - (ref->getData(1, 0) * ref->getData(0, 1));

		cout << "a = " << a << " b = " << b << " c = " << c << endl << endl;
		//cout << (-b + sqrt(pow(b, 2) - 4 * a * c)) / (2 * a) << endl;

		eigenval = new matrix(2, 1);

		eigenval->addData((-b + sqrt(pow(b, 2) - (4 * a * c))) / (2 * a), 0, 0);
		eigenval->addData((-b - sqrt(pow(b, 2) - (4 * a * c))) / (2 * a), 1, 0);
	}
	else
	{
		//ref->print();
		// Based on 3x3 determinant (y as lambda): -y3 +(i +a +e)y2 +(-a*i -e*i -a*e +g*c +h*f +d*b)*y +(a*e*i +b*f*g -g*e*c -h*f*a -i*d*b)
		// For matrix:
		// |a  b  c|
		// |d  e  f|
		// |g  h  i|
		double a = -1;
		double b = (ref->getData(2, 2) + ref->getData(0, 0) + ref->getData(1, 1));
		double c = (-(ref->getData(0, 0) * ref->getData(2, 2)) - (ref->getData(1, 1) * ref->getData(2, 2)) - (ref->getData(0, 0) * ref->getData(1, 1)) + (ref->getData(2, 0) * ref->getData(0, 2)) + (ref->getData(2, 1) * ref->getData(1, 2)) + (ref->getData(1, 0) * ref->getData(0, 1)));
		double d = ((ref->getData(0, 0) * ref->getData(1, 1) * ref->getData(2, 2)) + (ref->getData(0, 1) * ref->getData(1, 2) * ref->getData(2, 0)) + (ref->getData(0, 2) * ref->getData(1, 0) * ref->getData(2, 1)) - (ref->getData(2, 0) * ref->getData(1, 1) * ref->getData(0, 2)) - (ref->getData(2, 1) * ref->getData(1, 2) * ref->getData(0, 0)) - (ref->getData(2, 2) * ref->getData(1, 0) * ref->getData(0, 1)));
		
		cout << "a = " << a << " b = " << b << " c = " << c << " d = " << d << endl;

		double discriminant = 18 * a * b * c * d - 4 * pow(b, 3) * d + pow(b, 2) * pow(c, 2) - 4 * a * pow(c, 3) - 27 * pow(a, 2) * pow(d, 2);

		cout << "disc = " << discriminant << endl << endl;

		double disc0 = pow(b, 2) - (3 * a * c);
		double disc1 = (2 * pow(b, 3)) - (9 * a * b * c) + (27 * pow(a, 2) * d);
		double C = cbrt((disc1 + sqrt(pow(disc1, 2) - (4 * pow(disc0, 3)))) / 2);
		
		//double disc = (pow(disc1, 2) - (4 * pow(disc0, 3)) / (-27 * pow(a, 2)));
		//double C = cbrt((sqrt(pow(disc1, 2) - (4 * pow(disc0, 3))) + disc1) / 2);
		//double u = (-1 + sqrt(-3)) / 2;

		cout << "disc0 = " << disc0 << " disc1 = " << disc1 << " C = " << C << endl << endl;

		eigenval = new matrix(3, 1);

		for (int i = 1; i < 4; i++)
		{
			double root = -(((b + (pow((-1 + sqrt(-3)) / 2, i) * C) + disc0) / pow((-1 + sqrt(-3)) / 2, i)) / 3 * a);
			eigenval->addData(root, i - 1, 0);
		}
	}
	cout << "Eigenvalues:" << endl;
	eigenval->print();

	return eigenval;
}

// Eigenvectros calculation
matrix* eigenvectors(matrix* ref, matrix* eigenvalues)
{
	matrix* eigvector = new matrix(ref->getRows(), eigenvalues->getRows());
	
	if (ref->getRows() == 2)
	{
		for (int i = 0; i < eigenvalues->getRows(); i++)
		{
			eigvector->addData(ref->getData(0, 1) / (eigenvalues->getData(i, 0) - ref->getData(0, 0)), 0, i);
			eigvector->addData(1, 1, i);
		}
	}
	else
	{

	}

	cout << "Eigenvector:" << endl;
	eigvector->print();

	return eigvector;
}

// Best eigenvector selection
matrix* eigenvecBest(matrix* eigenvec, matrix* eigenval)
{
	int best = 0;
	
	for (int i = 1; i < eigenval->getRows(); i++)
	{
		if (eigenval->getData(i, 0) > eigenval->getData(best, 0))
		{
			best = i;
		}
	}

	matrix* bestvec = new matrix(eigenvec->getRows(), 1);

	for (int i = 0; i < eigenvec->getRows(); i++)
	{
		bestvec->addData(eigenvec->getData(i, best), i, 0);
	}

	return bestvec;
}

// PCA calculation
matrix* PCA(matrix* x, matrix* y)
{
	matrix* xy = new matrix(x->getRows(), x->getCols());
	
	for (int i = 0; i < xy->getRows(); i++)
	{
		for (int j = 0; j < x->getCols(); j++)
		{
			xy->addData(x->getData(i, j+1), i, j);
		}
		xy->addData(y->getData(i, 0), i, x->getCols()-1);
	}

	matrix* cov = covariance(xy);
	matrix* eigenval = eigenvalues(cov);
	matrix* eigenvec = eigenvectors(cov, eigenval);
	matrix* evBest = eigenvecBest(eigenvec, eigenval);
	matrix* evBestT = matTrans(evBest);
	matrix* xyT = matTrans(xy);
	matrix* result = matMult(evBestT, xyT);
	result = matTrans(result);

	cout << "Result:" << endl;
	result->print();

	return result;
}

// Alps Water Calculus
void AlpsCalc(void)
{
	cout << "------------------------  Calculating Alps Water  ------------------------" << endl << endl;

	matrix* x = new matrix(17, 2);
	for (int i = 0; i < x->getRows(); i++)
	{
		x->addData(1, i, 0);
		x->addData(alps[i][0], i, 1);
	}

	matrix* y = new matrix(17, 1);
	for (int i = 0; i < y->getRows(); i++)
	{
		y->addData(alps[i][1], i, 0);
	}

	//matrix* lin = linear(x, y);
	//matrix* qua = quad(x, y);
	//matrix* rob = robust(x, y, lin);

	matrix* pca = PCA(x, y);
}

// Books and grades calculus
void BooksCalc(void)
{
	cout << "------------------------  Calculating Books and Grades  ------------------------" << endl << endl;

	matrix* x = new matrix(40, 3);
	for (int i = 0; i < x->getRows(); i++)
	{
		x->addData(1, i, 0);
		x->addData(books[i][0], i, 1);
		x->addData(books[i][1], i, 2);
	}

	matrix* y = new matrix(40, 1);
	for (int i = 0; i < y->getRows(); i++)
	{
		y->addData(books[i][2], i, 0);
	}

	matrix* xy = new matrix(40, 3);
	for (int i = 0; i < xy->getRows(); i++)
	{
		xy->addData(x->getData(i, 1), i, 0);
		xy->addData(x->getData(i, 2), i, 1);
		xy->addData(y->getData(i, 0), i, 2);
	}

	matrix* cov = covariance(xy);
	matrix* eigenval = eigenvalues(cov);

	//matrix* lin = linear(x, y);
	//matrix* qua = quad(x, y);
	//matrix* rob = robust(x, y, lin);
}

// US census calculation
void CensusCalc(void)
{
	cout << "------------------------  Calculating US Census  ------------------------" << endl << endl;

	matrix* x = new matrix(11, 2);
	for (int i = 0; i < x->getRows(); i++)
	{
		x->addData(1, i, 0);
		x->addData(census[i][0], i, 1);
	}

	matrix* y = new matrix(11, 1);
	for (int i = 0; i < y->getRows(); i++)
	{
		y->addData(census[i][1], i, 0);
	}

	//matrix* lin = linear(x, y);
	//matrix* qua = quad(x, y);
	//matrix* rob = robust(x, y, lin);

	matrix* pca = PCA(x, y);
}

// Main loop
int main()
{
	AlpsCalc();

	BooksCalc();

	CensusCalc();

	cin.get();
}