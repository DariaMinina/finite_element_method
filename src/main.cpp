#include <math.h>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>

using Matrix = std::vector <std::vector<double>>;
using Vector = std::vector<double>;

Matrix initializeMatrix(const int &rows, const int &cols) {
    Matrix res(rows);
    for (int i = 0; i < rows; ++i) {
        res[i].resize(cols, 0);
    }
    return res;
}

void zero(Matrix &m) {
    for (int i = 0; i < m.size(); ++i) {
        std::fill(m[i].begin(), m[i].end(), 0);
    }
}

Vector initializeVector(int size) {
    return Vector(size, 0);
}

void zero(Vector &v) {
    std::fill(v.begin(), v.end(), 0);
}

// it breaks contents of m and b
int gauss(Matrix &m, Vector &b) {
    
    double *x;
    int N = m.size();
    x = new double[N];
    double d, s;

    for (int k = 0; k < N; k++) // прямой ход
    {
        for (int j = k + 1; j < N; j++)
        {       
            d = m[j][k] / m[k][k]; // формула (1)
            for (int i = k; i < N; i++)
            {
                m[j][i] = m[j][i] - d * m[k][i]; // формула (2)
            }
            b[j] = b[j] - d * b[k]; // формула (3)
        }
    }
    for (int k = N - 1; k >= 0; k--) // обратный ход
    {
        d = 0;
        for (int j = k + 1; j < N; j++)
        {
            s = m[k][j] * x[j]; // формула (4)
            d = d + s; // формула (4)
        }
        x[k] = (b[k] - d) / m[k][k]; // формула (4)
    }

    for (int k = 0; k < N; k++) 
    {
        b[k] = x[k];
    }
    return 0;
}

Vector solveFEM(const Matrix &localMatrix, const Vector &localVector, const int &n) {
    // assemble matrix
    auto matrix = initializeMatrix(n + 1, n + 1);
    for (int i = 0; i < n; ++i) {
        matrix[i][i] += localMatrix[0][0];
        matrix[i][i + 1] += localMatrix[0][1];
        matrix[i + 1][i] += localMatrix[1][0];
        matrix[i + 1][i + 1] += localMatrix[1][1];
    }

    // assemble vector
    auto vector = initializeVector(n + 1);
    for (int i = 0; i < n; ++i) {
        vector[i] += localVector[0];
        vector[i + 1] += localVector[1];
    }

    vector[0] -= 22 * 5;
    vector[n] = 10;
    matrix[n][n] = 1;
    matrix[n][n - 1] = 0;

    if (gauss(matrix, vector) != 0) {
        std::cerr << "Error: degraded matrix\n";
        exit(-1);
    }

    return vector;
}

void saveToFile(const Vector &vector, const std::string &name) {
    std::ofstream fout(name);
    if (!fout.is_open()) {
        std::cerr << "Error occurred during opening file " << name << '\n';
        exit(-1);
    }

    double L = 2.0 / (vector.size() - 1);
    for (auto i = 0; i < vector.size(); ++i) {
        fout << 3 + i * L << ' ' << vector[i] << '\n';
    }

    fout.close();
}

void saveErrorToFile(const std::vector <std::pair<int, double>> &vector, const std::string &name) {
    std::ofstream fout(name);
    if (!fout.is_open()) {
        std::cerr << "Error occurred during opening file " << name << '\n';
        exit(-1);
    }

    for (const auto &value: vector) {
        fout << value.first << ',' << value.second << '\n';
    }

    fout.close();
}

void compress(Matrix &matrix, Vector &vector) {
    for (int i = 1; i < 3; ++i) {
        for (int j = 0; j < 4; ++j) {
            if (fabs(matrix[j][i]) < 1e-10 || i == j) {
                continue;
            }
            double val = matrix[j][i] / matrix[i][i];
            vector[j] -= val * vector[i];
            for (int k = 0; k < 4; ++k) {
                matrix[j][k] -= val * matrix[i][k];
            }
        }
    }
}

double maxD(const Vector &a, const Vector &b) {
    double maxD = -10;
    for (int i = 0; i < a.size(); ++i) {
        auto d = fabs(a[i] - b[i]);
        if (d > maxD) {
            maxD = d;
        }
    }
    return maxD;
}

Vector solveAnalytical(const int &nodeAmount) {
    Vector result(nodeAmount);
    double L = 2.0 / (nodeAmount - 1);

    for (int i = 0; i < nodeAmount; ++i) {
        double x = 3 + L * i; // рассматриваем решение от x = 3
        double c1 = -3806.0*exp(4.0/11.0 + 3.0)/1369.0 + 310.0/37.0;
        double c2 = 3806.0 * exp((21.0 / 22.0) - 6.0) / 1369.0;
        result[i] = c1 + c2 * exp(37.0 / 22.0 * x) + 12.0 / 37.0 * x;
    }

    return result;
}

Vector solveLinear(const int &nodeAmount) {
    double L = 2.0 / (nodeAmount - 1);

    Matrix localMatrix = {
            {-37.0 / 2.0 + 22.0 / L, 37.0 / 2.0 - 22.0 / L},
            {-37.0 / 2.0 - 22.0 / L, 37.0 / 2.0 + 22.0 / L}
    };
    Vector localVector = {
            6.0 * L,
            6.0 * L
    };

    return solveFEM(localMatrix, localVector, nodeAmount - 1);
}

Vector solveCubic(const int &nodeAmount) {
    double L = 2.0 / (nodeAmount - 1);
    Matrix temp_matrix = {
            {
                    814.0 / (10.0 * L) - 37.0 / 2.0,
                    -4158.0 / (40.0 * L) + 57.0 * 37.0 / 80.0,
                    594.0 / (20.0 * L) - 111.0 / 10.0,
                    -286.0 / (40.0 * L) + 37.0 * 7.0 / 80.0
            },
            {
                    -4158.0 / (40.0 * L) - 57.0 * 37.0 / 80.0,
                    1188.0 / (5.0 * L),
                    -6534.0 / (40.0 * L) + 2997.0 / 80.0,
                    594.0 / (20.0 * L) - 111.0 / 10.0,
            },
            {
                    594.0 / (20.0 * L) + 111.0 / 10.0,
                    -6534.0 / (40.0 * L) - 2997.0 / 80.0,
                    1188.0 / (5.0 * L),
                    -4158.0 / (40.0 * L) + 57.0 * 37.0 / 80.0,
            },
            {
                    -286.0 / (40.0 * L) - 37.0 * 7.0 / 80.0,
                    594.0 / (20.0 * L) + 111.0 / 10.0,
                    -4158.0 / (40.0 * L) - 57.0 * 37.0 / 80.0,
                    814.0 / (10.0 * L) + 37.0 / 2.0,
            }
    };
    Vector temp_vector = {
            3.0 * L / 2.0,
            9.0 * L / 2.0,
            9.0 * L / 2.0,
            3.0 * L / 2.0
    };

    compress(temp_matrix, temp_vector);

    Matrix localMatrix = {
            {temp_matrix[0][0], temp_matrix[0][3]},
            {temp_matrix[3][0], temp_matrix[3][3]}
    };
    Vector localVector = {
            temp_vector[0],
            temp_vector[3]
    };
    return solveFEM(localMatrix, localVector, nodeAmount - 1);
}

int main() {
    const auto analyticSolution_20nodes = solveAnalytical(20);
    const auto analyticSolution_40nodes = solveAnalytical(40);
    const auto linearSolution_20nodes = solveLinear(20);
    const auto linearSolution_40nodes = solveLinear(40);
    const auto cubicSolution_20nodes = solveCubic(20);
    const auto cubicSolution_40nodes = solveCubic(40);

    saveToFile(analyticSolution_20nodes, "data/analyticSolution_20nodes.txt");
    saveToFile(analyticSolution_40nodes, "data/analyticSolution_40nodes.txt");
    saveToFile(linearSolution_20nodes, "data/linearSolution_20nodes.txt");
    saveToFile(linearSolution_40nodes, "data/linearSolution_40nodes.txt");
    saveToFile(cubicSolution_20nodes, "data/cubicSolution_20nodes.txt");
    saveToFile(cubicSolution_40nodes, "data/cubicSolution_40nodes.txt");
    

    // Погрешности
    const auto cubicSolution_40nodes_error = maxD(analyticSolution_40nodes, cubicSolution_40nodes);
    const auto cubicSolution_20nodes_error = maxD(analyticSolution_20nodes, cubicSolution_20nodes);
    const auto linearSolution_40nodes_error = maxD(analyticSolution_20nodes, linearSolution_20nodes);
    const auto linearSolution_20nodes_error = maxD(analyticSolution_40nodes, linearSolution_40nodes);

    std::cout << "Погрешность в линейном решении для 20 узлов:   " << linearSolution_20nodes_error << '\n';
    std::cout << "Погрешность в линейном решении для 40 узлов:   " << linearSolution_40nodes_error << '\n';
    std::cout << "Погрешность в кубическом решении для 20 узлов: " << cubicSolution_20nodes_error << '\n';
    std::cout << "Погрешность в кубическом решении для 40 узлов: " << cubicSolution_40nodes_error << '\n';


    // Вычисление количества узлов для линейного решения при котором точность будет равна кубическому решению при 20-ти узлах
    // TODO построить график
    int n = 100;
    std::vector <std::pair<int, double>> errors;
    while (true) {
        auto analyticSolution_Nnodes = solveAnalytical(n);
        auto linearSolution_Nnodes = solveLinear(n);
        auto error = maxD(analyticSolution_Nnodes, linearSolution_Nnodes);
        std::cout << "Для " << n << " узлов при линейном решении получили точность: " << error << '\n';

        errors.push_back(std::make_pair(n, error));
        n += 200;

        if (n > 4000) {
            break;
        }
    }
    /*
    saveErrorToFile(errors, "data/errors_linearSolution_Nnodes.txt");

    // Табличный вывод (TODO добавить в отчет)
    std::cout << "\n\tТабличные результаты для 20 узлов\n";
    double L = 2.0 / (20 - 1);
    for (int i = 0; i < 20; ++i) {
        double x = 3 + L * i;
        std::cout << x << " & " << analyticSolution_20nodes[i] << " & " << linearSolution_20nodes[i] << " & "
                  << cubicSolution_20nodes[i] <<
                  " & " << linearSolution_20nodes[i] - analyticSolution_20nodes[i] << " & "
                  << cubicSolution_20nodes[i] - analyticSolution_20nodes[i] << "\\\\ \\hline\n";
    }

    std::cout << "\n\tТабличные результаты для 40 узлов\n";
    L = 2.0 / (40 - 1);
    for (int i = 0; i < 40; ++i) {
        double x = 3 + L * i;
        std::cout << x << " & " << analyticSolution_40nodes[i] << " & " << linearSolution_40nodes[i] << " & "
                  << cubicSolution_40nodes[i] <<
                  " & " << linearSolution_40nodes[i] - analyticSolution_40nodes[i] << " & "
                  << cubicSolution_40nodes[i] - analyticSolution_40nodes[i] << "\\\\ \\hline \n";
    }
    */
    
    return 0;
}

/*
 *
 *
 *
 */
