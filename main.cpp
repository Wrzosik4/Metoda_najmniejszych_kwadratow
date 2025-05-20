#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>
#include <iomanip>

// Obliczanie potęg x i iloczynów z y
std::vector<std::vector<double>> x_powers(const std::vector<double>& x, const std::vector<double>& y) {
    std::vector<double> x0, x1, x2, x3, x4, x0y, xy, x2y;

    for (int i = 0; i < x.size(); i++) {
        x0.push_back(1.0);
        x1.push_back(x[i]);
        x2.push_back(x[i] * x[i]);
        x3.push_back(x[i] * x[i] * x[i]);
        x4.push_back(x[i] * x[i] * x[i] * x[i]);
    }

    for (int i = 0; i < y.size(); i++) {
        x0y.push_back(x0[i] * y[i]);
        xy.push_back(x1[i] * y[i]);
        x2y.push_back(x2[i] * y[i]);
    }

    return {x0, x1, x2, x3, x4, x0y, xy, x2y};
}

// Budowa macierzy A i wektora b
std::pair<std::vector<std::vector<double>>, std::vector<double>> build_matrix(const std::vector<double>& x, const std::vector<double>& y, int n) {
    std::vector<std::vector<double>> results = x_powers(x, y);

    std::vector<double> sums;
    for (const auto& row : results) {
        sums.push_back(std::accumulate(row.begin(), row.end(), 0.0));
    }

    std::vector<std::vector<double>> A(n+1, std::vector<double>(n+1));
    for (int i = 0; i <= n; i++) {
        for (int j = 0; j <= n; ++j) {
            A[i][j] = sums[i + j];
        }
    }

    std::vector<double> b(n+1);
    for (int i = 0; i <= n; ++i) {
        b[i] = sums[5 + i];
    }

    return {A, b};
}

// Rozwiązanie układu Ax = b metodą Gaussa
std::vector<double> gauss_solve(std::vector<std::vector<double>> A, std::vector<double> b) {
    int n = b.size();

    for (int i = 0; i < n; ++i) {
        A[i].push_back(b[i]);
    }

    for (int i = 0; i < n; ++i) {
        int maxRow = i;
        for (int k = i + 1; k < n; ++k) {
            if (std::abs(A[k][i]) > std::abs(A[maxRow][i])) {
                maxRow = k;
            }
        }
        std::swap(A[i], A[maxRow]);

        double diag = A[i][i];
        for (int j = i; j <= n; ++j) {
            A[i][j] /= diag;
        }

        for (int k = i + 1; k < n; ++k) {
            double factor = A[k][i];
            for (int j = i; j <= n; ++j) {
                A[k][j] -= factor * A[i][j];
            }
        }
    }

    std::vector<double> x(n);
    for (int i = n - 1; i >= 0; --i) {
        x[i] = A[i][n];
        for (int j = i + 1; j < n; ++j) {
            x[i] -= A[i][j] * x[j];
        }
    }

    return x;
}

// Obliczanie W(x) dla danego x
double evaluate_polynomial(const std::vector<double>& coeffs, double x_val) {
    double result = 0.0;
    for (int i = 0; i < coeffs.size(); ++i) {
        result += coeffs[i] * std::pow(x_val, i);
    }
    return result;
}

int main() {
    int n = 2;
    std::vector<double> x = {1, 2, 3, 4};
    std::vector<double> y = {6, 19, 40, 69};

    // Obliczenia potęg i sum
    std::vector<std::vector<double>> results = x_powers(x, y);
    std::vector<std::string> labels = {"x^0", "x^1", "x^2", "x^3", "x^4", "x^0*y", "x^1*y", "x^2*y"};

    std::cout << "Tabela wartości i sum:\n";
    for (int row = 0; row < results.size(); ++row) {
        std::cout << labels[row] << ": ";
        double sum = 0.0;
        for (double val : results[row]) {
            std::cout << std::setw(6) << val << " ";
            sum += val;
        }
        std::cout << " | Suma: " << sum << "\n";
    }

    // Budowa macierzy
    auto [A, b] = build_matrix(x, y, n);

    std::cout << "\nMacierz A:\n";
    for (const auto& row : A) {
        for (double val : row) {
            std::cout << std::setw(10) << val << " ";
        }
        std::cout << "\n";
    }

    std::cout << "\nWektor b:\n";
    for (double val : b) {
        std::cout << std::setw(10) << val << "\n";
    }

    // Rozwiązanie układu
    std::vector<double> coeffs = gauss_solve(A, b);

    std::cout << "\nWspółczynniki wielomianu:\n";
    for (int i = 0; i < coeffs.size(); ++i) {
        std::cout << "a" << i << " = " << coeffs[i] << "\n";
    }

    // Wyświetlenie wzoru W(x)
    std::cout << "\nW(x) = ";
    for (int i = 0; i < coeffs.size(); ++i) {
        if (i > 0) std::cout << " + ";
        std::cout << coeffs[i];
        if (i > 0) std::cout << "*x^" << i;
    }
    std::cout << "\n";

    // Obliczenie W(x) dla podanej wartości x
    double x_val;
    std::cout << "\nPodaj wartość x do obliczenia W(x): ";
    std::cin >> x_val;

    double result = evaluate_polynomial(coeffs, x_val);
    std::cout << "W(" << x_val << ") = " << result << "\n";

    return 0;
}
