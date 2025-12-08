#include "matrix.h"
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <stdexcept>

Matrix::Matrix(int n) : size(n) {
    if (n <= 0) {
        throw std::invalid_argument("Размер матрицы должен быть положительным");
    }
    data = std::vector<std::vector<double>>(n, std::vector<double>(n, 0.0));
}

Matrix::Matrix(int n, bool random) : Matrix(n) {
    if (random) {
        fillRandom();
    }
}

void Matrix::fillRandom() {
    std::srand(std::time(nullptr));
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            data[i][j] = static_cast<double>(std::rand()) / RAND_MAX * 10.0;
        }
    }
}

void Matrix::fillZero() {
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            data[i][j] = 0.0;
        }
    }
}

void Matrix::print() const {
    std::cout << "Матрица " << size << "x" << size << ":\n";
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            std::cout << std::fixed << std::setprecision(2) 
                     << std::setw(8) << data[i][j] << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";
}

bool Matrix::compare(const Matrix& other, double epsilon) const {
    if (size != other.size) return false;
    
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (std::fabs(data[i][j] - other.data[i][j]) > epsilon) {
                return false;
            }
        }
    }
    return true;
}