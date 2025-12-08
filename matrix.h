#ifndef MATRIX_H
#define MATRIX_H

#include <vector>

class Matrix {
private:
    int size;
    std::vector<std::vector<double>> data;
    
public:
    Matrix(int n);
    Matrix(int n, bool random);
    
    int getSize() const { return size; }
    double get(int i, int j) const { return data[i][j]; }
    void set(int i, int j, double value) { data[i][j] = value; }
    
    std::vector<double>& operator[](int i) { return data[i]; }
    const std::vector<double>& operator[](int i) const { return data[i]; }
    
    void fillRandom();
    void fillZero();
    void print() const;
    bool compare(const Matrix& other, double epsilon = 1e-6) const;
};

#endif