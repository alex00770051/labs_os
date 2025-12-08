#include "matrix.h"
#include <iostream>
#include <iomanip>
#include <sys/mman.h>
#include <sys/wait.h>
#include <unistd.h>
#include <chrono>
#include <vector>
#include <fstream>

using namespace std;
using namespace std::chrono;

const int NUM_MEASUREMENTS = 3;

// Последовательное умножение матриц
double sequentialMultiply(const Matrix& A, const Matrix& B, Matrix& C) {
    auto start = high_resolution_clock::now();
    
    int n = A.getSize();
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double sum = 0.0;
            for (int k = 0; k < n; k++) {
                sum += A[i][k] * B[k][j];
            }
            C[i][j] = sum;
        }
    }
    
    auto end = high_resolution_clock::now();
    return duration_cast<duration<double>>(end - start).count();
}

// Копирование матрицы в разделяемую память
void copyMatrixToShared(const Matrix& matrix, double* shared, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            shared[i * n + j] = matrix[i][j];
        }
    }
}

// Копирование матрицы из разделяемой памяти
void copyMatrixFromShared(Matrix& matrix, double* shared, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            matrix[i][j] = shared[i * n + j];
        }
    }
}

// Работа дочернего процесса
void childProcessWork(int pid, int num_procs, int n, int start_row, int end_row,
                     double* A, double* B, double* C) {
    for (int i = start_row; i < end_row; i++) {
        for (int j = 0; j < n; j++) {
            double sum = 0.0;
            for (int k = 0; k < n; k++) {
                sum += A[i * n + k] * B[k * n + j];
            }
            C[i * n + j] = sum;
        }
    }
    exit(0);
}

// Умножение с использованием процессов
double multiplyWithProcesses(const Matrix& A, const Matrix& B, Matrix& C, int num_procs) {
    int n = A.getSize();
    size_t size = n * n * sizeof(double);
    
    double* shared_A = (double*)mmap(NULL, size, PROT_READ | PROT_WRITE,
                                    MAP_SHARED | MAP_ANONYMOUS, -1, 0);
    double* shared_B = (double*)mmap(NULL, size, PROT_READ | PROT_WRITE,
                                    MAP_SHARED | MAP_ANONYMOUS, -1, 0);
    double* shared_C = (double*)mmap(NULL, size, PROT_READ | PROT_WRITE,
                                    MAP_SHARED | MAP_ANONYMOUS, -1, 0);
    
    copyMatrixToShared(A, shared_A, n);
    copyMatrixToShared(B, shared_B, n);
    
    auto start = high_resolution_clock::now();
    
    vector<pid_t> pids(num_procs);
    int rows_per = n / num_procs;
    int extra = n % num_procs;
    int current = 0;
    
    for (int i = 0; i < num_procs; i++) {
        int start_row = current;
        int end_row = start_row + rows_per;
        if (i < extra) end_row++;
        current = end_row;
        
        pid_t pid = fork();
        if (pid < 0) {
            cerr << "Ошибка создания процесса\n";
            return -1;
        }
        
        if (pid == 0) {
            childProcessWork(i, num_procs, n, start_row, end_row,
                           shared_A, shared_B, shared_C);
        } else {
            pids[i] = pid;
        }
    }
    
    for (int i = 0; i < num_procs; i++) {
        waitpid(pids[i], nullptr, 0);
    }
    
    auto end = high_resolution_clock::now();
    double time = duration_cast<duration<double>>(end - start).count();
    
    copyMatrixFromShared(C, shared_C, n);
    
    munmap(shared_A, size);
    munmap(shared_B, size);
    munmap(shared_C, size);
    
    return time;
}

// Главная функция
int main(int argc, char* argv[]) {
    if (argc != 3) {
        cout << "Использование: " << argv[0] << " <размер_матрицы> <количество_процессов>\n";
        cout << "Пример: " << argv[0] << " 1000 4\n";
        return 1;
    }
    
    int n = stoi(argv[1]);
    int num_procs = stoi(argv[2]);
    
    if (n < 500) {
        cout << "Размер матрицы должен быть не менее 500\n";
        return 1;
    }
    
    if (num_procs < 1) {
        cout << "Количество процессов должно быть не менее 1\n";
        return 1;
    }
    
    cout << "Умножение матриц " << n << "x" << n << " с " << num_procs << " процессами\n";
    cout << "Количество замеров: " << NUM_MEASUREMENTS << "\n\n";
    
    // Создаем матрицы
    Matrix A(n, true);
    Matrix B(n, true);
    Matrix C_seq(n, false);
    Matrix C_par(n, false);
    
    // Последовательное умножение
    cout << "Последовательное умножение (1 процесс):\n";
    double seq_time = sequentialMultiply(A, B, C_seq);
    cout << "Время: " << fixed << setprecision(2) << seq_time << " с\n\n";
    
    // Параллельное умножение с заданным количеством процессов
    cout << "Параллельное умножение с " << num_procs << " процессами:\n";
    
    vector<double> measurements;
    double total_time = 0;
    
    // Делаем 3 замера
    for (int meas = 0; meas < NUM_MEASUREMENTS; meas++) {
        Matrix temp(n, false);
        double time = multiplyWithProcesses(A, B, temp, num_procs);
        measurements.push_back(time);
        total_time += time;
        
        cout << "  Замер " << (meas + 1) << ": " 
             << fixed << setprecision(2) << time << " мс\n";
        
        // Сохраняем результат первого замера для проверки
        if (meas == 0) {
            for (int row = 0; row < n; row++) {
                for (int col = 0; col < n; col++) {
                    C_par[row][col] = temp[row][col];
                }
            }
        }
    }
    
    // Среднее время
    double avg_time = total_time / NUM_MEASUREMENTS;
    cout << "\n  Среднее время: " << fixed << setprecision(2) 
         << avg_time << " с\n";
    
    // Ускорение и эффективность (только для файла)
    double speedup = seq_time / avg_time;
    double efficiency = (speedup / num_procs) * 100.0;
    
    // Проверка корректности
    cout << "\n";
    if (C_seq.compare(C_par)) {
        cout << "Результаты корректны\n";
    } else {
        cout << "Обнаружены расхождения в результатах!\n";
    }
    
    // Сохранение результатов в файл
    ofstream file("results.csv");
    file << "Процессы,Время(с),Ускорение,Эффективность(%)\n";
    file << num_procs << "," 
         << fixed << setprecision(2) << avg_time << ","
         << fixed << setprecision(3) << speedup << ","
         << fixed << setprecision(1) << efficiency << "\n";
    file.close();
    
    cout << "\nРезультаты сохранены в файл results.csv\n";
    
    return 0;
}