
#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <stdlib.h>
#include <windows.h>
#include <xmmintrin.h> 
#include <immintrin.h> 
#include <malloc.h> 
using namespace std;

// 串行高斯消元函数
void serial_algorithm(std::vector<std::vector<float>>& A, std::vector<float>& b) {
    int n = A.size();
    // 消去过程
    for (int k = 0; k < n; ++k) {
        for (int i = k + 1; i < n; ++i) {
            float factor = A[i][k] / A[k][k];
            for (int j = k; j < n; ++j) {
                A[i][j] -= factor * A[k][j];
            }
            b[i] -= factor * b[k];
        }
    }
}

// 串行高斯消元函数 cache优化
void serial_algorithm_cache(std::vector<std::vector<float>>& A, std::vector<float>& b) {
    int n = A.size();
    int blockSize = 64; // 假设以64作为一个合适的块大小

    // 消去过程
    for (int k = 0; k < n; k += blockSize) {
        int K = ((blockSize) < (n - k) ? (blockSize) : (n - k));
        for (int i = k; i < k + K; ++i) {
            for (int j = i + 1; j < k + K; ++j) {
                float factor = A[j][i] / A[i][i];
                for (int l = i + 1; l < k + K; ++l) {
                    A[j][l] -= factor * A[i][l];
                }
                b[j] -= factor * b[i];
            }
        }

        // 更新剩下的行
        for (int i = k + K; i < n; ++i) {
            for (int j = k; j < k + K; ++j) {
                float factor = A[i][j] / A[j][j];
                for (int l = j + 1; l < k + K; ++l) {
                    A[i][l] -= factor * A[j][l];
                }
                b[i] -= factor * b[j];
            }
        }

        // 更新剩下的列
        for (int i = k; i < k + K; ++i) {
            for (int j = k + K; j < n; ++j) {
                float factor = A[j][i] / A[i][i];
                for (int l = i + 1; l < n; ++l) {
                    A[j][l] -= factor * A[i][l];
                }
                b[j] -= factor * b[i];
            }
        }
    }
}

// 使用SSE优化的高斯消元函数 对齐
void SSE_algorithm_align(std::vector<std::vector<float>>& A, std::vector<float>& b) {
    int n = A.size();
    for (int k = 0; k < n; ++k) {
        float pivot = A[k][k];
        if (pivot == 0) {
            cerr << "Pivot element is zero." << endl;
            return;
        }
        __m128 vaik;
        int j = 0;
        for (int i = k + 1; i < n; ++i) {
            float factor = A[i][k] / pivot;
            vaik = _mm_set1_ps(factor);
            for (j = k; (j < n) && ((reinterpret_cast<uintptr_t>(&A[i][j]) % 16) != 0); ++j) {
                A[i][j] -= factor * A[k][j];
            }
            for (; j + 4 <= n - 4; j += 4) {
                __m128 vakj = _mm_load_ps(&A[k][j]);
                __m128 vaij = _mm_load_ps(&A[i][j]);
                vaij = _mm_sub_ps(vaij, _mm_mul_ps(vakj, vaik));
                _mm_store_ps(&A[i][j], vaij);
            }
            for (; j < n; ++j) {
                A[i][j] -= factor * A[k][j];
            }
            b[i] -= factor * b[k];
        }
        A[k][k] = 1.0f; // Pivot is now normalized to 1
        for (j = k + 1; j < n; ++j) A[k][j] /= pivot; // Normalize the pivot row
        b[k] /= pivot;
    }
}
// 使用SSE优化的高斯消元函数 未对齐
void SSE_algorithm_un(std::vector<std::vector<float>>& A, std::vector<float>& b) {
    int n = A.size();
    for (int k = 0; k < n; ++k) {
        float pivot = A[k][k];
        if (pivot == 0) {
            cerr << "Pivot element is zero." << endl;
            return;
        }
        for (int i = k + 1; i < n; ++i) {
            float factor = A[i][k] / pivot;
            __m128 vaik = _mm_set1_ps(factor);
            int j = k;
            for (; j + 4 <= n - 4; j += 4) {
                __m128 vakj = _mm_loadu_ps(&A[k][j]);
                __m128 vaij = _mm_loadu_ps(&A[i][j]);
                vaij = _mm_sub_ps(vaij, _mm_mul_ps(vakj, vaik));
                _mm_storeu_ps(&A[i][j], vaij);
            }
            for (; j < n; ++j) {
                A[i][j] -= factor * A[k][j];
            }
            b[i] -= factor * b[k];
        }
        A[k][k] = 1.0f;
        for (int j = k + 1; j < n; ++j) A[k][j] /= pivot; 
        b[k] /= pivot;
    }
}

void AVX_algorithm_align(std::vector<std::vector<float>>& A, std::vector<float>& b) {
    int n = A.size();
    for (int k = 0; k < n; ++k) {
        float pivot = A[k][k];
        if (pivot == 0) {
            cerr << "Pivot element is zero." << endl;
            return;
        }
        __m256 vt = _mm256_set1_ps(pivot);
        for (int i = k + 1; i < n; ++i) {
            float factor = A[i][k] / pivot;
            __m256 vaik = _mm256_set1_ps(factor);
            int j = k;
            // 为了确保数据对齐到 32 字节边界
            while ((reinterpret_cast<uintptr_t>(&A[i][j]) % 32) != 0 && j < n) {
                A[i][j] -= factor * A[k][j];
                j++;
            }
            // 使用对齐的 AVX 指令
            for (; j + 8 <= n-8; j += 8) {
                __m256 vakj = _mm256_load_ps(&A[k][j]);
                __m256 vaij = _mm256_load_ps(&A[i][j]);
                vaij = _mm256_sub_ps(vaij, _mm256_mul_ps(vakj, vaik));
                _mm256_store_ps(&A[i][j], vaij);
            }
            // 处理剩余的部分
            for (; j < n; ++j) {
                A[i][j] -= factor * A[k][j];
            }
            b[i] -= factor * b[k];
            A[i][k] = 0;
        }
    }
}

void AVX_algorithm_un(std::vector<std::vector<float>>& A, std::vector<float>& b) {
    int n = A.size();
    for (int k = 0; k < n; ++k) {
        float pivot = A[k][k];
        if (pivot == 0) {
            cerr << "Pivot element is zero." << endl;
            return;
        }
        __m256 vt = _mm256_set1_ps(pivot);
        for (int i = k + 1; i < n; ++i) {
            float factor = A[i][k] / pivot;
            __m256 vaik = _mm256_set1_ps(factor);
            int j = k;
            for (; j + 8 <= n-8; j += 8) {
                __m256 vakj = _mm256_loadu_ps(&A[k][j]);
                __m256 vaij = _mm256_loadu_ps(&A[i][j]);
                vaij = _mm256_sub_ps(vaij, _mm256_mul_ps(vakj, vaik));
                _mm256_storeu_ps(&A[i][j], vaij);
            }
            for (; j < n; ++j) {
                A[i][j] -= factor * A[k][j];
            }
            b[i] -= factor * b[k];
            A[i][k] = 0;
        }
    }
}

// 生成测试矩阵的函数
void generateMatrix(std::vector<std::vector<float>>& A, std::vector<float>& b, int N) {
    A.resize(N, std::vector<float>(N, 0));
    b.resize(N);

    srand(time(0)); // 初始化随机数生成器

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            A[i][j] = (i >= j) ? 1.0f + rand() % 100 : 0;
        }
        b[i] = rand() % 100;
    }

    // 使矩阵为上三角形
    for (int k = 0; k < N; k++) {
        for (int i = k + 1; i < N; i++) {
            for (int j = 0; j < N; j++) {
                A[i][j] += A[k][j];
            }
            b[i] += b[k];
        }
    }
}
int main() {

    //for (int N = 2; N <= 2048; N *= 2) {}
        int N = 1024;
        long long head1, tail1, freq;
        QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
        std::vector<std::vector<float>> A;
        std::vector<float> b;
        generateMatrix(A, b, N);

        QueryPerformanceCounter((LARGE_INTEGER*)&head1);
        serial_algorithm(A, b);
        QueryPerformanceCounter((LARGE_INTEGER*)&tail1);
        cout << "问题规模为" << N << endl;
        cout << "串行: " << (tail1 - head1) * 1000.0 / freq << "ms" << endl;

        /*long long head6, tail6;
        QueryPerformanceCounter((LARGE_INTEGER*)&head6);
        serial_algorithm_cache(A, b);
        QueryPerformanceCounter((LARGE_INTEGER*)&tail6);
        cout << "串行 cache: " << (tail6 - head6) * 1000.0 / freq << "ms" << endl;*/

        /*long long head2, tail2;
        QueryPerformanceCounter((LARGE_INTEGER*)&head2);
        SSE_algorithm_align(A, b);
        QueryPerformanceCounter((LARGE_INTEGER*)&tail2);
        cout << "问题规模为" << N << endl;
        cout << "SSE 对齐: " << (tail2 - head2) * 1000.0 / freq << "ms" << endl;*/

        /*long long head3, tail3;
        QueryPerformanceCounter((LARGE_INTEGER*)&head3);
        SSE_algorithm_un(A, b);
        QueryPerformanceCounter((LARGE_INTEGER*)&tail3);
        cout << "SSE 未对齐: " << (tail3 - head3) * 1000.0 / freq << "ms" << endl;
        cout << endl;*/

        /*long long head4, tail4;
        QueryPerformanceCounter((LARGE_INTEGER*)&head4);
        AVX_algorithm_align(A, b);
        QueryPerformanceCounter((LARGE_INTEGER*)&tail4);
        cout << "问题规模为" << N << endl;
        cout << "AVX 对齐: " << (tail4 - head4) * 1000.0 / freq << "ms" << endl;*/


        /*generateMatrix(A, b, N);
        long long head5, tail5;
        QueryPerformanceCounter((LARGE_INTEGER*)&head5);
        AVX_algorithm_un(A, b);
        QueryPerformanceCounter((LARGE_INTEGER*)&tail5);
        cout << "AVX 未对齐: " << (tail5 - head5) * 1000.0 / freq << "ms" << endl;*/

        cout << endl;
    
        
    return 0;
}
//
