//#include <iostream>
//#include <vector>
//#include <cstdlib>
//#include <ctime>
//#include <chrono>
//#include <arm_neon.h>
//
//// 使用 ARM NEON 优化的高斯消元算法
//void gaussian_elimination(std::vector<std::vector<float>>& A, std::vector<float>& b) {
//    int n = A.size();
//    for (int k = 0; k < n; ++k) {
//        float32x4_t vk = vdupq_n_f32(A[k][k]);
//        for (int i = k + 1; i < n; ++i) {
//            float32x4_t factor = vdupq_n_f32(A[i][k] / A[k][k]);
//            for (int j = k; j + 4 <= n; j += 4) {
//                float32x4_t vrowi = vld1q_f32(&A[i][j]);
//                float32x4_t vrowk = vld1q_f32(&A[k][j]);
//                vrowi = vmlsq_f32(vrowi, factor, vrowk);
//                vst1q_f32(&A[i][j], vrowi);
//            }
//            for (int j = (n / 4) * 4 + k; j < n; j++) {
//                A[i][j] -= A[i][k] / A[k][k] * A[k][j];
//            }
//            b[i] -= A[i][k] / A[k][k] * b[k];
//            A[i][k] = 0;
//        }
//    }
//}
//
//// 串行高斯消元函数
//void serial_algorithm(std::vector<std::vector<float>>& A, std::vector<float>& b) {
//    int n = A.size();
//    for (int k = 0; k < n; ++k) {
//        for (int i = k + 1; i < n; ++i) {
//            float factor = A[i][k] / A[k][k];
//            for (int j = k; j < n; ++j) {
//                A[i][j] -= factor * A[k][j];
//            }
//            b[i] -= factor * b[k];
//        }
//    }
//}
//
//// 生成随机上三角矩阵和向量 b
//void generate_random_matrix(std::vector<std::vector<float>>& A, std::vector<float>& b, int n) {
//    srand(time(nullptr)); // Seed random number generator
//    A.resize(n, std::vector<float>(n, 0.0f));
//    b.resize(n, 0.0f);
//    for (int i = 0; i < n; i++) {
//        for (int j = 0; j < n; j++) {
//            A[i][j] = (i >= j) ? static_cast<float>(rand() % 100 + 1) : 0.0f;
//        }
//        b[i] = static_cast<float>(rand() % 100 + 1);
//    }
//}
//
//int main() {
//    std::vector<int> sizes = { 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048 };
//
//    for (int n : sizes) {
//        std::cout << "Testing matrix size: " << n << std::endl;
//
//        std::vector<std::vector<float>> A;
//        std::vector<float> b;
//
//        // 测试 NEON 优化算法
//        generate_random_matrix(A, b, n);
//        auto start_time_neon = std::chrono::high_resolution_clock::now();
//        gaussian_elimination(A, b);
//        auto end_time_neon = std::chrono::high_resolution_clock::now();
//        std::chrono::duration<double, std::milli> total_time_neon = end_time_neon - start_time_neon;
//        std::cout << "NEON elimination time for size " << n << ": " << total_time_neon.count() << " milliseconds\n";
//
//        // 测试串行算法
//        generate_random_matrix(A, b, n); // 重新生成随机矩阵和向量
//        auto start_time_serial = std::chrono::high_resolution_clock::now();
//        serial_algorithm(A, b);
//        auto end_time_serial = std::chrono::high_resolution_clock::now();
//        std::chrono::duration<double, std::milli> total_time_serial = end_time_serial - start_time_serial;
//        std::cout << "Serial elimination time for size " << n << ": " << total_time_serial.count() << " milliseconds\n";
//    }
//
//    return 0;
//}
