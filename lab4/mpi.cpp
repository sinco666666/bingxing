#include <iostream>
#include <mpi.h>
#include <cmath>
#include <chrono>


using namespace std;
using namespace std::chrono;


int N;
const int L = 100;
int LOOP = 1;
float** origin_data;
float** matrix = nullptr;

int NUM_THREADS = 8;

void init_data() {
    origin_data = new float* [N], matrix = new float* [N];
    auto* tmp = new float[N * N];
    for (int i = 0; i < N; i++) {
        origin_data[i] = new float[N], matrix[i] = tmp + i * N;
    }
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            matrix[i][j] = 0, origin_data[i][j] = 0;
        }
    }
    for (int i = 0; i < N; i++) {
        for (int j = i; j < N; j++) {
            origin_data[i][j] = rand() * 1.0 / RAND_MAX * L;
        }
    }
    for (int i = 0; i < N - 1; i++) {
        for (int j = i + 1; j < N; j++) {
            for (int k = 0; k < N; k++) {
                origin_data[j][k] += origin_data[i][k];
            }
        }
    }
}

void init_matrix() {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            matrix[i][j] = origin_data[i][j];
        }
    }
}

void calculate_serial() {
    for (int k = 0; k < N; k++) {
        for (int j = k + 1; j < N; j++) {
            matrix[k][j] = matrix[k][j] / matrix[k][k];
        }
        matrix[k][k] = 1;
        for (int i = k + 1; i < N; i++) {
            for (int j = k + 1; j < N; j++) {
                matrix[i][j] = matrix[i][j] - matrix[i][k] * matrix[k][j];
            }
            matrix[i][k] = 0;
        }
    }
}
// MPI_block �����㷨
double calculate_MPI_block() {
    double start_time, end_time;

    int rank;
    int size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    // ֻ����0�Ž��̣��Ž��г�ʼ������
    if (rank == 0) {
        init_matrix();
    }
    start_time = MPI_Wtime();
    int task_num = ceil(N * 1.0 / size);
    // 0�Ž��̸�������ĳ�ʼ�ַ�����
    if (rank == 0) {
        for (int i = 1; i < size; i++) {
            int start = i * task_num;
            int end = (i + 1) * task_num;
            if (i == size - 1) {
                end = N;
            }
            MPI_Send(&matrix[start][0], (end - start) * N, MPI_FLOAT, i, 0, MPI_COMM_WORLD);
        }
    }
    // ��0�Ž��̸�������Ľ��չ���
    else {
        if (rank != size - 1) {
            MPI_Recv(&matrix[rank * task_num][0], task_num * N, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        else {
            MPI_Recv(&matrix[rank * task_num][0], (N - rank * task_num) * N, MPI_FLOAT, 0, 0, MPI_COMM_WORLD,
                MPI_STATUS_IGNORE);
        }
    }

    // ����Ԫ����
    int start = rank * task_num;
    int end = (rank + 1) * task_num < N ? (rank + 1) * task_num : N;
    for (int k = 0; k < N; k++) {
        // ������������Ǳ����̸�������񣬲�����������㲥
        if (k >= start && k < end) {
            for (int j = k + 1; j < N; j++) {
                matrix[k][j] /= matrix[k][k];
            }
            matrix[k][k] = 1;
            for (int p = 0; p < size; p++) {
                if (p != rank) {
                    MPI_Send(&matrix[k][0], N, MPI_FLOAT, p, 0, MPI_COMM_WORLD);
                }
            }
        }
        // ������̽��ճ����еĽ��
        else {
            MPI_Recv(&matrix[k][0], N, MPI_FLOAT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        // ������Ԫ����
        for (int i = max(k + 1, start); i < end; i++) {
            for (int j = k + 1; j < N; j++) {
                matrix[i][j] = matrix[i][j] - matrix[i][k] * matrix[k][j];
            }
            matrix[i][k] = 0;
        }
    }
    end_time = MPI_Wtime();
    return (end_time - start_time) * 1000;
}
// MPI_cycle �����㷨
double calculate_MPI_cycle() {
    double start_time, end_time;

    int rank;
    int size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    // ֻ����0�Ž��̣��Ž��г�ʼ������
    if (rank == 0) {
        init_matrix();
    }
    start_time = MPI_Wtime();
    int task_num = rank < N % size ? N / size + 1 : N / size;
    // 0�Ž��̸�������ĳ�ʼ�ַ�����
    auto* buff = new float[task_num * N];
    if (rank == 0) {
        for (int p = 1; p < size; p++) {
            for (int i = p; i < N; i += size) {
                for (int j = 0; j < N; j++) {
                    buff[i / size * N + j] = matrix[i][j];
                }
            }
            int count = p < N % size ? N / size + 1 : N / size;
            MPI_Send(buff, count * N, MPI_FLOAT, p, 0, MPI_COMM_WORLD);
        }
    }
    // ��0�Ž��̸�������Ľ��չ���
    else {
        MPI_Recv(&matrix[rank][0], task_num * N, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for (int i = 0; i < task_num; i++) {
            for (int j = 0; j < N; j++) {
                matrix[rank + i * size][j] = matrix[rank + i][j];
            }
        }
    }
    // ����Ԫ����
    for (int k = 0; k < N; k++) {
        // ������������Ǳ����̸�������񣬲�����������㲥
        if (k % size == rank) {
            for (int j = k + 1; j < N; j++) {
                matrix[k][j] /= matrix[k][k];
            }
            matrix[k][k] = 1;
            for (int p = 0; p < size; p++) {
                if (p != rank) {
                    MPI_Send(&matrix[k][0], N, MPI_FLOAT, p, 1, MPI_COMM_WORLD);
                }
            }
        }
        // ������̽��ճ����еĽ��
        else {
            MPI_Recv(&matrix[k][0], N, MPI_FLOAT, k % size, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        // ������Ԫ����
        int begin = N / size * size + rank < N ? N / size * size + rank : N / size * size + rank - size;
        for (int i = begin; i > k; i -= size) {
            for (int j = k + 1; j < N; j++) {
                matrix[i][j] = matrix[i][j] - matrix[i][k] * matrix[k][j];
            }
            matrix[i][k] = 0;
        }
    }
    end_time = MPI_Wtime();
    return (end_time - start_time) * 1000;
}

// MPI pipeline �����㷨
double calculate_MPI_pipeline() {
    double start_time, end_time;

    int rank;
    int size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    // ֻ����0�Ž��̣��Ž��г�ʼ������
    if (rank == 0) {
        init_matrix();
    }
    start_time = MPI_Wtime();
    int task_num = rank < N % size ? N / size + 1 : N / size;
    // 0�Ž��̸�������ĳ�ʼ�ַ�����
    auto* buff = new float[task_num * N];
    if (rank == 0) {
        for (int p = 1; p < size; p++) {
            for (int i = p; i < N; i += size) {
                for (int j = 0; j < N; j++) {
                    buff[i / size * N + j] = matrix[i][j];
                }
            }
            int count = p < N % size ? N / size + 1 : N / size;
            MPI_Send(buff, count * N, MPI_FLOAT, p, 0, MPI_COMM_WORLD);
        }
    }
    // ��0�Ž��̸�������Ľ��չ���
    else {
        MPI_Recv(&matrix[rank][0], task_num * N, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for (int i = 0; i < task_num; i++) {
            for (int j = 0; j < N; j++) {
                matrix[rank + i * size][j] = matrix[rank + i][j];
            }
        }
    }
    // ����Ԫ����
    int pre_proc = (rank + (size - 1)) % size;
    int next_proc = (rank + 1) % size;
    for (int k = 0; k < N; k++) {
        // ������������Ǳ����̸�������񣬲�����������㲥
        if (k % size == rank) {
            for (int j = k + 1; j < N; j++) {
                matrix[k][j] /= matrix[k][k];
            }
            matrix[k][k] = 1;
            MPI_Send(&matrix[k][0], N, MPI_FLOAT, next_proc, 1, MPI_COMM_WORLD);
        }
        // ������̽��ճ����еĽ��
        else {
            MPI_Recv(&matrix[k][0], N, MPI_FLOAT, pre_proc, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            if (next_proc != k % size) {
                MPI_Send(&matrix[k][0], N, MPI_FLOAT, next_proc, 1, MPI_COMM_WORLD);
            }
        }
        // ������Ԫ����
        int begin = N / size * size + rank < N ? N / size * size + rank : N / size * size + rank - size;
        for (int i = begin; i > k; i -= size) {
            for (int j = k + 1; j < N; j++) {
                matrix[i][j] = matrix[i][j] - matrix[i][k] * matrix[k][j];
            }
            matrix[i][k] = 0;
        }
    }
    end_time = MPI_Wtime();
    return (end_time - start_time) * 1000;
}



void print_matrix() {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            printf("%.2f ", matrix[i][j]);
        }
        printf("\n");
    }
}

void test(int n) {
    N = n;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        cout << N << endl;
    }
    double time = 0;
    init_data();

    if (rank == 0) {
        time = 0;
        for (int i = 0; i < LOOP; i++) {
            init_matrix();
            auto start = high_resolution_clock::now();
            calculate_serial();
            auto end = high_resolution_clock::now();
            time += duration_cast<milliseconds>(end - start).count();
        }
        cout << "serial:" << time / LOOP << "ms" << endl;
    }

    time = 0;
    for (int i = 0; i < LOOP; i++) {
        time += calculate_MPI_block();
    }
    if (rank == 0) {
        cout << "MPI_block:" << time / LOOP << "ms" << endl;
    }

    time = 0;
    for (int i = 0; i < LOOP; i++) {
        time += calculate_MPI_cycle();
    }
    if (rank == 0) {
        cout << "MPI_cycle:" << time / LOOP << "ms" << endl;
    }

    time = 0;
    for (int i = 0; i < LOOP; i++) {
        time += calculate_MPI_pipeline();
    }
    if (rank == 0) {
        cout << "MPI_pipeline:" << time / LOOP << "ms" << endl;
    }
}

int main() {
    MPI_Init(nullptr, nullptr);
    // ���Ժ�������֣���֤���в��Բ��ֶ���ִ��
    LOOP = 50;
    for (int i = 256; i < 1024; i += 256)
        test(i);
    LOOP = 1;
    for (int i = 1024; i <= 3840; i += 256)
        test(i);

    test(10); // �����Ĳ��ԣ����ڴ�ӡ���

    MPI_Finalize();
    return 0;
}
