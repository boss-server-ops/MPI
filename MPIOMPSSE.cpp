//#include<iostream>
//#include<Windows.h>
//#include<fstream>
//#include<mpi.h>
//#include<emmintrin.h>
//#include<omp.h>
//using namespace std;
//const int N = 1000;
//const int task = 1;
//float m[N][N];
//void eliminate(float m[][N], int rank, int num_proc)
//{
//    __m128 t1, t2, t3;
//    int seg = task * num_proc;
//#pragma omp parallel num_threads(thread_count)
//    for (int k = 0; k < N; k++)
//    {
//        if (int((k % seg) / task) == rank)
//        {
//            float temp1[4] = { m[k][k], m[k][k], m[k][k], m[k][k] };
//            t1 = _mm_loadu_ps(temp1);
//            int j = k + 1;
//#pragma omp for schedule(guided, 20)
//            for (j; j < N - 3; j += 4)
//            {
//                t2 = _mm_loadu_ps(m[k] + j);
//                t3 = _mm_div_ps(t2, t1);
//                _mm_storeu_ps(m[k] + j, t3);
//            }
//#pragma omp for schedule(guided, 20)
//            for (j; j < N; j++)
//            {
//                m[k][j] = m[k][j] / m[k][k];
//            }
//            m[k][k] = 1.0;
//
//            for (int p = 0; p < num_proc; p++)
//                if (p != rank)
//                    MPI_Send(&m[k], N, MPI_FLOAT, p, 2, MPI_COMM_WORLD);
//        }
//        else
//        {
//            MPI_Recv(&m[k], N, MPI_FLOAT, int((k % seg) / task), 2,
//                MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//        }
//        for (int i = k + 1; i < N; i++)
//        {
//            if (int((i % seg) / task) == rank)
//            {
//                float temp2[4] = { m[i][k], m[i][k], m[i][k], m[i][k] };
//                t1 = _mm_loadu_ps(temp2);
//                int j = k + 1;
//#pragma omp for schedule(guided, 20)
//                for (j; j <= N - 3; j += 4)
//                {
//                    t2 = _mm_loadu_ps(m[i] + j);
//                    t3 = _mm_loadu_ps(m[k] + j);
//                    t3 = _mm_mul_ps(t1, t3);
//                    t2 = _mm_sub_ps(t2, t3);
//                    _mm_storeu_ps(m[i] + j, t2);
//                }
//#pragma omp for schedule(guided, 20)
//                for (j; j < N; j++)
//                    m[i][j] = m[i][j] - m[i][k] * m[k][j];
//                m[i][k] = 0;
//            }
//        }
//    }
//}
//
//void run(int argc, char* argv[])
//{
//
//    long long head, tail, freq; // timers
//    timeval t_start;
//    timeval t_end;
//
//    int num_proc;
//    int rank;
//    MPI_Init(&argc, &argv);
//    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
//    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//
//    int seg = task * num_proc;
//    if (rank == 0)
//    {
//        ifstream infile("F://example.txt");
//        for (int i = 0; i < N; i++)
//        {
//            for (int j = 0; j < N; j++)
//            {
//                char c;
//                infile >> m[i][j];
//
//            }
//        }
//        infile.close();
//        cout << 'z' << m[0][1] << endl;
//        // similar to CLOCKS_PER_SEC
//        QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
//        // start time
//        QueryPerformanceCounter((LARGE_INTEGER*)&head);
//        for (int i = 0; i < N; i++)
//        {
//            int flag = (i % seg) / task;
//            if (flag == rank)
//                continue;
//            else
//                MPI_Send(&m[i], N, MPI_FLOAT, flag, 0, MPI_COMM_WORLD);
//        }
//        eliminate(m, rank, num_proc);
//        for (int i = 0; i < N; i++)
//        {
//            int flag = (i % seg) / task;
//            if (flag == rank)
//                continue;
//            else
//                MPI_Recv(&m[i], N, MPI_FLOAT, flag, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//        }
//        QueryPerformanceCounter((LARGE_INTEGER*)&tail);
//        cout << "serialCol: " << (tail - head) * 1000.0 / freq << "ms" << endl;
//    }
//    else
//    {
//        for (int i = task * rank; i < N; i += seg)
//        {
//            for (int j = 0; j < task && i + j < N; j++)
//                MPI_Recv(&m[i + j], N, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//        }
//        eliminate(m, rank, num_proc);
//        for (int i = task * rank; i < N; i += seg)
//        {
//            for (int j = 0; j < task && i + j < N; j++)
//                MPI_Send(&m[i + j], N, MPI_FLOAT, 0, 1, MPI_COMM_WORLD);
//        }
//    }
//    MPI_Finalize();
//}
//
//int main(int argc, char* argv[])
//{
//
//    run(argc, argv);
//
//}
