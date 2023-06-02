//#include<iostream>
//#include<Windows.h>
//#include<fstream>
//#include<mpi.h>
//using namespace std;
//const int N = 2000;
//const int task = 1;
//float m[N][N];
//void eliminate(float m[][N], int rank, int num_proc)
//{
//    //    ���н��̽���1�ε����ļ�������
//    int seg = task * num_proc;
//    for (int k = 0; k < N; k++)
//    {
//        //        �жϵ�ǰ���Ƿ����Լ�������
//        if (int((k % seg) / task) == rank)
//        {
//            for (int j = k + 1; j < N; j++)
//                m[k][j] = m[k][j] / m[k][k];
//            m[k][k] = 1.0;
//            //            ��ɼ�������������̷�����Ϣ
//            for (int p = 0; p < num_proc; p++)
//                if (p != rank)
//                    MPI_Send(&m[k], N, MPI_FLOAT, p, 2, MPI_COMM_WORLD);
//        }
//        else
//        {
//            //            �����ǰ�в����Լ������񣬽������Ե�ǰ�д�����̵���Ϣ
//            MPI_Recv(&m[k], N, MPI_FLOAT, int((k % seg) / task), 2,
//                MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//        }
//        for (int i = k + 1; i < N; i++)
//        {
//            if (int((i % seg) / task) == rank)
//            {
//                for (int j = k + 1; j < N; j++)
//                    m[i][j] = m[i][j] - m[i][k] * m[k][j];
//                m[i][k] = 0.0;
//            }
//        }
//    }
//}
//
//void run(int argc, char* argv[])
//{
//
//    long long head, tail, freq; // timers
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
//        // similar to clocks_per_sec
//        QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
//        // start time
//        QueryPerformanceCounter((LARGE_INTEGER*)&head);
//        //        ��0�Ž��̽������񻮷�
//        for (int i = 0; i < N; i++)
//        {
//            int flag = (i % seg) / task;
//            if (flag == rank)
//                continue;
//            else
//                MPI_Send(&m[i], N, MPI_FLOAT, flag, 0, MPI_COMM_WORLD);
//        }
//        eliminate(m, rank, num_proc);
//        //        ������0�Ž����Լ��������������������̴���֮��Ľ��
//        for (int i = 0; i < N; i++)
//        {
//            int flag = (i % seg) / task;
//            if (flag == rank)
//                continue;
//            else
//                MPI_Recv(&m[i], N, MPI_FLOAT, flag, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//        }
//        QueryPerformanceCounter((LARGE_INTEGER*)&tail);
//        cout << "RecycleCol: " << (tail - head) * 1000.0 / freq << "ms" << endl;
//    }
//    else
//    {
//        //        ��0�Ž����Ƚ�������
//        for (int i = task * rank; i < N; i += seg)
//        {
//            for (int j = 0; j < task && i + j < N; j++)
//                MPI_Recv(&m[i + j], N, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//        }
//        eliminate(m, rank, num_proc);
//        //        �����������Ž��̷��ؽ��
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
