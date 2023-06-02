//#include<iostream>
//#include<sys/time.h>
//#include<fstream>
//#include<mpi.h>
//using namespace std;
//const int N = 2000;
//const int task = 1;
//float m[N][N];
//void eliminate(float m[][N], int rank, int num_proc)
//{
//    int block = N / num_proc;
//    //    δ���������ֵ�ʣ�ಿ��
//    int remain = N % num_proc;
//
//    int begin = rank * block;
//    //    ��ǰ����Ϊ���һ������ʱ���账��ʣ�ಿ��
//    int end = rank != num_proc - 1 ? begin + block : begin + block + remain;
//    for (int k = 0; k < N; k++)
//    {
//        //        �жϵ�ǰ���Ƿ����Լ�������
//        if (k >= begin && k < end)
//        {
//            for (int j = k + 1; j < N; j++)
//                m[k][j] = m[k][j] / m[k][k];
//            m[k][k] = 1.0;
//            //            ��֮��Ľ��̷�����Ϣ
//            for (int p = rank + 1; p < num_proc; p++)
//                MPI_Send(&m[k], N, MPI_FLOAT, p, 2, MPI_COMM_WORLD);
//        }
//        else
//        {
//            int cur_p = k / block;
//            //            �����������ڵ�ǰ����ǰһ���̵������������Ϣ
//            if (cur_p < rank)
//                MPI_Recv(&m[k], N, MPI_FLOAT, cur_p, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//        }
//        for (int i = begin; i < end && i < N; i++)
//        {
//            if (i >= k + 1)
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
//    timeval start, finish;
//    int num_proc;
//    int rank;
//    MPI_Init(&argc, &argv);
//    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
//    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//
//    int block = N / num_proc;
//    int remain = N % num_proc;
//    if (rank == 0)
//    {
//        //        ��0�Ž��̽������񻮷�
//        ifstream infile("example.txt");
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
//        gettimeofday(&start, NULL);
//        for (int i = 1; i < num_proc; i++)
//        {
//            if (i != num_proc - 1)
//            {
//                for (int j = 0; j < block; j++)
//                    MPI_Send(&m[i * block + j], N, MPI_FLOAT, i, 0, MPI_COMM_WORLD);
//            }
//            else
//            {
//                for (int j = 0; j < block + remain; j++)
//                    MPI_Send(&m[i * block + j], N, MPI_FLOAT, i, 0, MPI_COMM_WORLD);
//            }
//        }
//        eliminate(m, rank, num_proc);
//        //        ������0�Ž����Լ��������������������̴���֮��Ľ��
//        for (int i = 1; i < num_proc; i++)
//        {
//            if (i != num_proc - 1)
//            {
//                for (int j = 0; j < block; j++)
//                    MPI_Recv(&m[i * block + j], N, MPI_FLOAT, i, 1,
//                        MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//            }
//            else
//            {
//                for (int j = 0; j < block + remain; j++)
//                    MPI_Recv(&m[i * block + j], N, MPI_FLOAT, i, 1,
//                        MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//            }
//        }
//            gettimeofday(&finish, NULL);
//    cout  << ((finish.tv_sec - start.tv_sec) * 1000000.0 + finish.tv_usec - start.tv_usec) / 1000.0 << endl;
//    
//
//    }
//    else
//    {
//        //        ��0�Ž����Ƚ�������
//        if (rank != num_proc - 1)
//        {
//            for (int j = 0; j < block; j++)
//                MPI_Recv(&m[rank * block + j], N, MPI_FLOAT, 0, 0,
//                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//        }
//        else
//        {
//            for (int j = 0; j < block + remain; j++)
//                MPI_Recv(&m[rank * block + j], N, MPI_FLOAT, 0, 0,
//                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//        }
//        eliminate(m, rank, num_proc);
//        //        �����������Ž��̷��ؽ��
//        if (rank != num_proc - 1)
//        {
//            for (int j = 0; j < block; j++)
//                MPI_Send(&m[rank * block + j], N, MPI_FLOAT, 0, 1, MPI_COMM_WORLD);
//        }
//        else
//        {
//            for (int j = 0; j < block + remain; j++)
//                MPI_Send(&m[rank * block + j], N, MPI_FLOAT, 0, 1, MPI_COMM_WORLD);
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
