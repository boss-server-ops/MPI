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
//    //    未能整除划分的剩余部分
//    int remain = N % num_proc;
//
//    int begin = rank * block;
//    //    当前进程为最后一个进程时，需处理剩余部分
//    int end = rank != num_proc - 1 ? begin + block : begin + block + remain;
//    for (int k = 0; k < N; k++)
//    {
//        //        判断当前行是否是自己的任务
//        if (k >= begin && k < end)
//        {
//            for (int j = k + 1; j < N; j++)
//                m[k][j] = m[k][j] / m[k][k];
//            m[k][k] = 1.0;
//            //            向之后的进程发送消息
//            for (int p = rank + 1; p < num_proc; p++)
//                MPI_Send(&m[k], N, MPI_FLOAT, p, 2, MPI_COMM_WORLD);
//        }
//        else
//        {
//            int cur_p = k / block;
//            //            当所处行属于当前进程前一进程的任务，需接收消息
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
//        //        在0号进程进行任务划分
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
//        //        处理完0号进程自己的任务后需接收其他进程处理之后的结果
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
//        //        非0号进程先接收任务
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
//        //        处理完后向零号进程返回结果
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
