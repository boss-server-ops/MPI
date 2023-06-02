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
//    //    所有进程进行1次迭代的计算行数
//    int seg = task * num_proc;
//    for (int k = 0; k < N; k++)
//    {
//        //        判断当前行是否是自己的任务
//        if (int((k % seg) / task) == rank)
//        {
//            for (int j = k + 1; j < N; j++)
//                m[k][j] = m[k][j] / m[k][k];
//            m[k][k] = 1.0;
//            //            完成计算后向其他进程发送消息
//            for (int p = 0; p < num_proc; p++)
//                if (p != rank)
//                    MPI_Send(&m[k], N, MPI_FLOAT, p, 2, MPI_COMM_WORLD);
//        }
//        else
//        {
//            //            如果当前行不是自己的任务，接收来自当前行处理进程的消息
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
//        //        在0号进程进行任务划分
//        for (int i = 0; i < N; i++)
//        {
//            int flag = (i % seg) / task;
//            if (flag == rank)
//                continue;
//            else
//                MPI_Send(&m[i], N, MPI_FLOAT, flag, 0, MPI_COMM_WORLD);
//        }
//        eliminate(m, rank, num_proc);
//        //        处理完0号进程自己的任务后需接收其他进程处理之后的结果
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
//        //        非0号进程先接收任务
//        for (int i = task * rank; i < N; i += seg)
//        {
//            for (int j = 0; j < task && i + j < N; j++)
//                MPI_Recv(&m[i + j], N, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//        }
//        eliminate(m, rank, num_proc);
//        //        处理完后向零号进程返回结果
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
