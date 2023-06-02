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
//    int seg = task * num_proc;
//    //    ���㵱ǰ���̵�ǰһ���̼���һ����
//    int pre_proc = (rank + (num_proc - 1)) % num_proc;
//    int next_proc = (rank + 1) % num_proc;
//    for (int k = 0; k < N; k++)
//    {
//        //        �жϵ�ǰ���Ƿ����Լ�������
//        if (int((k % seg) / task) == rank)
//        {
//            for (int j = k + 1; j < N; j++)
//                m[k][j] = m[k][j] / m[k][k];
//            m[k][k] = 1.0;
//            //            �������Լ������������һ���̷�����Ϣ
//            MPI_Send(&m[k], N, MPI_FLOAT, next_proc, 2, MPI_COMM_WORLD);
//        }
//        else
//        {
//            //            �����ǰ�в��ǵ�ǰ���̵����������ǰһ���̵���Ϣ
//            MPI_Recv(&m[k], N, MPI_FLOAT, pre_proc, 2,
//                MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//            //            �����ǰ�в�����һ���̵������轫��Ϣ���д���
//            if (int((k % seg) / task) != next_proc)
//                MPI_Send(&m[k], N, MPI_FLOAT, next_proc, 2, MPI_COMM_WORLD);
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
//    long long head, tail, freq; // timers
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
//
//        //        ��0�Ž��̽������񻮷�
//        ifstream infile("F://example.txt");
//        for (int i = 0; i < N; i++)
//        {
//        for (int j = 0; j < N; j++)
//        {
//        char c;
//        infile >> m[i][j];
//
//        }
//        }
//        infile.close();
//        cout << 'z' << m[0][1] << endl;
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
//        cout << "serialCol: " << (tail - head) * 1000.0 / freq << "ms" << endl;
//
//        
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
//    run(argc,argv);
//   
//}
