//#include<iostream>
//#include<Windows.h>
//#include<fstream>
//#include<mpi.h>
//using namespace std;
//const int n = 2000;
//float m[n][n];
//
//
//int main(int argc, char* argv[])
//{
//
//
//    long long head, tail, freq; // timers
//
//    
//
//int myid, numprocs;
//MPI_Init(&argc, &argv);
//MPI_Comm_rank(MPI_COMM_WORLD, &myid);
//MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
////int distributerow=n/(numprocs-1);
//int distributecol = n / (numprocs);
////0号进程首先完成初始化的工作，再将按行划分的每一行传给不同的进程
//if (myid == 0)
//{
//        ifstream infile("F://example.txt");
//       for (int i = 0; i < n; i++)
//       {
//           for (int j = 0; j < n; j++)
//           {
//               char c;
//               infile >> m[i][j];
//
//           }
//       }
//       infile.close();
//       cout << 'z' << m[0][1] << endl;
//
//
//       // similar to CLOCKS_PER_SEC
//       QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
//       // start time
//       QueryPerformanceCounter((LARGE_INTEGER*)&head);
//    //自定义数据类型
//    for (int i = 1; i < numprocs; i++)
//    {
//        int begin = i * distributecol;
//        int end = begin + distributecol;
//        if (i == numprocs - 1)
//            end = n;
//        MPI_Datatype block;
//        MPI_Type_vector(n, (end - begin), n, MPI_FLOAT, &block);
//        MPI_Type_commit(&block);
//        MPI_Send((void*)(m[0] + begin), 1, block, i, 0, MPI_COMM_WORLD);
//    }
//    //printA();
//}
//else//接受消息后并更新对应的矩阵
//{
//    int begin = myid * distributecol;
//    int end = begin + distributecol;
//    if (myid == numprocs - 1)
//        end = n;
//    MPI_Datatype block;
//    MPI_Type_vector(n, (end - begin), n, MPI_FLOAT, &block);
//    MPI_Type_commit(&block);
//    MPI_Recv((void*)(m[0] + begin), 1, block, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//}
//
////开始进行消去
//int begin = (myid)*distributecol;
//int end = begin + distributecol;
//if (myid == numprocs - 1)
//end = n;
//for (int k = 0; k < n; k++)
//{
//    int source = k / distributecol;//注意source
//    if (source >= numprocs)
//        source = numprocs - 1;
//    MPI_Datatype temcol;
//    MPI_Type_vector(n - k, 1, n, MPI_FLOAT, &temcol);
//    MPI_Type_commit(&temcol);
//    MPI_Bcast((void*)(m[k] + k), 1, temcol, source, MPI_COMM_WORLD);//将用于矩阵更新的列传出
//
//    for (int j = (begin >= (k + 1) ? begin : (k + 1)); j < end; j++)
//        m[k][j] = m[k][j] / m[k][k];
//    m[k][k] = 1;
//
//    for (int j = k + 1; j < n; j++)
//    {
//        for (int i = (begin >= (k + 1) ? begin : (k + 1)); i < end; i++)
//        {
//            m[j][i] = m[j][i] - m[j][k] * m[k][i];
//        }
//        m[j][k] = 0;
//    }
//
//}
//
//
//if (myid != 0)
//{
//    //将每个进程更新后的结果传回
//    //MPI_Send((void *)m[begin],count,MPI_FLOAT,0,1,MPI_COMM_WORLD);
//    MPI_Datatype block;
//    MPI_Type_vector(n, (end - begin), n, MPI_FLOAT, &block);
//    MPI_Type_commit(&block);
//    MPI_Send((void*)(m[0] + begin), 1, block, 0, 1, MPI_COMM_WORLD);
//}
//else
//{
//    for (int i = 1; i < numprocs; i++)
//    {
//        int begin = i * distributecol;
//        int end = begin + distributecol;
//        if (i == numprocs - 1)
//            end = n;
//        MPI_Datatype block;
//        MPI_Type_vector(n, (end - begin), n, MPI_FLOAT, &block);
//        MPI_Type_commit(&block);
//        MPI_Recv((void*)(m[0] + begin), 1, block, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//       
//    }
//    QueryPerformanceCounter((LARGE_INTEGER*)&tail);
//    cout << "serialCol: " << (tail - head) * 1000.0 / freq << "ms" << endl;
//}
//
//MPI_Finalize();
//
//}
