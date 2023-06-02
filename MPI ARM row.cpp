//#include<iostream>
//#include<sys/time.h>
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
//    timeval start, finish;
//
//
//
//    int myid, numprocs;
//    MPI_Init(&argc, &argv);
//    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
//    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
//    int distributerow = n / (numprocs - 1);
//
//
//    if (myid == 0)
//    {
//        ifstream infile("example.txt");
//        for (int i = 0; i < n; i++)
//        {
//            for (int j = 0; j < n; j++)
//            {
//                char c;
//                infile >> m[i][j];
//
//            }
//        }
//        infile.close();
//        cout <<'z'<< m[0][1] << endl;
//       
//
//        gettimeofday(&start, NULL);
//        for (int i = 1; i < numprocs; i++)
//        {
//            int begin = (i - 1) * distributerow;
//            int end = begin + distributerow;
//            if (i == numprocs - 1)
//                end = n;
//            int count = (end - begin) * n;
//
//            MPI_Send((void*)m[begin], count, MPI_FLOAT, i, 0, MPI_COMM_WORLD);
//        }
//
//    }
//    else
//    {
//        
//
//        int begin = (myid - 1) * distributerow;
//        int end = begin + distributerow;
//        if (myid == numprocs - 1)
//            end = n;
//        int count = (end - begin) * n;
//        MPI_Recv((void*)m[begin], count, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//    }
//    int begin = (myid - 1) * distributerow;
//    int end = begin + distributerow;
//    if (myid == numprocs - 1)
//        end = n;
//    int count = (end - begin) * n;
//    for (int k = 0; k < n; k++)
//    {
//        if (myid == 0)
//        {
//            if (k != 0)
//            {
//                int source = (k / distributerow + 1) < (numprocs - 1) ? (k / distributerow + 1) : (numprocs - 1);
//                MPI_Recv((void*)(m[k] + k), n - k, MPI_FLOAT, source, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//            }
//            for (int j = k + 1; j < n; j++)
//                m[k][j] = m[k][j] / m[k][k];
//            
//            m[k][k] = 1;
//        }
//        
//        MPI_Bcast((void*)(m[k] + k), n - k, MPI_FLOAT, 0, MPI_COMM_WORLD);
//
//      
//        if (myid != 0)
//        {
//            for (int j = (begin > (k + 1) ? begin : k + 1); j < end; j++) 
//            {
//                for (int i = k + 1; i < n; i++)
//                {
//                    m[j][i] = m[j][i] - m[j][k] * m[k][i];
//                }
//                m[j][k] = 0;
//            }
//            if ((k + 1 < n) && (k + 1) >= begin && (k + 1) < end)
//            {
//                MPI_Send((void*)(m[k + 1] + k + 1), n - (k + 1), MPI_FLOAT, 0, 2, MPI_COMM_WORLD);
//            }
//        }
//
//    }
//    if (myid != 0)
//    {
//        MPI_Send((void*)m[begin], count, MPI_FLOAT, 0, 1, MPI_COMM_WORLD);
//    }
//    else
//    {
//        for (int i = 1; i < numprocs; i++)
//        {
//            int begin = (i - 1) * distributerow;
//            int end = begin + distributerow;
//            if (i == numprocs - 1)
//                end = n;
//            int count = (end - begin) * n; 
//            MPI_Recv((void*)m[begin], count, MPI_FLOAT, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//        }
//     gettimeofday(&finish, NULL);
//    cout  << ((finish.tv_sec - start.tv_sec) * 1000000.0 + finish.tv_usec - start.tv_usec) / 1000.0 << endl;
//    }
//    MPI_Finalize();
//
//
//
//
//
//
//
//}
