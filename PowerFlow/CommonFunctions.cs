using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace PowerFlow
{
    public static class CommonFunctions
    {
        public static void FindMaxAbs(List<double> data, out double max)
        {
            //  find the maximum absolute value from list data and return the maximum value.

            try
            {
                max = Math.Abs(data[0]);
                foreach (double value in data)
                {
                    if (Math.Abs(value) > max)
                    {
                        max = Math.Abs(value);
                    }
                }
            }
            catch (Exception ex)
            {
                max = 0.0;
                MessageBox.Show("Error in Finding Maximum: " + ex.Message);
            }
        }

        public static void Ybus(List<List<string>> bus, List<List<string>> branch, out double[,] Gbus, out double[,] Bbus)
        {
            // return Gbus, Bbus (i.e. Ybus) matrix from bus data and branch data.

            #region Useful infomation from bus and branch

            List<int> BusNO = bus[0].ToList().Select(int.Parse).ToList();
            List<int> BusType = bus[6].ToList().Select(int.Parse).ToList();
            List<double> VolMag = bus[7].ToList().Select(double.Parse).ToList();
            List<double> VolAng = bus[8].ToList().Select(double.Parse).ToList();
            List<double> LoadMW = bus[9].ToList().Select(double.Parse).ToList();
            List<double> LoadMVAR = bus[10].ToList().Select(double.Parse).ToList();
            List<double> GenMW = bus[11].ToList().Select(double.Parse).ToList();
            List<double> GenMVAR = bus[12].ToList().Select(double.Parse).ToList();
            List<double> busG = bus[17].ToList().Select(double.Parse).ToList();
            List<double> busB = bus[18].ToList().Select(double.Parse).ToList();

            List<int> fromBus = branch[0].ToList().Select(int.Parse).ToList();
            List<int> toBus = branch[1].ToList().Select(int.Parse).ToList();
            List<double> LineR = branch[6].ToList().Select(double.Parse).ToList();
            List<double> LineX = branch[7].ToList().Select(double.Parse).ToList();
            List<double> LineB = branch[8].ToList().Select(double.Parse).ToList();

            #endregion

            // initialize.
            Gbus = new double[BusNO.Count, BusNO.Count];
            Bbus = new double[BusNO.Count, BusNO.Count];

            //  generate Gbus and Bbus.
            try
            {
                for (int i = 0; i < fromBus.Count(); ++i)
                {
                    int from = fromBus[i] - 1;
                    int to = toBus[i] - 1;
                    Gbus[from, from] = Gbus[from, from] + LineR[i] / (LineR[i] * LineR[i] + LineX[i] * LineX[i]);
                    Bbus[from, from] = Bbus[from, from] - LineX[i] / (LineR[i] * LineR[i] + LineX[i] * LineX[i]) + 0.5 * LineB[i];

                    Gbus[to, to] = Gbus[to, to] + LineR[i] / (LineR[i] * LineR[i] + LineX[i] * LineX[i]);
                    Bbus[to, to] = Bbus[to, to] - LineX[i] / (LineR[i] * LineR[i] + LineX[i] * LineX[i]) + 0.5 * LineB[i];

                    Gbus[from, to] = Gbus[from, to] - LineR[i] / (LineR[i] * LineR[i] + LineX[i] * LineX[i]);
                    Bbus[from, to] = Bbus[from, to] + LineX[i] / (LineR[i] * LineR[i] + LineX[i] * LineX[i]);

                    Gbus[to, from] = Gbus[to, from] - LineR[i] / (LineR[i] * LineR[i] + LineX[i] * LineX[i]);
                    Bbus[to, from] = Bbus[to, from] + LineX[i] / (LineR[i] * LineR[i] + LineX[i] * LineX[i]);
                }                

                for(int j = 0; j < BusNO.Count; ++j)
                {
                    if(busG[j] != 0)
                    {
                        Gbus[j, j] = Gbus[j, j] + busG[j];
                    }
                    if(busB[j] != 0)
                    {
                        Bbus[j, j] = Bbus[j, j] + busB[j];
                    }
                }
            }
            catch (Exception ex)
            {
                MessageBox.Show("Error in forming Ybus: " + ex.Message);
            }
        }

        public static void Jacobian(List<int> BusNO, List<int> BusType, List<double> VolAng, List<double> VolMag, List<double> Pi, List<double> Qi, double[,] Gbus, double[,] Bbus, out double[,] Jacobian, out List<double> dPQ)
        {
            // return the Jacobian matrix from Gbus and Bbus matrix.

            // find the bus type index.
            List<int> PQ, PV, PVPQ, slack;
            FindBusType(BusNO, BusType, out PQ, out PV, out PVPQ, out slack);

            // initialize.
            Jacobian = new double[PVPQ.Count + PQ.Count, PVPQ.Count + PQ.Count];
            List<double> dP = new List<double>();
            List<double> dQ = new List<double>();
            dPQ = new List<double>();
            double[,] J11 = new double[PVPQ.Count, PVPQ.Count];
            double[,] J12 = new double[PVPQ.Count, PQ.Count];
            double[,] J21 = new double[PQ.Count, PVPQ.Count];
            double[,] J22 = new double[PQ.Count, PQ.Count];

            // deltaP & Jacobian J11 i/i
            foreach (int i in PVPQ)
            {
                double sigema = 0.0, sigemaJ = 0.0;
                for(int j = 0; j < BusNO.Count; ++j)
                {
                    sigema = sigema + VolMag[j]*(Gbus[i,j]*(Math.Cos(VolAng[i]-VolAng[j])) + Bbus[i,j]*(Math.Sin(VolAng[i]-VolAng[j])));
                    sigemaJ = sigemaJ + VolMag[j]*(Gbus[i,j]*((-1)*Math.Sin(VolAng[i]-VolAng[j])) + Bbus[i,j]*(Math.Cos(VolAng[i]-VolAng[j])));
                }
                dP.Add(Pi[i] - VolMag[i]*sigema);
                dPQ.Add(Pi[i] - VolMag[i] * sigema);
                J11[i - 1, i - 1] = VolMag[i] * sigemaJ - VolMag[i] * VolMag[i] * Bbus[i, i];
                Jacobian[i - 1, i - 1] = J11[i - 1, i - 1];
            }

            // Jacobian J11 i/j
            foreach (int i in PVPQ)
            {
                foreach (int j in PVPQ)
                {
                    if (i != j)
                    {
                        int indexi = PVPQ.FindIndex(x => x==i);
                        int indexj = PVPQ.FindIndex(x => x==j);
                        J11[indexi, indexj] = VolMag[i]*VolMag[j]*(Gbus[i,j]*(Math.Sin(VolAng[i]-VolAng[j]))-Bbus[i,j]*(Math.Cos(VolAng[i]-VolAng[j])));
                        Jacobian[indexi, indexj] = J11[indexi, indexj];
                    }
                }
            }

            // Jacobian J12 i/i  &  J21 i/i
            foreach (int i in PQ)
            {
                double sigemaJ = 0.0;
                for(int j = 0; j < BusNO.Count; ++j)
                {
                    sigemaJ = sigemaJ + VolMag[j]*(Gbus[i,j]*(Math.Cos(VolAng[i]-VolAng[j])) + Bbus[i,j]*(Math.Sin(VolAng[i]-VolAng[j])));
                }
                int indexi = PVPQ.FindIndex(x => x==i);
                int indexj = PQ.FindIndex(x => x==i);
                J12[indexi, indexj] = sigemaJ + VolMag[i] * Gbus[i, i];
                J21[indexj, indexi] = VolMag[i] * sigemaJ - (VolMag[i] * VolMag[i]) * Gbus[i, i];
                Jacobian[indexi, indexj + PVPQ.Count] = J12[indexi, indexj];
                Jacobian[indexj + PVPQ.Count, indexi] = J21[indexj, indexi];
            }

            // Jacobian J12 i/j
            foreach (int i in PVPQ)
            {
                foreach (int j in PQ)
                {
                    if (i != j)
                    {
                        int indexi = PVPQ.FindIndex(x => x == i);
                        int indexj = PQ.FindIndex(x => x == j);
                        J12[indexi, indexj] = VolMag[i] * (Gbus[i, j] * (Math.Cos(VolAng[i] - VolAng[j])) + Bbus[i, j] * (Math.Sin(VolAng[i] - VolAng[j])));
                        Jacobian[indexi, indexj + PVPQ.Count] = J12[indexi, indexj];
                    }
                }
            }

            // Jacobian J21 i/j
            foreach (int i in PQ)
            {
                foreach (int j in PVPQ)
                {
                    if (i != j)
                    {
                        int indexi = PQ.FindIndex(x => x == i);
                        int indexj = PVPQ.FindIndex(x => x == j);
                        J21[indexi, indexj] = VolMag[i] * VolMag[j] * (-Gbus[i, j] * (Math.Cos(VolAng[i] - VolAng[j])) - Bbus[i, j] * (Math.Sin(VolAng[i] - VolAng[j])));
                        Jacobian[indexi + PVPQ.Count, indexj] = J21[indexi, indexj];
                    }
                }
            }

            // deltaQ & Jacobian J22 i/i
            foreach (int i in PQ)
            {
                double sigema = 0.0;
                for (int j = 0; j < BusNO.Count; ++j)
                {
                    sigema = sigema + VolMag[j] * (Gbus[i, j] * (Math.Sin(VolAng[i] - VolAng[j])) - Bbus[i, j] * (Math.Cos(VolAng[i] - VolAng[j])));
                }
                dQ.Add(Qi[i] - VolMag[i] * sigema); 
                dPQ.Add(Qi[i] - VolMag[i] * sigema);
                int indexi = PQ.FindIndex(x => x == i);
                J22[indexi, indexi] = sigema - VolMag[i] * Bbus[i, i];
                Jacobian[indexi + PVPQ.Count, indexi + PVPQ.Count] = J22[indexi, indexi];
            }

            // Jacobian J22 i/j
            foreach (int i in PQ)
            {
                foreach (int j in PQ)
                {
                    if (i != j)
                    {
                        int indexi = PQ.FindIndex(x => x == i);
                        int indexj = PQ.FindIndex(x => x == j);
                        J22[indexi, indexj] = VolMag[i] * (Gbus[i, j] * (Math.Sin(VolAng[i] - VolAng[j])) - Bbus[i, j] * (Math.Cos(VolAng[i] - VolAng[j])));
                        Jacobian[indexi + PVPQ.Count, indexj + PVPQ.Count] = J22[indexi, indexj];
                    }
                }
            }
        }

        public static void FindBusType(List<int> BusNO, List<int> BusType, out List<int> PQ, out List<int> PV, out List<int> PVPQ, out List<int> slack)
        {
            // find the index of PQ, PV and slack buses.

            PQ = new List<int>();
            PV = new List<int>();
            PVPQ = new List<int>();
            slack = new List<int>();

            for (int i = 0; i < BusNO.Count; ++i)
            {
                int type = BusType[i];
                switch (type)
                {
                    case 0:
                        PQ.Add(i);
                        PVPQ.Add(i);
                        break;
                    case 1:
                        PQ.Add(i);
                        PVPQ.Add(i);
                        break;
                    case 2:
                        PV.Add(i);
                        PVPQ.Add(i);
                        break;
                    case 3:
                        slack.Add(i);
                        break;
                    default:
                        MessageBox.Show("Error in finding bus type.");
                        break;
                }
            }
        }

        public static void LUFact(double[,] A, List<double> b, out List<double> x)
        {
            // use LU factorization to get the inverse of a matrix.

            int row = A.GetLength(0);

            //  get Q matrix using Crout's algorithm
            double[,] Q;
            crout(A, out Q);

            //  forward substitution to get y
            List<double> y = new List<double>();
            for (int i = 0; i < row; ++i)
            {
                y.Add(0.0);
            }
            for (int k = 0; k < row; ++k)
            {
                double sigmay = 0.0;
                for (int j = 0; j < k; ++j)
                {
                    sigmay = sigmay + Q[k, j] * y[j];
                }
                y[k] = (1 / Q[k, k]) * (b[k] - sigmay);
            }

            //  backward substitution to get x            
            x = new List<double>();
            for (int i = 0; i < row; ++i)
            {
                x.Add(0.0);
            }
            for (int k = row-1; k >= 0; --k)
            {
                double sigmax = 0.0;
                for (int j = (k+1); j < row; ++j)
                {
                    sigmax = sigmax + Q[k,j]*x[j];
                }
                x[k] = y[k] - sigmax;
            }
        }

        public static void crout(double[,] A, out double[,] Q)
        {
            //  use Crout's algorithm to compute the element of Q.

            int row = A.GetLength(0);
            int col = A.GetLength(1);
            Q = new double[row, col];

            for (int j = 0; j < col; ++j)
            {
                // the jth column
                for (int k = j; k < row; ++k)
                {
                    double sigmacol = 0;
                    for (int i = 0; i < j; ++i)
                    {
                        sigmacol = sigmacol + Q[k,i]*Q[i,j];
                    }
                    Q[k,j] = A[k,j] - sigmacol;
                }
                // the jth row
                if (Q[j,j] != 0)
                {
                    for (int k = (j+1); k < col; ++k)
                    {
                        double sigmarow = 0;
                        for (int i = 0; i < j; ++i)
                        {
                            sigmarow = sigmarow + Q[j,i]*Q[i,k];
                        }
                        Q[j,k] = (1/Q[j,j])*(A[j,k]-sigmarow);
                    }
                }
            }
        }            
            
    }
}