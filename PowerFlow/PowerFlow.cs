//////////////////////////////////////////////////////////
// Copyright (c) 2016 All Rights Reserved
// Author:          Yang Zheng 
// Last Update:     10/18/2016
// Summary:         Power flow using Newton Raphson method
// Contact:         yang.zheng@wsu.edu
//////////////////////////////////////////////////////////

using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using System.IO;

namespace PowerFlow
{
    public partial class PowerFlow : Form
    {
        public List<List<string>> bus = new List<List<string>>();
        public List<List<string>> branch = new List<List<string>>();
        public static double stoppingError = 0.01;

        public PowerFlow()
        {
            InitializeComponent();
        }

        private void readBus_Click(object sender, EventArgs e)
        {
            //  read bus file and generate a string list of list "bus".

            OpenFileDialog dlg = new OpenFileDialog();
            dlg.Filter = "CSV Files (*.csv)|*.csv";
            if (dlg.ShowDialog() == DialogResult.OK)
            {
                try
                {
                    string[] lines = File.ReadAllLines(dlg.FileName);
                    var sep = new[] { "," };
                    for (int i = 0; i < lines.Count(); i++)
                    {
                        string[] cols = lines[i].Split(sep, StringSplitOptions.RemoveEmptyEntries);
                        bus.Add(new List<string>());
                        for (int j = 0; j < cols.Count(); j++)
                        {
                            bus[i].Add(cols[j]);
                        }
                    }
                    MessageBox.Show("Bus data imported.");
                }
                catch (Exception ex)
                {
                    MessageBox.Show("Could not read bus data. Original Error: " + ex.Message);
                }
            }
        }
        
        private void readBranch_Click(object sender, EventArgs e)
        {
            //  read branch file and generate a string list of list "branch".

            OpenFileDialog dlg = new OpenFileDialog();
            dlg.Filter = "CSV Files (*.csv)|*.csv";
            if (dlg.ShowDialog() == DialogResult.OK)
            {
                try
                {
                    string[] lines = File.ReadAllLines(dlg.FileName);
                    var sep = new[] { "," };
                    for (int i = 0; i < lines.Count(); i++)
                    {
                        string[] cols = lines[i].Split(sep, StringSplitOptions.RemoveEmptyEntries);
                        branch.Add(new List<string>());
                        for (int j = 0; j < cols.Count(); j++)
                        {
                            branch[i].Add(cols[j]);
                        }
                    }
                    MessageBox.Show("Branch data imported.");
                }
                catch (Exception ex)
                {
                    MessageBox.Show("Could not read branch data. Original Error: " + ex.Message);
                }
            }
        }

        private void NewtonRaphson_Click(object sender, EventArgs e)
        {
            //  start Newton-Raphson power flow calculation when user clicks the button.

            try
            {
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

                //  generate the Gbus and Bbus matrix.
                double[,] Gbus, Bbus;
                CommonFunctions.Ybus(bus, branch, out Gbus, out Bbus);

                //  find the PQ, PV, PVPQ, slack bus index.
                List<int> PQ, PV, PVPQ, slack;
                CommonFunctions.FindBusType(BusNO, BusType, out PQ, out PV, out PVPQ, out slack);

                //  initialize Pi, Qi, VolMag, VolAng according to bus type.
                List<double> Pi = new List<double>();
                List<double> Qi = new List<double>();
                for (int i = 0; i < BusNO.Count; ++i)
                {
                    if (BusType[i] == 3)
                    {
                        //  slack bus, set GenMW = 0, GenMVAR = 0;
                        GenMW[i] = 0;
                        GenMVAR[i] = 0;
                    }
                    if (BusType[i] == 2)
                    {
                        //  PV bus, set GenMVAR = 0, VolAng = 0;
                        GenMVAR[i] = 0;
                        VolAng[i] = 0;
                    }
                    if (BusType[i] == 1 || BusType[i] == 0)
                    {
                        //  PQ bus, set VolMag = 1, VolAng = 0;
                        VolMag[i] = 1;
                        VolAng[i] = 0;
                    }
                    //  Pinjection = Pgen - Pload, must convert to p.u.!
                    Pi.Add((GenMW[i] - LoadMW[i]) / 100);
                    Qi.Add((GenMVAR[i] - LoadMVAR[i]) / 100);
                }

                //  initialize the iteration count and error.
                int iter = 0;
                double error = 1.0;

                //  Newton Raphson iteration.
                while (error > stoppingError)
                {
                    // generate Jacobian matrix and mismatch deltaP and deltaQ.
                    double[,] Jacobian;
                    List<double> dPQ;
                    CommonFunctions.Jacobian(BusNO, BusType, VolAng, VolMag, Pi, Qi, Gbus, Bbus, out Jacobian, out dPQ);

                    // compute [dVol] by inverse of Jacobian multiply [dPQ].
                    List<double> dVol;
                    CommonFunctions.LUFact(Jacobian, dPQ, out dVol);

                    // update the VolAng and VolMag, 
                    // start another iteration until stopping error is satisfied.
                    foreach (int i in PVPQ)
                    {
                        int index = PVPQ.FindIndex(x => x == i);
                        VolAng[i] = VolAng[i] + dVol[index];
                    }
                    foreach (int i in PQ)
                    {
                        int index = PQ.FindIndex(x => x == i) + PVPQ.Count;
                        VolMag[i] = VolMag[i] + dVol[index];
                    }

                    iter++;

                    // update the error.
                    CommonFunctions.FindMaxAbs(dPQ, out error);
                }

                //  display the power flow results.
                UpdateResults(VolAng, VolMag, iter);

                //  results save to a txt file.
                SaveResults(VolAng, VolMag);
            }
            catch (Exception ex)
            {
                MessageBox.Show("Please check data format! " + ex.Message);
            }
        }

        private void UpdateResults(List<double> VolAng, List<double> VolMag, int iter)
        {
            //  display the results in the main window.
            //  only display the first 10 voltage magnitudes and angles.

            textBuses.Text = (bus[0].Count).ToString();
            textIter.Text = (iter).ToString();
            textBox1.Text = (VolAng[0]).ToString("0.00");
            textBox2.Text = (VolAng[1]).ToString("0.00");
            textBox3.Text = (VolAng[2]).ToString("0.00");
            textBox4.Text = (VolAng[3]).ToString("0.00");
            textBox5.Text = (VolAng[4]).ToString("0.00");
            textBox6.Text = (VolAng[5]).ToString("0.00");
            textBox7.Text = (VolAng[6]).ToString("0.00");
            textBox8.Text = (VolAng[7]).ToString("0.00");
            textBox9.Text = (VolAng[8]).ToString("0.00");
            textBox10.Text = (VolAng[9]).ToString("0.00");
            textBox11.Text = (VolMag[0]).ToString("0.00");
            textBox12.Text = (VolMag[1]).ToString("0.00");
            textBox13.Text = (VolMag[2]).ToString("0.00");
            textBox14.Text = (VolMag[3]).ToString("0.00");
            textBox15.Text = (VolMag[4]).ToString("0.00");
            textBox16.Text = (VolMag[5]).ToString("0.00");
            textBox17.Text = (VolMag[6]).ToString("0.00");
            textBox18.Text = (VolMag[7]).ToString("0.00");
            textBox19.Text = (VolMag[8]).ToString("0.00");
            textBox20.Text = (VolMag[9]).ToString("0.00");
        }

        private void SaveResults(List<double> VolAng, List<double> VolMag)
        {
            //  save the full results to a txt file named "results.stg".

            StreamWriter settingfile = new StreamWriter(".\\results.stg");
            settingfile.WriteLine("Voltage Angle:");
            for (int i = 0; i < VolAng.Count; ++i)
            {
                settingfile.WriteLine(VolAng[i]);
            }
            settingfile.WriteLine("Voltage Magnitude:");
            for (int i = 0; i < VolMag.Count; ++i)
            {
                settingfile.WriteLine(VolMag[i]);
            }
            settingfile.Close();
        }

        private void clear_Click(object sender, EventArgs e)
        {
            //  clear bus data and branch data when user clicks the button.

            bus = new List<List<string>>();
            branch = new List<List<string>>();
            MessageBox.Show("Data is cleared.");
        }

        private void updateSetting_Click(object sender, EventArgs e)
        {
            //  set the stopping error according to user's input.

            StreamWriter settingfile = new StreamWriter(".\\settings.stg");
            try
            {
                settingfile.WriteLine(double.Parse(errorSetting.Text));
                stoppingError = double.Parse(errorSetting.Text);
            }
            catch
            {
                settingfile.WriteLine("0.01");
                stoppingError = 0.01;
            }
            settingfile.Close();
        }

        private void defaultSetting_Click(object sender, EventArgs e)
        {
            //  the default value for stopping error is 0.01.

            errorSetting.Text = "0.01";
        }

        private void FastDecoupled_Click(object sender, EventArgs e)
        {
            MessageBox.Show("This part is under construction. Thank you for your patience.");
        }
        
    }
}
