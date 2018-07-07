using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace MMMST.Utils
{
    public class Visualization
    {
        public static void VisualizationClusters(string clusterDirectory, int numClusters)
        {
            //create the MATLAB instance 
            MLApp.MLApp matlab = new MLApp.MLApp();

            //Change to the direcoty where the function is located
            matlab.Execute(@"cd D:\Project\Matlab"); 

            DirectoryInfo dir = new DirectoryInfo(clusterDirectory);
            FileInfo[] files = dir.GetFiles();

            double[,] hists = new double[numClusters, files.Count()];
            int count = 0;
            foreach (var f in files)
            {
                //Define the output
                object values = null;

                //Call the MATLAB function VisualizationClusters
                matlab.Feval("VisualizeClusters", 3, out values, f.FullName);

                object[] val = values as object[];
                double[,] value = val[1] as double[,];
                double[,] hist = val[2] as double[,];

                List<double> listValue = new List<double>();
                foreach(var v in value)
                {
                    listValue.Add(v);
                }
                List<double> listHist = new List<double>();
                foreach(var h in hist)
                {
                    listHist.Add(h);
                }

                for (int i = 0; i < numClusters; i++)
                {
                    if(listValue.Contains(i))
                    {
                        int index = listValue.IndexOf(i);
                        hists[i, count] = listHist[index];
                    }
                    else
                    {
                        hists[i, count] = 0;
                    }
                }
                count++;
            }

            Object bar_handle = null;
            Object result0, result1, result2, result3 = null;
            matlab.Feval("bar", 1, out bar_handle, hists, "grouped");
            Object[] handle = bar_handle as object[];
            double[,] handles = handle[0] as double[,];
            matlab.Feval("set", 0, out result0, handles[0, 0], "FaceColor", "b");
            matlab.Feval("set", 0, out result1, handles[0, 1], "FaceColor", "g");
            matlab.Feval("set", 0, out result2, handles[0, 2], "FaceColor", "r");
            matlab.Feval("set", 0, out result3, handles[0, 3], "FaceColor", "c");
        }
    }
}
