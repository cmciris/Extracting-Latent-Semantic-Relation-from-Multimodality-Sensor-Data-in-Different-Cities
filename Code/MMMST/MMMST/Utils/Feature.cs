using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace MMMST.Utils
{
    public class Feature
    {
        /// <summary>
        /// normalize the feature matrix
        /// </summary>
        /// <param name="featureMatrix"></param>specifies the feature matrix to be normalized
        /// <param name="normType"></param>specifies the type of normalization 0 - standardization 1 - max-min scaling 2- normalization 
        /// <returns></returns>
        public static List<double>[] NormalizeFeature(List<double>[] featureMatrix, int normType)
        {
            int nSamples = featureMatrix.Length;
            int dim = 0;
            for (int i = 0; i < nSamples; i++)
            {
                if (featureMatrix[i] != null)
                {
                    if (featureMatrix[i].Count() > 1)
                        dim = featureMatrix[i].Count();
                }
                else
                {
                    featureMatrix[i] = new List<double>();
                    featureMatrix[i].Add(double.MaxValue);
                }
            }

            List<double>[] featureMatrixTrans = new List<double>[dim];
            for (int i = 0; i < dim; i++)
            {
                featureMatrixTrans[i] = new List<double>();
            }

            for (int i = 0; i < nSamples; i++)
            {
                if (featureMatrix[i].Sum() == 0)
                {
                    featureMatrix[i].Clear();
                    featureMatrix[i].Add(Double.MaxValue);
                }
                double norm = 0.0; //normalization
                if (featureMatrix[i].Count() > 1)
                {
                    for (int j = 0; j < dim; j++)
                    {
                        featureMatrixTrans[j].Add(featureMatrix[i][j]);
                        norm += Math.Pow(featureMatrix[i][j], 2);
                    }
                    norm = Math.Sqrt(norm);
                    if (normType == 2)
                    {
                        for (int j = 0; j < dim; j++)
                        {
                            if (norm > 0)
                            {
                                featureMatrix[i][j] /= norm;
                            }
                            else
                            {
                                featureMatrix[i][j] = 0;
                            }
                        }
                    }
                }
            }

            //read in the mean and std
            StreamReader srMean = new StreamReader(@"D:\Dropbox\Dropbox\Project\MSRA\Result\Features\mean.txt");
            StreamReader srStd = new StreamReader(@"D:\Dropbox\Dropbox\Project\MSRA\Result\Features\std.txt");
            List<double> meanR = new List<double>();
            List<double> stdR = new List<double>();
            while (srMean.Peek() > 0)
            {
                string[] items = srMean.ReadLine().Split(new char[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
                if (items.Count() == dim)
                {
                    foreach (var it in items)
                    {
                        meanR.Add(Convert.ToDouble(it));
                    }
                }
            }
            while (srStd.Peek() > 0)
            {
                string[] items = srStd.ReadLine().Split(new char[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
                if (items.Count() == dim)
                {
                    foreach (var it in items)
                    {
                        stdR.Add(Convert.ToDouble(it));
                    }
                }
            }
            srMean.Close();
            srStd.Close();

            for (int i = 0; i < dim; i++)
            {
                double max = featureMatrixTrans[i].Max();
                double min = featureMatrixTrans[i].Min();
                //double mean = featureMatrixTrans[i].Average();
                //double std = CalculateStdDev(featureMatrixTrans[i]);

                double mean = meanR[i];
                double std = stdR[i];


                if (featureMatrixTrans[i].Count() > 1)
                {
                    for (int j = 0; j < nSamples; j++)
                    {
                        if (featureMatrix[j].Count() > 1)
                        {
                            if (normType == 0)
                            {
                                if (std > 0)
                                {
                                    featureMatrix[j][i] = (featureMatrix[j][i] - mean) / std;
                                }
                                else
                                {
                                    featureMatrix[j][i] = 0;
                                }
                            }
                            else if (normType == 1)
                            {
                                if ((max - min) > 0)
                                {
                                    featureMatrix[j][i] = (featureMatrix[j][i] - min) / (max - min);
                                }
                                else
                                {
                                    featureMatrix[j][i] = 0;
                                }
                            }
                        }
                    }
                }
            }

            return featureMatrix;
        }

        private static double CalculateStdDev(IEnumerable<double> values)
        {
            double ret = 0;
            if (values.Count() > 0)
            {
                //Compute the Average      
                double avg = values.Average();
                //Perform the Sum of (value-avg)_2_2      
                double sum = values.Sum(d => Math.Pow(d - avg, 2));
                //Put it all together      
                ret = Math.Sqrt((sum) / (values.Count() - 1));
            }
            return ret;
        }
    }
}
