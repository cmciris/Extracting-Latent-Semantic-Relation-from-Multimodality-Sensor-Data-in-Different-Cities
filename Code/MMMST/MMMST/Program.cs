using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MMMST.Preprocess;
using System.IO;

using System.Diagnostics;

namespace MMMST
{
    public class Program
    {
        static void Main(string[] args)
        {
            #region Step 1: Generate Features of Different Modalities from Raw Data
            //Methods.GenerateFeature(@"E:\Dropbox\Project\MSRA\Data\Beijing\roadnetwork\Road_Network_2014.txt", @"D:\Dropbox\Dropbox\Project\MSRA\Result\Features3\Beijing\roadNetworkFeature_0.015.txt",
            //                        @"E:\Dropbox\Project\MSRA\Data\Beijing\poi\poi.csv", @"D:\Dropbox\Dropbox\Project\MSRA\Result\Features3\Beijing\cityPOIFeature_0.015.txt", @"D:\Dropbox\Dropbox\Project\MSRA\Data\Beijing\roadnetwork\road_remap_2014.txt",
            //                        @"E:\Dropbox\Project\MSRA\Data\Beijing\mobility\speed", @"D:\Dropbox\Dropbox\Project\MSRA\Data\Beijing\mobility\PickDropPoints", @"D:\Dropbox\Dropbox\Project\MSRA\Result\Features3\Beijing\mobilityFeature\original\timeSlot24\0.015\",
            //                        @"E:\Dropbox\Project\MSRA\Data\Beijing\new meterolgy\meterology.csv", @"D:\Dropbox\Dropbox\Project\MSRA\Data\Beijing\new meterolgy\station.csv", @"D:\Dropbox\Dropbox\Project\MSRA\Result\Features3\Beijing\meterologyFeature\mode0\0.015\",
            //                        @"E:\Dropbox\Project\MSRA\Data\Beijing\aqi\airquality.csv", @"D:\Dropbox\Dropbox\Project\MSRA\Data\Beijing\aqi\aqi.txt", @"D:\Dropbox\Dropbox\Project\MSRA\Data\Beijing\aqi\aqi station.txt", @"D:\Dropbox\Dropbox\Project\MSRA\Result\Features3\Beijing\aqiFeature\0.015\", 0);
            #endregion

            #region Step 2: Construct the Graph between the Generated Features of Differet Modalities
            //for (int knn = 50; knn <= 50; knn = knn + 10)
            //{
            //    for (double sigma = 5; sigma <= 5; sigma += 1)
            //    {
            //        for (int timeSlot = 0; timeSlot <= 23; timeSlot++)
            //        {
            //            string roadNetworkFilename = @"E:\Dropbox\Project\MSRA\Data\Beijing\roadnetwork\Road_Network_2014.txt";
            //            string roadNetworkFeature = @"E:\Dropbox\Project\MSRA\Result\Features\Beijing\roadNetworkFeature_0.015.txt";
            //            string poiFeature = @"E:\Dropbox\Project\MSRA\Result\Features\Beijing\cityPOIFeature_0.015.txt";
            //            string meterologyFeature = @"E:\Dropbox\Project\MSRA\Result\Features\Beijing\meterologyFeature\mode0\0.015\";
            //            string mobilityFeature = @"E:\Project\MSRA\Result\Features\Beijing\mobilityFeature\original\timeSlot24\0.015\";
            //            string aqiFeature = @"E:\Project\MSRA\Result\Features\Beijing\aqiFeature\0.015\";
            //            string graphDirectory = @"E:\Project\MSRA\Result\Temporary\Graph\100\inter-1\" + knn.ToString() + @"\" + sigma.ToString();

            //            string dir = Path.Combine(graphDirectory, timeSlot.ToString("D2"));
            //            if (!Directory.Exists(dir))
            //            {
            //                Directory.CreateDirectory(dir);
            //            }
            //            Methods.ConstructGraph(roadNetworkFilename, roadNetworkFeature, poiFeature, meterologyFeature, mobilityFeature, aqiFeature, Path.Combine(dir, "graph_labeled.txt"), 100, knn, timeSlot, sigma);
            //        }
            //    }
            //}
            #endregion

            #region Coupled dictionary learning : adjust the # of clusters
            //StreamWriter swObjective = new StreamWriter(@"E:\Dropbox\Project\MSRA\Result\Temporary\Objective\010\inter-1\objective_0.txt");
            //for (int timeSlot = 11; timeSlot <= 11; timeSlot++) //the specific timeslot of graph to perform graph clustering
            //{
            //    for (int knn = 50; knn <= 50; knn = knn + 10) //k-nearest neighbours
            //    {
            //        for (double sigma = 5; sigma <= 5; sigma += 1) //sigma
            //        {
            //            for (int nc = 100; nc <=100; nc = nc + 100) //number of clusters
            //            {
            //                for (double lambda = 500; lambda <= 500; lambda = lambda * 10)
            //                {
            //                    for (double gamma = 100; gamma <= 100; gamma = gamma * 10)
            //                    {
            //                        for (double mu = 100; mu <= 100; mu = mu * 10)
            //                        {
            //                            string graphDirectory = Path.Combine(@"E:\Dropbox\Project\MSRA\Result\Temporary\Graph\100\inter-1\", knn.ToString(), sigma.ToString(), timeSlot.ToString("D2"));
            //                            string clusterDirectory = @"E:\Dropbox\Project\MSRA\Result\Temporary\Cluster\Labeled\0.015\100\inter-1\";
            //                            string dictDirectory = Path.Combine(@"E:\Dropbox\Project\MSRA\Result\Temporary\Dict\Labeled\0.015\100\inter-1\", "knn_" + knn.ToString(), "sigma_" + sigma.ToString(), "nc_" + nc.ToString(), "lambda_" + lambda.ToString(), "gamma_" + gamma.ToString(), "mu_" + mu.ToString(), timeSlot.ToString("D2"));

            //                            string dir = Path.Combine(clusterDirectory, "knn_" + knn.ToString(), "sigma_" + sigma.ToString(), "nc_" + nc.ToString(), "lambda_" + lambda.ToString(), "gamma_" + gamma.ToString(), "mu_" + mu.ToString(), timeSlot.ToString("D2"));

            //                            if (!Directory.Exists(dir))
            //                            {
            //                                Directory.CreateDirectory(dir);
            //                            }
            //                            if (!Directory.Exists(dictDirectory))
            //                            {
            //                                Directory.CreateDirectory(dictDirectory);
            //                            }
            //                            double entropy, purity, balance, diversity, objective;
            //                            Methods.CoupledDictionaryLearning(Path.Combine(graphDirectory, "graph_labeled.txt"), dir, nc, lambda, gamma, mu, out entropy, out purity, out balance, out diversity, out objective);
            //                            CoupledDL.LazyGreedy.InferDictionary(dir, nc, dictDirectory);
            //                            swObjective.WriteLine("knn = " + knn + " nc = " + nc + " sigma = " + sigma + " lambda" + lambda + " gamma" + gamma + " : " + mu + " : " + entropy + " " + purity + " " + balance + " " + objective + " " + (entropy + balance + purity + diversity));
            //                        }
            //                    }
            //                }
            //            }
            //        }
            //    }
            //}
            //swObjective.Close();
            #endregion

            ////Utils.Visualization.VisualizationClusters(@"D:\Dropbox\Dropbox\Project\MSRA\Result\Temporary\Dict\inter-1\", 30);

            #region Retrieve Only Labeled Subset of Generated Features
            //for (int ts = 0; ts < 24; ts++)
            //{
            //    LabelSubset.RetrieveLabelSubset(@"D:\Dropbox\Dropbox\Project\MSRA\Data\Tianjin", @"D:\Dropbox\Dropbox\Project\MSRA\Result\Features\Tianjin", ts, 1);
            //}
            #endregion
        }
    }

    public class Methods
    {
        /// <summary>
        /// generate all the features of all the data sources 
        /// mode = 0, meterology data is found by knn; mode = 1, only find the stations' region
        /// </summary>
        public static void GenerateFeature(string roadNetworkFilename, string roadNetworkOutputFilename,
                                    string cityPOIFilename, string cityPOIOutputFilename,
                                    string roadMapFilename, string speedDirectory, string pickDropPointsDirectory, string mobilityOutputDirectory,
                                    string meterologyFilename, string stationLocationFilename, string meterologyFeatureDirectory,
                                    string airQualityFilename, string aqiFilename, string airStationFilename, string aqiFeatureDirectory, int mode)
        {
            RoadNetwork roadNetwork = IO.ReadRoadNetwork(roadNetworkFilename);
            CityPOI cityPOI = IO.ReadCityPOI(cityPOIFilename);

            #region Solution: just use the grid of the road network
            Grid grid = new Grid(roadNetwork.MinPoint, roadNetwork.MaxPoint, 0.015);
            #endregion 

            #region roadnetwork feature
            roadNetwork.RoadFeatures(grid, roadNetworkOutputFilename);
            #endregion
            #region cityPOI
            cityPOI.CityPOIFeature(grid, cityPOIOutputFilename);
            #endregion
            #region totalMobility feature
            DirectoryInfo dirSpeed = new DirectoryInfo(speedDirectory);
            FileInfo[] files1 = dirSpeed.GetFiles();
            DirectoryInfo dirPoints = new DirectoryInfo(pickDropPointsDirectory);
            FileInfo[] files2 = dirPoints.GetFiles();
            int countDay = files1.Count();

            for (int i = 0; i < countDay; i++)
            {
                Console.WriteLine(i + "day is readed.");
                using (ToTalMobility totalMobility = IO.ReadMobility(roadNetwork, roadMapFilename, files1[i].FullName, files2[i].FullName, files1[i].Name))
                {
                    totalMobility.TotalMobilityFeature(grid, mobilityOutputDirectory);
                }
            }
            #endregion
            #region meterology feature
            TotalMeterology totalMeterology = IO.ReadMeterology(meterologyFilename, stationLocationFilename);
            totalMeterology.TotalMeterologyFeature(grid, meterologyFeatureDirectory, mode);// mode = 0, meterology data is found by knn; mode = 1, only find the stations' region
            #endregion
            #region AQI label
            TotalAirQuality totalAirQuality = IO.ReadAirQuality(airQualityFilename, aqiFilename, airStationFilename);
            totalAirQuality.TotalAirQualityFeature(grid, aqiFeatureDirectory);
            #endregion
        }
   
        public static void ConstructGraph(string roadFilename, string roadNetworkFilename, string poiFilename, string meterologyDirectory, string mobilityDirectory, string aqiDirectory, string outputGraphFilename, int numDays, int knn, int timeSlot, double s)
        {
            RoadNetwork roadNetwork = IO.ReadRoadNetwork(roadFilename);
            Grid grid = new Grid(roadNetwork.MinPoint, roadNetwork.MaxPoint, 0.015);

            CoupledDL.GraphInput graphInput = new CoupledDL.GraphInput(timeSlot, numDays, 4, roadNetworkFilename, poiFilename, meterologyDirectory, mobilityDirectory, aqiDirectory, grid);
            CoupledDL.Graph.HyperGraph graph = graphInput.ReadGraphFromFiles(1);
            double[] sigma = new double[] { 10, 2, 4, 4 };
            graphInput.BuildEdgesOfGraph(graph, knn, 9, sigma, 0, outputGraphFilename);
        }
    
        public static void CoupledDictionaryLearning(string graphFilename, string outputClusterDirectory, int numClusters, double lambda, double gamma, double mu, out double entropy, out double purity, out double balance, out double diversity, out double objective)
        {
            //read in the graph
            CoupledDL.Graph.HyperGraph inGraph = CoupledDL.Graph.HyperGraph.ReadGraph(graphFilename);
            //perform clustering
            List<CoupledDL.Graph.CompositeNode> clusters = CoupledDL.LazyGreedy.Clustering(inGraph, lambda, gamma, mu, numClusters, outputClusterDirectory, out entropy, out purity, out balance, out diversity, out objective);
            //check the clustering result
            //CoupledDL.LazyGreedy.OutputClusterResult(clusters, inGraph, outputClusterDirectory, 1);
        }
    }
}
