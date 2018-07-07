using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace MMMST.CoupledDL
{
    public class GraphInput
    {
        public int TimeSlot;
        public int NumModality;
        public string RoadNetworkFilename;
        public string CityPOIFilename;
        public string MeterologyDirectory;
        public string MobilityDirectory;
        public string AQIDirectory;
        public Preprocess.Grid Grid;
        public int GridNum;
        public List<string> Days = new List<string>();
        public int NumDays;

        public GraphInput(int ts, int nd, int nm, string rnFilename, string poiFilename, string mtd, string md, string ad, Preprocess.Grid g)
        {
            this.TimeSlot = ts;
            this.NumDays = nd;
            this.NumModality = nm;
            this.RoadNetworkFilename = rnFilename;
            this.CityPOIFilename = poiFilename;
            this.MeterologyDirectory = mtd;
            this.MobilityDirectory = md;
            this.AQIDirectory = ad;

            this.Grid = new Preprocess.Grid(g.Min, g.Max, g.GridSize);
            this.GridNum = g.GridNumLat * g.GridNumLon;
            SetNumDays(MeterologyDirectory);
        }

        public Graph.HyperGraph ReadGraphFromFiles(int normType)
        {
            //read all the modality nodes in
            //1. read in the road network
            List<double>[] roadNetworkFeature = new List<double>[this.GridNum];
            roadNetworkFeature = ReadFeature(this.RoadNetworkFilename);
            roadNetworkFeature = Utils.Feature.NormalizeFeature(roadNetworkFeature, normType);
            //2. read in the poi feature
            List<double>[] cityPOIFeature = new List<double>[this.GridNum];
            cityPOIFeature = ReadFeature(this.CityPOIFilename);
            cityPOIFeature = Utils.Feature.NormalizeFeature(cityPOIFeature, normType);
            //3. read in the meterology feature in a number of days
            DirectoryInfo dir = new DirectoryInfo(this.MeterologyDirectory);
            DirectoryInfo[] subDirs = dir.GetDirectories();
            List<double>[] meterologyFeature = new List<double>[this.GridNum * NumDays];
            //4. read in the mobility feature in a number of days
            DirectoryInfo dir2 = new DirectoryInfo(this.MobilityDirectory);
            DirectoryInfo[] subDirs2 = dir2.GetDirectories();
            List<double>[] mobilityFeature = new List<double>[this.GridNum * NumDays];
            //5. read in the aqi label in a number of days 
            DirectoryInfo dir3 = new DirectoryInfo(this.AQIDirectory);
            DirectoryInfo[] subDirs3 = dir3.GetDirectories();
            List<double>[] aqiLabel = new List<double>[this.GridNum * NumDays];

            int count1 = 0, count2 = 0;
            while (true)
            {
                string day = this.Days[count1];

                int indexDir1 = 0, indexDir2 = 0, indexDir3 = 0;
                bool find1 = false, find2 = false, find3 = false;
                foreach (var n in subDirs)
                {
                    if (n.Name == day)
                    {
                        find1 = true;
                        break;
                    }
                    indexDir1++;
                }

                foreach (var n in subDirs2)
                {
                    if (n.Name == day)
                    {
                        find2 = true;
                        break;
                    }
                    indexDir2++;
                }

                foreach (var n in subDirs3)
                {
                    if (n.Name == day)
                    {
                        find3 = true;
                        break;
                    }
                    indexDir3++;
                }

                if (find1 && find2 && find3)
                {
                    FileInfo[] files = subDirs[indexDir1].GetFiles();
                    FileInfo[] files2 = subDirs2[indexDir2].GetFiles();
                    FileInfo[] files3 = subDirs3[indexDir3].GetFiles();

                    foreach (var f in files)
                    {
                        if (Convert.ToInt32(f.Name) == this.TimeSlot)
                        {
                            List<double>[] tmp = new List<double>[this.GridNum];
                            tmp = ReadFeature(f.FullName);
                            for (int j = 0; j < this.GridNum; j++)
                            {
                                meterologyFeature[count2 * this.GridNum + j] = tmp[j]; //if = null indicates that the timeslot does not have any data.
                            }
                        }
                    }

                    foreach (var f in files2)
                    {
                        if (Convert.ToInt32(f.Name) == this.TimeSlot)
                        {
                            List<double>[] tmp = new List<double>[this.GridNum];
                            tmp = ReadFeature(f.FullName);
                            for (int j = 0; j < this.GridNum; j++)
                            {
                                mobilityFeature[count2 * this.GridNum + j] = tmp[j];
                                if((count2 * this.GridNum + j) == 8896)
                                {
                                    Console.WriteLine();
                                }
                            }
                        }
                    }

                    foreach (var f in files3)
                    {
                        if (Convert.ToInt32(f.Name) == this.TimeSlot)
                        {
                            List<double>[] tmp = new List<double>[this.GridNum];
                            tmp = ReadFeature(f.FullName);
                            for (int j = 0; j < this.GridNum; j++)
                            {
                                aqiLabel[count2 * this.GridNum + j] = tmp[j];
                            }
                        }
                    }
                    count2++;
                }
                count1++;
                if (count2 >= this.NumDays)
                    break;
            }

            meterologyFeature = Utils.Feature.NormalizeFeature(meterologyFeature, normType);
            mobilityFeature = Utils.Feature.NormalizeFeature(mobilityFeature, normType);

            int[] aqi = AQILabelInference(aqiLabel);

            List<Graph.SingleNode>[] nodes = new List<Graph.SingleNode>[this.NumModality];
            for (int i = 0; i < this.NumModality; i++)
            {
                nodes[i] = new List<Graph.SingleNode>();
            }
            int nodeID = 0;
            for (int i = 0; i < this.GridNum; i++)
            {
                if (roadNetworkFeature[i].Count() > 1) // exclude those regions which do not have road network features
                {
                    nodes[0].Add(new Graph.SingleNode(nodeID++, roadNetworkFeature[i], -1, 1, i, -1)); //modality = 1 - roadnetwork
                }
            }
            for (int i = 0; i < this.GridNum; i++)
            {
                if (cityPOIFeature[i].Count() > 1)
                {
                    nodes[1].Add(new Graph.SingleNode(nodeID++, cityPOIFeature[i], -1, 2, i, -1)); //modality = 2 - poi
                }
            }
            for (int i = 0; i < this.NumDays * this.GridNum; i++)
            {
                if (meterologyFeature[i].Count() > 1 && aqi[i] > 0)
                {
                    nodes[2].Add(new Graph.SingleNode(nodeID++, meterologyFeature[i], aqi[i], 3, i % this.GridNum, i / this.GridNum)); //modality = 3 -meterology
                }
            }
            for (int i = 0; i < this.NumDays * this.GridNum; i++)
            {
                if (mobilityFeature[i].Count() > 1 && aqi[i] > 0)
                {
                    nodes[3].Add(new Graph.SingleNode(nodeID++, mobilityFeature[i], aqi[i], 4, i % this.GridNum, i / this.GridNum)); // modality = 4 - mobility
                }
            }
            List<Graph.Edge> edges = new List<Graph.Edge>();

            return new Graph.HyperGraph(nodes, edges);
        }


        public Graph.HyperGraph ReadGraphFromFiles2(int normType)
        {
            //read all the modality nodes in
            //1. read in the road network
            List<double>[] roadNetworkFeature = new List<double>[this.GridNum];
            roadNetworkFeature = ReadFeature(this.RoadNetworkFilename);
            roadNetworkFeature = Utils.Feature.NormalizeFeature(roadNetworkFeature, normType);
            //2. read in the poi feature
            List<double>[] cityPOIFeature = new List<double>[this.GridNum];
            cityPOIFeature = ReadFeature(this.CityPOIFilename);
            cityPOIFeature = Utils.Feature.NormalizeFeature(cityPOIFeature, normType);
            //3. read in the meterology feature in a number of days
            DirectoryInfo dir = new DirectoryInfo(this.MeterologyDirectory);
            DirectoryInfo[] subDirs = dir.GetDirectories();
            List<double>[] meterologyFeature = new List<double>[this.GridNum * NumDays];

            //5. read in the aqi label in a number of days 
            DirectoryInfo dir3 = new DirectoryInfo(this.AQIDirectory);
            DirectoryInfo[] subDirs3 = dir3.GetDirectories();
            List<double>[] aqiLabel = new List<double>[this.GridNum * NumDays];

            int count1 = 0, count2 = 0;
            while (true)
            {
                string day = this.Days[count1];

                int indexDir1 = 0, indexDir3 = 0;
                bool find1 = false,  find3 = false;
                foreach (var n in subDirs)
                {
                    if (n.Name == day)
                    {
                        find1 = true;
                        break;
                    }
                    indexDir1++;
                }

                foreach (var n in subDirs3)
                {
                    if (n.Name == day)
                    {
                        find3 = true;
                        break;
                    }
                    indexDir3++;
                }

                if (find1 && find3)
                {
                    FileInfo[] files = subDirs[indexDir1].GetFiles();
                    FileInfo[] files3 = subDirs3[indexDir3].GetFiles();

                    foreach (var f in files)
                    {
                        if (Convert.ToInt32(f.Name) == this.TimeSlot)
                        {
                            List<double>[] tmp = new List<double>[this.GridNum];
                            tmp = ReadFeature(f.FullName);
                            for (int j = 0; j < this.GridNum; j++)
                            {
                                meterologyFeature[count2 * this.GridNum + j] = tmp[j]; //if = null indicates that the timeslot does not have any data.
                            }
                        }
                    }

                    foreach (var f in files3)
                    {
                        if (Convert.ToInt32(f.Name) == this.TimeSlot)
                        {
                            List<double>[] tmp = new List<double>[this.GridNum];
                            tmp = ReadFeature(f.FullName);
                            for (int j = 0; j < this.GridNum; j++)
                            {
                                aqiLabel[count2 * this.GridNum + j] = tmp[j];
                            }
                        }
                    }
                    count2++;
                }
                count1++;
                if (count2 >= this.NumDays)
                    break;
            }

            meterologyFeature = Utils.Feature.NormalizeFeature(meterologyFeature, normType);

            int[] aqi = AQILabelInference(aqiLabel);

            List<int> indexLabel = new List<int>();
            for (int i = 0; i < this.GridNum * NumDays; i++)
            {
                if (aqi[i] != -1 && !indexLabel.Contains(i % this.GridNum))
                    indexLabel.Add(i % this.GridNum);
            }


            List<Graph.SingleNode>[] nodes = new List<Graph.SingleNode>[this.NumModality];
            for (int i = 0; i < this.NumModality; i++)
            {
                nodes[i] = new List<Graph.SingleNode>();
            }
            int nodeID = 0;
            for (int i = 0; i < this.GridNum; i++)
            {
                if (roadNetworkFeature[i].Count() > 1 && indexLabel.Contains(i)) // exclude those regions which do not have road network features
                {
                    nodes[0].Add(new Graph.SingleNode(nodeID++, roadNetworkFeature[i], -1, 1, i, -1)); //modality = 1 - roadnetwork
                }
            }
            for (int i = 0; i < this.GridNum; i++)
            {
                if (cityPOIFeature[i].Count() > 1 && indexLabel.Contains(i))
                {
                    nodes[1].Add(new Graph.SingleNode(nodeID++, cityPOIFeature[i], -1, 2, i, -1)); //modality = 2 - poi
                }
            }
            for (int i = 0; i < this.NumDays * this.GridNum; i++)
            {
                if (meterologyFeature[i].Count() > 1 && indexLabel.Contains(i % this.GridNum) && aqi[i] > 0)
                {
                    nodes[2].Add(new Graph.SingleNode(nodeID++, meterologyFeature[i], aqi[i], 3, i % this.GridNum, i / this.GridNum)); //modality = 3 -meterology
                }
            }

            List<Graph.Edge> edges = new List<Graph.Edge>();

            return new Graph.HyperGraph(nodes, edges);
        }


        public Graph.HyperGraph ReadGraphFromFiles3(int normType)
        {
            //read all the modality nodes in
            //1. read in the road network
            List<double>[] roadNetworkFeature = new List<double>[this.GridNum];
            roadNetworkFeature = ReadFeature(this.RoadNetworkFilename);
                roadNetworkFeature = Utils.Feature.NormalizeFeature(roadNetworkFeature, normType);
            //2. read in the poi feature
            List<double>[] cityPOIFeature = new List<double>[this.GridNum];
            cityPOIFeature = ReadFeature(this.CityPOIFilename);
                cityPOIFeature = Utils.Feature.NormalizeFeature(cityPOIFeature, normType);
            //3. read in the meterology feature in a number of days
            DirectoryInfo dir = new DirectoryInfo(this.MeterologyDirectory);
            DirectoryInfo[] subDirs = dir.GetDirectories();
            List<double>[] meterologyFeature = new List<double>[this.GridNum * NumDays];
            //4. read in the mobility feature in a number of days
            DirectoryInfo dir2 = new DirectoryInfo(this.MobilityDirectory);
            DirectoryInfo[] subDirs2 = dir2.GetDirectories();
            List<double>[] mobilityFeature = new List<double>[this.GridNum * NumDays];
            //5. read in the aqi label in a number of days 
            DirectoryInfo dir3 = new DirectoryInfo(this.AQIDirectory);
            DirectoryInfo[] subDirs3 = dir3.GetDirectories();
            List<double>[] aqiLabel = new List<double>[this.GridNum * NumDays];

            int count1 = 0, count2 = 0;
            while (true)
            {
                string day = this.Days[count1];

                int indexDir1 = 0, indexDir2 = 0, indexDir3 = 0;
                bool find1 = false, find2 = false, find3 = false;
                foreach (var n in subDirs)
                {
                    if (n.Name == day)
                    {
                        find1 = true;
                        break;
                    }
                    indexDir1++;
                }

                foreach (var n in subDirs2)
                {
                    if (n.Name == day)
                    {
                        find2 = true;
                        break;
                    }
                    indexDir2++;
                }

                foreach (var n in subDirs3)
                {
                    if (n.Name == day)
                    {
                        find3 = true;
                        break;
                    }
                    indexDir3++;
                }

                if (find1 && find2 && find3)
                {
                    FileInfo[] files = subDirs[indexDir1].GetFiles();
                    FileInfo[] files2 = subDirs2[indexDir2].GetFiles();
                    FileInfo[] files3 = subDirs3[indexDir3].GetFiles();

                    foreach (var f in files)
                    {
                        if (Convert.ToInt32(f.Name) == this.TimeSlot)
                        {
                            List<double>[] tmp = new List<double>[this.GridNum];
                            tmp = ReadFeature(f.FullName);
                            for (int j = 0; j < this.GridNum; j++)
                            {
                                meterologyFeature[count2 * this.GridNum + j] = tmp[j]; //if = null indicates that the timeslot does not have any data.
                            }
                        }
                    }

                    foreach (var f in files2)
                    {
                        if (Convert.ToInt32(f.Name) == this.TimeSlot)
                        {
                            List<double>[] tmp = new List<double>[this.GridNum];
                            tmp = ReadFeature(f.FullName);
                            for (int j = 0; j < this.GridNum; j++)
                            {
                                mobilityFeature[count2 * this.GridNum + j] = tmp[j];
                                if ((count2 * this.GridNum + j) == 8896)
                                {
                                    Console.WriteLine();
                                }
                            }
                        }
                    }

                    foreach (var f in files3)
                    {
                        if (Convert.ToInt32(f.Name) == this.TimeSlot)
                        {
                            List<double>[] tmp = new List<double>[this.GridNum];
                            tmp = ReadFeature(f.FullName);
                            for (int j = 0; j < this.GridNum; j++)
                            {
                                aqiLabel[count2 * this.GridNum + j] = tmp[j];
                            }
                        }
                    }
                    count2++;
                }
                count1++;
                if (count2 >= this.NumDays)
                    break;
            }

                meterologyFeature = Utils.Feature.NormalizeFeature(meterologyFeature, normType);
                mobilityFeature = Utils.Feature.NormalizeFeature(mobilityFeature, normType);
            

            int[] aqi = AQILabelInference(aqiLabel);
            List<int> indexLabel = new List<int>();
            for (int i = 0; i < this.GridNum * NumDays; i++)
            {
                if (aqi[i] != -1 && !indexLabel.Contains(i % this.GridNum))
                    indexLabel.Add(i % this.GridNum);
            }

            List<Graph.SingleNode>[] nodes = new List<Graph.SingleNode>[this.NumModality];
            for (int i = 0; i < this.NumModality; i++)
            {
                nodes[i] = new List<Graph.SingleNode>();
            }
            int nodeID = 0;
            for (int i = 0; i < this.GridNum; i++)
            {
                if (roadNetworkFeature[i].Count() > 1 && indexLabel.Contains(i)) // exclude those regions which do not have road network features
                {
                    nodes[0].Add(new Graph.SingleNode(nodeID++, roadNetworkFeature[i], -1, 1, i, -1)); //modality = 1 - roadnetwork
                }
            }
            for (int i = 0; i < this.GridNum; i++)
            {
                if (cityPOIFeature[i].Count() > 1 && indexLabel.Contains(i))
                {
                    nodes[1].Add(new Graph.SingleNode(nodeID++, cityPOIFeature[i], -1, 2, i, -1)); //modality = 2 - poi
                }
            }
            for (int i = 0; i < this.NumDays * this.GridNum; i++)
            {
                if (meterologyFeature[i].Count() > 1 && indexLabel.Contains(i % this.GridNum) && aqi[i] > 0 && mobilityFeature[i].Count() > 1)
                {
                    nodes[2].Add(new Graph.SingleNode(nodeID++, meterologyFeature[i], aqi[i], 3, i % this.GridNum, i / this.GridNum)); //modality = 3 -meterology
                    nodes[3].Add(new Graph.SingleNode(nodeID++, mobilityFeature[i], aqi[i], 4, i % this.GridNum, i / this.GridNum)); // modality = 4 - mobility
                }
            }

            List<Graph.Edge> edges = new List<Graph.Edge>();

            return new Graph.HyperGraph(nodes, edges);
        }


        public int[] AQILabelInference(List<double>[] inAQI)
        {
            int[] aqiLabel = new int[this.NumDays * this.GridNum];
            for (int i = 0; i < this.NumDays * this.GridNum; i++)
            {
                if (inAQI[i][0] != Double.MaxValue) //there exists a label
                {
                    int aqi = Convert.ToInt32(inAQI[i][0]);
                    if (aqi <= 50)
                        aqiLabel[i] = 1; //Good
                    else if (aqi > 50 && aqi <= 100)
                        aqiLabel[i] = 2; //Moderate
                    else if (aqi > 100 && aqi <= 150)
                        aqiLabel[i] = 3; //Unhealthy for sensitive groups
                    else if (aqi > 150 && aqi <= 200)
                        aqiLabel[i] = 4; //Unhealthy
                    else if (aqi > 200 && aqi <= 300)
                        aqiLabel[i] = 5; //Very Unhealthy
                    else
                        aqiLabel[i] = 6; //Hazardous
                }
                else
                {
                    aqiLabel[i] = -1;
                }

            }
            return aqiLabel;
        }

        public void SetNumDays(string MeterologyDirectory)
        {
            DirectoryInfo dir = new DirectoryInfo(MeterologyDirectory);
            DirectoryInfo[] dirs = dir.GetDirectories();
            foreach(var d in dirs)
            {
                this.Days.Add(d.Name);
            }

        }

        public List<double>[] ReadFeature(string filename)
        {
            StreamReader sr = new StreamReader(filename);
            List<double>[] feature = new List<double>[this.GridNum];
            int indexRegion = 0;
            while (sr.Peek() > 0)
            {
                string[] items = sr.ReadLine().Split(new char[] { ' ' });
                feature[indexRegion] = new List<double>();

                if (items.Count() > 1)
                {
                    foreach (var item in items)
                    {
                        feature[indexRegion].Add(Convert.ToDouble(item));
                    }
                }
                else
                {
                    feature[indexRegion].Add(double.MaxValue);
                }
                indexRegion++;
            }
            sr.Close();
            return feature;
        }

        public Graph.HyperGraph BuildEdgesOfGraph(Graph.HyperGraph inGraph, int k1, int k2, double[] sigma, int radius, string outputGraphFilename)
        {
            //start to include the edges 

            for (int i = 0; i < inGraph.NumModality; i++)
            {
                List<Graph.SingleNode> mNodes = inGraph.Nodes[i];

                if (mNodes.Count() > 0)
                {
                    //there are two types of edges
                    //first type: intra-link type = 0
                    int nodeCount = mNodes.Count();
                    int nodeDim = mNodes[0].Data.Count();
                    double[,] inputData = new double[nodeCount, nodeDim];
                    for (int n = 0; n < nodeCount; n++)
                    {
                        if (mNodes[n].Data.Count() > 1)
                        {
                            for (int d = 0; d < nodeDim; d++)
                            {
                                inputData[n, d] = mNodes[n].Data[d];
                            }
                        }
                    }
                    int[] tags = new int[nodeCount];
                    for (int t = 0; t < nodeCount; t++)
                    {
                        tags[t] = t;
                    }
                    //build the k-d tree
                    alglib.kdtree kdt;
                    alglib.kdtreebuildtagged(inputData, tags, nodeDim, 0, 2, out kdt);
                    //find k-nn
                    int countNode = 0;
                    Dictionary<Preprocess.Pair<int, int>, double> mutual = new Dictionary<Preprocess.Pair<int, int>, double>();
                    foreach (var node in inGraph.Nodes[i])
                    {
                        Console.WriteLine(i + "th modality: " + countNode++ + "th node");
                        if (node.Data.Count() == 1)
                            break;
                        double[] x = new double[nodeDim];
                        for (int j = 0; j < nodeDim; j++)
                        {
                            x[j] = node.Data[j];
                        }
                        int k = alglib.kdtreequeryknn(kdt, x, k1);
                        double[,] result = new double[0, 0];
                        alglib.kdtreequeryresultsx(kdt, ref result);
                        double[] distances = new double[0];
                        alglib.kdtreequeryresultsdistances(kdt, ref distances);
                        int[] index = new int[0];
                        alglib.kdtreequeryresultstags(kdt, ref index);
                        List<int> foundBuffer = new List<int>();
                        List<double> css = new List<double>();
                        List<double> euDist = new List<double>();
                        for (int j = 1; j < k; j++) // j =0 is the data itself do not link itself
                        {
                            //List<double> tmp = new List<double>();
                            //for (int l = 0; l < nodeDim; l++)
                            //{
                            //    tmp.Add(result[j, l]);
                            //}

                            Graph.SingleNode neighbor = new Graph.SingleNode();
                            neighbor = mNodes[index[j]];
                            //foreach (var n in inGraph.Nodes[i])
                            //{
                            //    bool equal = true;
                            //    for (int l = 0; l < nodeDim; l++)
                            //    {
                            //        if (n.Data[l] != result[j, l])
                            //        {
                            //            equal = false;
                            //            break;
                            //        }
                            //    }
                            //    if (equal && n.NodeID != node.NodeID && !foundBuffer.Contains(n.NodeID))
                            //    {
                            //        neighbor = n;
                            //        foundBuffer.Add(n.NodeID);
                            //        break;
                            //    }
                            //}

                            double weight = Function.CalculateGaussianSimilarity(distances[j], i, sigma);
                            euDist.Add(weight);
                            //css.Add(Function.CalculateGaussianSimilarity(Function.GetCosineSimilarity(node.Data, neighbor.Data), i, sigma));
                            css.Add(Function.GetCosineSimilarity(node.Data, neighbor.Data));
                            int nodeID1 = node.NodeID;
                            int nodeID2 = neighbor.NodeID;

                            bool success = false;
                            if (mutual.ContainsKey(new Preprocess.Pair<int, int>(nodeID1, nodeID2)))
                            {
                                mutual[new Preprocess.Pair<int, int>(nodeID1, nodeID2)] = -1;
                            }
                            else if (mutual.ContainsKey(new Preprocess.Pair<int, int>(nodeID2, nodeID1)))
                            {
                                mutual[new Preprocess.Pair<int, int>(nodeID2, nodeID1)] = -1;
                            }
                            else
                            {
                                success = inGraph.SetEdge(node, neighbor, weight, 0);
                            }

                            if (success)
                            {
                                mutual.Add(new Preprocess.Pair<int, int>(nodeID1, nodeID2), weight);
                            }
                            //else
                            //{
                            //    if (mutual.ContainsKey(new Preprocess.Pair<int, int>(nodeID1, nodeID2)))
                            //    {
                            //        mutual[new Preprocess.Pair<int, int>(nodeID1, nodeID2)] = -1;
                            //    }
                            //    else
                            //    {
                            //        mutual[new Preprocess.Pair<int, int>(nodeID2, nodeID1)] = -1;
                            //    }
                            //}
                        }
                    }


                    List<Graph.Edge> delEdges = new List<Graph.Edge>();
                    //mutual knn
                    foreach (var m in mutual)
                    {
                        //if(m.Value == 1)
                        //{
                        //    inGraph.DelEdge(m.Key.Value1, m.Key.Value2);
                        //}
                        if (m.Value > 0)
                        {
                            delEdges.Add(new Graph.Edge(m.Key.Value1, m.Key.Value2, m.Value, 0));
                        }
                    }
                    inGraph.Edges = inGraph.Edges.Except(delEdges).ToList();

                    //second type: inter-link type = 1
                    for (int j = i + 1; j < inGraph.NumModality; j++)
                    {
                        List<Graph.SingleNode> nNodes = inGraph.Nodes[j];
                        //inGraph.PrintGraph(@"D:\Result\Temporary\Graph\graph.txt");
                        foreach (var mNode in mNodes)
                        {
                            int countNN = 0;
                            foreach (var nNode in nNodes)
                            {
                                //constraint of day
                                bool dayConstriant = mNode.Day == nNode.Day || mNode.Day == -1 || nNode.Day == -1;
                                if (radius > 0)
                                {
                                    if (this.Grid.IsAdjacent(mNode.GridIndex, nNode.GridIndex, radius) && dayConstriant) //DO NOT LINK DIFFERENT DAYS OF A REGION
                                    {
                                        inGraph.SetEdge(mNode, nNode, radius == 0 ? 1 : 1 / radius, 1);
                                        if (countNN >= k2)
                                            break;
                                        countNN++;
                                    }
                                }
                                else
                                {
                                    if (mNode.GridIndex == nNode.GridIndex && dayConstriant)
                                    {
                                        inGraph.SetEdge(mNode, nNode, 1, 1);
                                        break;
                                    }
                                }

                            }
                        }
                    }
                }
            }
            inGraph.PrintGraph(outputGraphFilename);
            return inGraph;
        }
    }
}
