using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace MMMST.CoupledDL
{
    public class LazyGreedy
    {
        public static void OutputClusterResult(List<Graph.CompositeNode> clusters, Graph.HyperGraph inGraph, string outputClusterDirectory, int threshold)
        {
            int numClusters = clusters.Count();
            int[] sortID = new int[numClusters];
            int id = 0;
            for (int l = 1; l <= 6; l++)
            {
                for (int i = 0; i < numClusters; i++)
                {
                    if (clusters[i].MaxLabel == l)
                        sortID[i] = id++;
                }
            }

            int numModality = inGraph.NumModality;
            for (int i = 0; i < numModality; i++)
            {
                List<double>[] results = new List<double>[inGraph.NumNodes];
                for (int j = 0; j < inGraph.NumNodes; j++)
                {
                    results[j] = new List<double>();
                }
                int count = 0;

                for (int l = 1; l <= 6; l++)
                {
                    foreach (var node in inGraph.Nodes)
                    {
                        foreach (var n in node)
                        {
                            if (n.Modality == i + 1)// && n.Label == l)
                            {
                                int clusterID = Function.FindIndex(inGraph, n.NodeID);
                                if (clusters[clusterID].SingleNodes.Count() > threshold && clusters[clusterID].MaxLabel == l)
                                {
                                    results[count].AddRange(n.Data);
                                    results[count].Add(sortID[clusterID]); //the last dimension is the cluster ID
                                    results[count].Add(clusters[clusterID].MaxLabel);
                                    //results[count].Add(n.Label);
                                    if (clusters[clusterID].LabelCount == 0)
                                        results[count].Add(0.0);
                                    else
                                        results[count].Add((double)clusters[clusterID].MaxLabelCount / clusters[clusterID].LabelCount);
                                    count++;
                                }
                            }
                        }
                    }
                }
                Preprocess.IO.WriteFeature(results, System.IO.Path.Combine(outputClusterDirectory, i.ToString("D2") + ".txt"));
            }
        }

        public static void InferDictionary(string clusterDirectory, int numClusters, string outputDictDirectory)
        {
            DirectoryInfo dir = new DirectoryInfo(clusterDirectory);
            FileInfo[] files = dir.GetFiles();

            int countFile = 0;
            foreach (var f in files)
            {
                List<List<double>>[] clusters = new List<List<double>>[numClusters];
                List<double>[] dict = new List<double>[numClusters];
                for (int j = 0; j < numClusters; j++)
                {
                    clusters[j] = new List<List<double>>();
                    dict[j] = new List<double>();
                }
                StreamReader srCluster = new StreamReader(f.FullName);
                while (srCluster.Peek() > 0)
                {
                    string[] items = srCluster.ReadLine().Split(new char[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
                    int count = items.Count();
                    if (count > 1)
                    {
                        List<double> tmp = new List<double>();
                        for (int i = 0; i < count; i++)
                        {
                            tmp.Add(Convert.ToDouble(items[i]));
                        }
                        clusters[Convert.ToInt32(items[count - 3])].Add(tmp);
                    }
                }
                srCluster.Close();

                
                for (int j = 0; j < numClusters; j++)
                {
                    if (clusters[j].Count() > 1)
                    {
                        int dim = clusters[j][0].Count();
                        for (int m = 0; m < dim; m++)
                        {
                            double avg = 0.0;
                            foreach (var d in clusters[j])
                            {
                                avg += d[m];
                            }
                            avg /= clusters[j].Count();
                            dict[j].Add(avg);
                        }
                    }
                }
                Preprocess.IO.WriteFeature(dict, Path.Combine(outputDictDirectory, f.Name));
                countFile++;
                if (countFile > 3)
                    break;
            }
        }

        public static List<Graph.CompositeNode> Clustering(Graph.HyperGraph inGraph, double lambda, double gamma, double mu, int numClusters, string outputDirectory, out double entropy, out double purity, out double balance, out double diversity, out double objective)
        {
            #region initialize
            //when begin, each cluster contains a single node
            List<Graph.CompositeNode> clusters = new List<Graph.CompositeNode>();
            foreach(var nodei in inGraph.Nodes)
            {
                foreach(var node in nodei)
                {
                    List<int> tmp = new List<int>();
                    List<int> tmp2 = new List<int>();
                    List<int> tmp3 = new List<int>();

                    tmp.Add(node.NodeID);
                    tmp2.Add(node.Label);
                    tmp3.Add(node.Modality);
                    clusters.Add(new Graph.CompositeNode(tmp, tmp2, tmp3));
                }
            }

            //for debugging
            //CoupledDL.LazyGreedy.OutputClusterResult(clusters, inGraph, @"D:\Result\Temporary\Cluster\");

            //calculate the self loop weight and the total of weight of each vertex
            double[] loop = Function.CalculateLoopWeight(inGraph);
            double totalWeight = Function.CalculateTotalWeight(loop);
            Function.NormalizeWeight(ref loop, ref inGraph, totalWeight, 1);
            //calculate the initial objective gain and decide the trade-off parameter
            double[] hGain = new double[inGraph.NumEdges];
            double[] pGain = new double[inGraph.NumEdges];
            double[] lGain = new double[inGraph.NumEdges];
            double[] mGain = new double[inGraph.NumEdges];

            double maxHGain = 0, maxPGain = 0, maxLGain = 0, maxMGain = 0;

            for (int i = 0; i < inGraph.NumEdges; i++)
            {
                Console.WriteLine(i + "th edge's gain is initialized.");
                hGain[i] = Function.CalculateHGain(inGraph.Edges[i].Weight, loop[inGraph.Edges[i].StartID] - inGraph.Edges[i].Weight, loop[inGraph.Edges[i].EndID] - inGraph.Edges[i].Weight);
                int clusterIndex1 = Function.FindIndex(inGraph, inGraph.Edges[i].StartID);
                int clusterIndex2 = Function.FindIndex(inGraph, inGraph.Edges[i].EndID);
                if (clusterIndex1 != clusterIndex2)
                {
                    pGain[i] = Function.CalculatePGain(clusters, inGraph.NumLabelNodes, clusterIndex1, clusterIndex2);
                    lGain[i] = Function.CalculateLGain(clusters, inGraph.NumLabelNodes, inGraph.NumUnLabelNodes, clusterIndex1, clusterIndex2);
                    mGain[i] = Function.CalculateMGain(clusters, inGraph.NumEachModality, clusterIndex1, clusterIndex2);
                }
                if (hGain[i] > maxHGain)
                    maxHGain = hGain[i];
                if (pGain[i] > maxPGain)
                    maxPGain = pGain[i];
                if (lGain[i] > maxLGain)
                    maxLGain = lGain[i];
                if (mGain[i] > maxMGain)
                    maxMGain = mGain[i];
            }

            double adjustedLambda = lambda * maxHGain / maxPGain;
            double adjustedGamma = gamma * maxHGain / maxLGain;
            double adjustedMu = mu * maxHGain / maxMGain;

            List<Graph.ClusteringEdge> edges = new List<Graph.ClusteringEdge>();
            for (int i = 0; i < inGraph.NumEdges; i++)
            {
                double gain = hGain[i] + adjustedLambda * pGain[i] + adjustedGamma * lGain[i];
                edges.Add(new Graph.ClusteringEdge(inGraph.Edges[i].StartID, inGraph.Edges[i].EndID, inGraph.Edges[i].Weight, inGraph.Edges[i].Type, gain));
            }
            hGain = null; pGain = null; lGain = null; mGain = null;
            #endregion

            #region build the heap
            Heap.SubmodularHeap heap = new Heap.SubmodularHeap(edges);
            heap.CheckMaxHeap(); // for debugging
            #endregion

            #region greedily add the edge and perform the clustering
            Graph.ClusteringEdge bestEdge;
            int countCluster = inGraph.NumNodes;
            List<double> addedEdges = new List<double>();
            while (countCluster > numClusters)
            {
                if (heap.IsEmpty())
                {
                    Console.WriteLine("Heap is empty!");
                    objective = Function.CalculateObjective(clusters, inGraph, addedEdges, loop, adjustedLambda, adjustedGamma, adjustedMu, out entropy, out purity, out balance, out diversity);
                    loop = null;
                    return clusters;
                }

                Console.WriteLine(countCluster + "->" + numClusters + "clusters remained");
                //find the best edge to add
                bestEdge = heap.ExtractMax();
                addedEdges.Add(bestEdge.Weight);

                //merge the clusters which the best edge connects
                int clusterID1 = Function.FindIndex(inGraph, bestEdge.StartID);
                int clusterID2 = Function.FindIndex(inGraph, bestEdge.EndID);
                if (clusterID1 != clusterID2)
                {
                    int remove = Function.MergeTwoClusters(ref clusters, clusterID1, clusterID2);
                    int stay = clusterID1 == remove ? clusterID2 : clusterID1;
                    foreach(var c in clusters[stay].SingleNodes)
                    {
                        inGraph.FindNode(c).ClusterID = stay;
                    }
                    countCluster--;
                    for (int c = remove; c < countCluster; c++)
                    {
                        foreach(var n in clusters[c].SingleNodes)
                        {
                            inGraph.FindNode(n).ClusterID--;
                        }
                    }

                    loop[bestEdge.StartID] -= bestEdge.Weight;
                    loop[bestEdge.EndID] -= bestEdge.Weight;
                }

                heap.UpdateSubmodularHeap(clusters, inGraph, loop, inGraph.NumLabelNodes, inGraph.NumUnLabelNodes, adjustedLambda, adjustedGamma, adjustedMu);
            }

            objective = Function.CalculateObjective(clusters, inGraph, addedEdges, loop, adjustedLambda, adjustedGamma, adjustedMu, out entropy, out purity, out balance, out diversity);

            //string clusterFilename = Path.Combine(outputDirectory, "cluster.txt");
            //string heapFilename = Path.Combine(outputDirectory, "heap.txt");
            //string paraFilename = Path.Combine(outputDirectory, "parameter.txt");

            //Function.WriteClusters(clusters, clusterFilename);
            //Function.WriteParameters(loop, totalWeight, adjustedLambda, adjustedGamma, paraFilename);
            //heap.PrintHeap(heapFilename);

            loop = null;
            return clusters;
            #endregion
        }  
    }
}
