using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MMMST.CoupledDL
{
    class Function
    {
        public static double GetCosineSimilarity(List<double> V1, List<double> V2)
        {
            double sim = 0.0d;
            int N = 0;
            N = ((V2.Count < V1.Count) ? V2.Count : V1.Count);
            double dot = 0.0d;
            double mag1 = 0.0d;
            double mag2 = 0.0d;
            for (int n = 0; n < N; n++)
            {
                dot += V1[n] * V2[n];
                mag1 += Math.Pow(V1[n], 2);
                mag2 += Math.Pow(V2[n], 2);
            }

            return dot / (Math.Sqrt(mag1) * Math.Sqrt(mag2));
        }

        public static void WriteParameters(double[] loop, double totalWeight, double lambda, double gamma, string filename)
        {
            System.IO.StreamWriter swPara = new System.IO.StreamWriter(filename);
            swPara.WriteLine(totalWeight);
            swPara.WriteLine(lambda);
            swPara.WriteLine(gamma);
            foreach(var l in loop)
            {
                swPara.WriteLine(l);
            }
            swPara.Close();
        }

        public static void ReadParameters(string filename, out double[] loop, out double totalWeight, out double lambda, out double gamma)
        {
            System.IO.StreamReader srPara = new System.IO.StreamReader(filename);
            totalWeight = Convert.ToDouble(srPara.ReadLine());
            lambda = Convert.ToDouble(srPara.ReadLine());
            gamma = Convert.ToDouble(srPara.ReadLine());
            int count = 0;
            List<double> loops = new List<double>();
            while (srPara.Peek() > 0)
            {
                loops.Add(Convert.ToDouble(srPara.ReadLine()));
                count++;
            }
            loop = new double[count];
            for (int i = 0; i < count; i++)
            {
                loop[i] = loops[i];
            }
            srPara.Close();
        }

        public static void WriteClusters(List<Graph.CompositeNode> clusters, string filename)
        {
            System.IO.StreamWriter swClusters = new System.IO.StreamWriter(filename);
            foreach(var c in clusters)
            {
                swClusters.WriteLine(c.ToString());
            }
            swClusters.Close();
        }

        public static List<Graph.CompositeNode> ReadClusters(string filename)
        {
            System.IO.StreamReader srClusters = new System.IO.StreamReader(filename);
            List<Graph.CompositeNode> clusters = new List<Graph.CompositeNode>();
            while (srClusters.Peek() > 0)
            {
                string content = srClusters.ReadLine();
                string[] items = content.Split(new char[] { '\t' }, StringSplitOptions.RemoveEmptyEntries);
                string[] nodes = items[0].Split(new char[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
                string[] labels = items[1].Split(new char[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
                List<int> node = new List<int>();
                foreach(var n in nodes)
                {
                    node.Add(Convert.ToInt32(n));
                }
                List<int> label = new List<int>();
                foreach(var l in labels)
                {
                    label.Add(Convert.ToInt32(l));
                }
                List<int> modality = new List<int>();
                foreach(var m in modality)
                {
                    modality.Add(Convert.ToInt32(m));
                }
                clusters.Add(new Graph.CompositeNode(node, label, modality));
            }
            srClusters.Close();
            return clusters;
        }

        public static Heap.SubmodularHeap ReadSubmodularHeap(string filename)
        {
            System.IO.StreamReader srHeap = new System.IO.StreamReader(filename);
            List<Graph.ClusteringEdge> data = new List<Graph.ClusteringEdge>();
            while (srHeap.Peek() > 0)
            {
                string content = srHeap.ReadLine();
                string[] items = content.Split(new char[] { ' ' });
                Graph.ClusteringEdge ce = new Graph.ClusteringEdge(Convert.ToInt32(items[0]), Convert.ToInt32(items[1]), Convert.ToDouble(items[2]), Convert.ToInt32(items[3]), Convert.ToDouble(items[4]));
                data.Add(ce);
            }
            srHeap.Close();
            Heap.SubmodularHeap subHeap = new Heap.SubmodularHeap(data);

            return subHeap;
        }

        public static double CalculateObjective(List<Graph.CompositeNode> clusters, Graph.HyperGraph inGraph, List<double> addedEdges, double[] loop, double adjustedLambda, double adjustedGamma, double adjustedMu, out double entropy, out double purity, out double balance, out double diversity)
        {
            entropy = 0.0;
            //first part, calculate the entropy rate 
            foreach(var e in addedEdges)
            {
                double weight = e;
                entropy += -Xlogx(weight);
            }
            for (int i = 0; i < inGraph.NumNodes; i++)
            {
                double loopWeight = loop[i];
                entropy += -Xlogx(loopWeight); 
            }
            purity = 0.0;
            //second part, calculate the total purity
            int numClusters = clusters.Count();
            for (int i = 0; i < numClusters; i++)
            {
                purity += clusters[i].MaxLabelCount;
            }
            purity /= inGraph.NumLabelNodes;
            purity -= numClusters;
            //third part, calculate the balancing part
            balance = 0.0;
            for (int i = 0; i < numClusters; i++)
            {
                double p1 = (double)clusters[i].LabelCount / inGraph.NumLabelNodes;
                balance += -Xlogx(p1);
                double p2 = (double)clusters[i].UnLabelCount / inGraph.NumUnLabelNodes;
                balance += -Xlogx(p2);
            }
            balance -= 2 * numClusters;

            //fourth part, calculate the modality part
            diversity = 0.0; adjustedMu = 0.0;
            for (int m = 1; m <= inGraph.NumModality; m++)
            {
                double p = 0.0;
                for (int i = 0; i < numClusters; i++)
                {
                    if (clusters[i].NumModalities.ContainsKey(m))
                    {
                        double tmp = (double)clusters[i].NumModalities[m] / inGraph.NumEachModality[m];
                        p += Xlogx(tmp);
                    }
                }
                p *= (double)inGraph.NumEachModality[m] / inGraph.NumNodes;
                diversity += p;
            }
            diversity -= numClusters;

            return entropy + adjustedLambda * purity + adjustedGamma * balance + adjustedMu * diversity;
        }


        public static List<double> AverageList(List<double>[] lists)
        {
            List<double> result = new List<double>();
            int count = lists.Length;
            int dim = lists[0].Count;

            for (int i = 0; i < dim; i++)
            {
                double r = 0.0;
                for (int j = 0; j < count; j++)
                {
                    r += lists[j][i];
                }
                r /= count;
                result.Add(r);
            }
            return result;
        }

        public static int MergeTwoClusters(ref List<Graph.CompositeNode> clusters, int iid, int jid)
        {
            List<int> merged = new List<int>();
            List<int> mergedLabels = new List<int>();
            List<int> mergedModalities = new List<int>();

            merged.AddRange(clusters[iid].SingleNodes);
            merged.AddRange(clusters[jid].SingleNodes);
            mergedLabels.AddRange(clusters[iid].Labels);
            mergedLabels.AddRange(clusters[jid].Labels);
            mergedModalities.AddRange(clusters[iid].Modalities);
            mergedModalities.AddRange(clusters[jid].Modalities);
            if (iid <= jid)
            {
                clusters.RemoveAt(jid);
                clusters[iid] = new Graph.CompositeNode(merged, mergedLabels, mergedModalities);
                return jid;
            }
            else
            {
                clusters.RemoveAt(iid);
                clusters[jid] = new Graph.CompositeNode(merged, mergedLabels, mergedModalities);
                return iid;
            }
        }

        public static int FindIndex(Graph.HyperGraph inGraph, int nodeID)
        {
            //foreach(var c in clusters)
            //{
            //    foreach(var n in c.SingleNodes)
            //    {
            //        if (n == nodeID)
            //            return clusters.IndexOf(c);
            //    }
            //}
            //return -1;

            Graph.SingleNode v = inGraph.FindNode(nodeID);
            return v.ClusterID;
        }


        public static double Xlogx(double x)
        {
            if (x > 0)
                return x * Math.Log(x);
            else
                return 0;
        }

        public static double CalculateMGain(List<Graph.CompositeNode> clusters, Dictionary<int, int> numEachModality, int iid, int jid)
        {
            double mGain = 0;
            int numAll = 0;
            foreach(var numM in numEachModality.Keys)
            {
                if (clusters[iid].NumModalities.ContainsKey(numM) && clusters[jid].NumModalities.ContainsKey(numM))
                {
                    double mi = clusters[iid].NumModalities[numM] * 1.0 / numEachModality[numM];
                    double mj = clusters[jid].NumModalities[numM] * 1.0 / numEachModality[numM];
                    mGain += numEachModality[numM] * (-Xlogx(mi + mj) + Xlogx(mi) + Xlogx(mj));
                }
                numAll += numEachModality[numM];
            }
            mGain /= numAll;
            mGain += 1.0;
            return mGain;
        }

        public static double CalculateLGain(List<Graph.CompositeNode> clusters, int numLabelNodes, int numUnLabelNodes, int iid, int jid)
        {
            double li = clusters[iid].LabelCount * 1.0 / numLabelNodes;
            double lj = clusters[jid].LabelCount * 1.0 / numLabelNodes;
            double ui = clusters[iid].UnLabelCount * 1.0 / numUnLabelNodes;
            double uj = clusters[jid].UnLabelCount * 1.0 / numUnLabelNodes;


            double lGain = -Xlogx(li + lj) + Xlogx(li) + Xlogx(lj) + 1.0
                           - Xlogx(ui + uj) + Xlogx(ui) + Xlogx(uj) + 1.0;
            return lGain;
        }

        public static double CalculatePGain(List<Graph.CompositeNode> clusters, int numLabelNodes, int iid, int jid)
        {
            int maxCounti = clusters[iid].MaxLabelCount;
            int maxCountj = clusters[jid].MaxLabelCount;

            List<Graph.CompositeNode> tmpClusters = new List<Graph.CompositeNode>();
            tmpClusters.AddRange(clusters);
            MergeTwoClusters(ref tmpClusters, iid, jid);

            int index = iid <= jid ? iid : jid;
            int maxCountMerged = tmpClusters[index].MaxLabelCount;

            double pGain = (maxCountMerged - maxCounti - maxCountj) / ((double)numLabelNodes) + 1;
            return pGain;
        }

        public static double CalculateHGain(double wij, double li, double lj)
        {
            double hGain = Xlogx(wij + li) + Xlogx(wij + lj) - Xlogx(li) - Xlogx(lj) - 2 * Xlogx(wij); // /log(2.0) can be considered later
            return hGain;
        }

        public static void NormalizeWeight(ref double[] loop, ref Graph.HyperGraph inGraph, double totalWeight, int init)
        {
            if (init == 1)
            {
                for (int i = 0; i < inGraph.NumNodes; i++)
                {
                    loop[i] /= totalWeight;
                }
            }
            foreach(var e in inGraph.Edges)
            {
                e.Weight /= totalWeight;
            }
        }

        public static double CalculateTotalWeight(double[] loop)
        {
            return loop.Sum();
        }

        public static double[] CalculateLoopWeight(Graph.HyperGraph inGraph)
        {
            double[] loop = new double[inGraph.NumNodes];
            foreach(var e in inGraph.Edges)
            {
                loop[e.StartID] += e.Weight;
                loop[e.EndID] += e.Weight;
            }
            return loop;
        }

        public static double CalculateGaussianSimilarity(double euDist, int m, double[] sigmas)
        {
            double sigma = sigmas[m];
            double twoSigmaSquare = 2 * sigma * sigma;
            return Math.Exp(-euDist * euDist / twoSigmaSquare);
        }
    }
}
