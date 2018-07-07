using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace MMMST.CoupledDL
{
    public class Graph
    {
        /// <summary>
        /// A single node which is pure 
        /// </summary>
        /// <typeparam name="T"></typeparam>
        public class SingleNode 
        {
            public int NodeID;
            public List<double> Data;
            public int Label; // Label = 1-6 according to the aqi condition. Label = -1 means that there is no label for this node
            public int Modality; // Modality = 1-4 according to which dataset the node comes from.
            public int GridIndex; // indicate the node is in which grid
            public int Day;
            public int ClusterID;

            public SingleNode() { }

            public SingleNode(int nid, List<double> v, int l, int m, int gi, int day)
            {
                this.NodeID = nid;
                this.Data = v;
                this.Label = l;
                this.Modality = m;
                this.GridIndex = gi;
                this.Day = day;
                this.ClusterID = NodeID;
            }

            public override string ToString()
            {
                string result = this.Modality.ToString();
                result += ": " + this.NodeID;
                foreach(var d in this.Data)
                {
                    result = result + " " + d;
                }
                result = result + " " + this.GridIndex.ToString() + " " + this.Label.ToString() + " " + this.Day.ToString();
                return result;
            }
        }

        /// <summary>
        /// A composite node which contains dozens of single nodes. It can also be regarded as a cluster.
        /// </summary>
        /// <typeparam name="T"></typeparam>
        public class CompositeNode
        {
            public List<int> SingleNodes;
            public List<int> Labels;
            public List<int> Modalities;
            public int NumNodes;
            public int MaxLabel;
            public int MaxLabelCount;
            public int LabelCount;
            public int UnLabelCount;
            public Dictionary<int, int> NumModalities;

            public CompositeNode(List<int> singles, List<int> labels, List<int> modalities)
            {
                this.SingleNodes = new List<int>();
                foreach(var s in singles)
                {
                    this.SingleNodes.Add(s);
                }
                this.NumNodes = this.SingleNodes.Count();

                var distinctLabels = labels.Distinct().ToList();
                List<int> countDistinct = new List<int>();
                foreach(var dl in distinctLabels)
                {
                    int count = 0;
                    foreach(var l in labels)
                    {
                        if (l == dl)
                        {
                            count++;
                        }
                    }
                    countDistinct.Add(count);
                }
                this.Labels = new List<int>();
                this.Labels.AddRange(labels);

                this.LabelCount = 0;
                this.UnLabelCount = 0;
                int numDistinct = countDistinct.Count();
                for (int i = 0; i < NumNodes; i++)
                {
                    if (this.Labels[i] == -1)
                        this.UnLabelCount++;
                    else
                        this.LabelCount++;
                }

                int indexUnlabel = distinctLabels.IndexOf(-1); //remove the unlabeled one.
                if (indexUnlabel == -1)
                {
                    this.MaxLabelCount = countDistinct.Max();
                    this.MaxLabel = distinctLabels[countDistinct.IndexOf(this.MaxLabelCount)];
                }
                else
                {
                    distinctLabels.RemoveAt(indexUnlabel);
                    countDistinct.RemoveAt(indexUnlabel);

                    if (distinctLabels.Count > 0)
                    {
                        this.MaxLabelCount = countDistinct.Max();
                        this.MaxLabel = distinctLabels[countDistinct.IndexOf(this.MaxLabelCount)];
                    }
                    else
                    {
                        this.MaxLabelCount = 0;
                        this.MaxLabel = -1;
                    }
                }

                this.Modalities = new List<int>();
                foreach(var m in modalities)
                {
                    this.Modalities.Add(m);
                }

                this.NumModalities = new Dictionary<int, int>();
                foreach(var m in modalities)
                {
                    if(!this.NumModalities.ContainsKey(m))
                    {
                        this.NumModalities.Add(m, 1);
                    }
                    else
                    {
                        this.NumModalities[m]++;
                    }
                }

            }

            public override string ToString()
            {
                string cnode = "";
                foreach(var sn in this.SingleNodes)
                {
                    cnode = cnode + sn + " ";
                }
                cnode += "\t";
                foreach(var l in this.Labels)
                {
                    cnode = cnode + l + " ";
                }
                return cnode;
            }
        }

        public interface IGraph
        {
            //get the number of vertexes in the graph.
            int GetNumOfVertex();

            //get the number of edges in the graph.
            int GetNumOfEdge();

            //set the edge between two composite nodes (containing the situation where a composite node contains only a single node)
            bool SetEdge(SingleNode v1, SingleNode v2, double weight, int type); 

            //delete the edge between two composite nodes (containing the situation where a composite node contains only one singel node)
            void DelEdge(SingleNode v1, SingleNode v2);

            //see if there exist any edges between two composite nodes
            bool IsEdge(SingleNode v1, SingleNode v2);
        }

        public class ClusteringEdge : Edge, IComparable
        {
            public double Gain;

            public ClusteringEdge(int sid, int eid, double w, int type, double gain)
            {
                this.StartID = sid;
                this.EndID = eid;
                this.Weight = w;
                this.Type = type;
                this.Gain = gain;
            }

            public int CompareTo(object obj)
            {
                ClusteringEdge edge = obj as ClusteringEdge;
                if (this.Gain > edge.Gain)
                    return 1;
                else if (this.Gain == edge.Gain)
                    return 0;
                else
                    return -1;
            }

            public override string ToString()
            {
                return String.Format(this.StartID + " " + this.EndID + " " + this.Weight + " " + this.Type + " " + this.Gain);
            }
        }

        public class Edge : IEquatable<Edge>
        {
            public int StartID;
            public int EndID;
            public double Weight;
            public int Type; // Type = 0 intra link and Type = 1 inter link

            public Edge() { }

            public Edge(int sid, int eid, double w, int type)
            {
                this.StartID = sid;
                this.EndID = eid;
                this.Weight = w;
                this.Type = type;
            }

            public override string ToString()
            {
                return String.Format(this.StartID + " " + this.EndID + " " + this.Weight + " " + this.Type);
            }

            public override bool Equals(object obj)
            {
                if (obj == null) return false;
                Edge objAsEdge = obj as Edge;
                if (objAsEdge == null) return false;
                else return Equals(objAsEdge);
            }
            public override int GetHashCode()
            {
                return this.StartID ^ this.EndID ^ this.Weight.GetHashCode() ^ this.Type.GetHashCode();
            }
            public bool Equals(Edge other)
            {
                if (other == null) return false;
                return (this.StartID == other.StartID && this.EndID == other.EndID && this.Weight == other.Weight && this.Type == other.Type);
            }
        }


        public class HyperGraph : IGraph
        {
            public List<SingleNode>[] Nodes;
            public List<Edge> Edges;
            public Dictionary<Preprocess.Pair<int, int>, Edge> EdgeNodes;  
            public Dictionary<int, int> UniqAQI;
            public int NumModality;
            public int NumNodes;
            public int NumEdges;
            public int NumLabelNodes;
            public int NumUnLabelNodes;
            public Dictionary<int, int> NumEachModality;

            public HyperGraph() { }

            public HyperGraph(List<SingleNode>[] nodes, List<Edge> edges)
            {
                this.NumModality = nodes.Count();
                this.Nodes = new List<SingleNode>[this.NumModality];
                this.NumNodes = 0;
                this.NumLabelNodes = 0;
                this.NumUnLabelNodes = 0;

                this.NumEachModality = new Dictionary<int, int>();

                List<int> labels = new List<int>();
                this.UniqAQI = new Dictionary<int, int>();
                for (int i = 0; i < this.NumModality; i++)
                {
                    this.Nodes[i] = new List<SingleNode>();
                    foreach(var node in nodes[i])
                    {
                        this.Nodes[i].Add(node);
                        if (!this.NumEachModality.ContainsKey(i + 1))
                        {
                            this.NumEachModality.Add(i + 1, 1);
                        }
                        else
                        {
                            this.NumEachModality[i + 1] += 1;
                        }

                        if (node.Label == -1)
                            this.NumUnLabelNodes++;
                        else
                        {
                            this.NumLabelNodes++;
                            labels.Add(node.Label);
                            if (!this.UniqAQI.ContainsKey(node.Label))
                            {
                                this.UniqAQI.Add(node.Label, 1);
                            }
                            else
                            {
                                this.UniqAQI[node.Label]++;
                            }

                        }
                    }
                    this.NumNodes += nodes[i].Count();
                }
                this.Edges = new List<Edge>();
                this.EdgeNodes = new Dictionary<Preprocess.Pair<int, int>, Edge>();
                foreach(var e in edges)
                {
                    this.Edges.Add(e);
                    this.EdgeNodes.Add(new Preprocess.Pair<int, int>(e.StartID, e.EndID), e);
                }
                this.NumEdges = edges.Count();
            }

            public SingleNode FindNode(int nodeID)
            {
                //foreach (var node in this.Nodes)
                //{
                //    foreach (var n in node)
                //    {
                //        if (n.NodeID == nodeID)
                //            return n;
                //    }
                //}
                int a = 0;
                int b = 0;
                for (int i = 0; i < this.NumModality; i++)
                {
                    b += this.Nodes[i].Count();
                    if (nodeID >= a && nodeID < b)
                    {
                        return this.Nodes[i][nodeID - a];
                    }
                    a = b;
                }

                return null;
            }

            public int GetNumOfVertex()
            {
                return this.NumNodes;
            }

            public int GetNumOfEdge()
            {
                return this.NumEdges;
            }

            public int GetIndexVertex(SingleNode v)
            {
                return v.NodeID;
            }

            public int GetIndexEdge(Edge e)
            {
                return this.Edges.IndexOf(e);
            }

            public bool IsNode(SingleNode v1)
            {
                if (v1.NodeID >= 0 && v1.NodeID < this.NumNodes)
                    return true;
                else
                    return false;
            }

            public bool IsEdge(SingleNode v1, SingleNode v2)
            {
                //foreach (var e in this.Edges)
                //{
                //    if ((e.StartID == v1.NodeID && e.EndID == v2.NodeID) || (e.EndID == v1.NodeID && e.StartID == v2.NodeID))
                //    {
                //        return true;
                //    }
                //}
                
                if (this.EdgeNodes.ContainsKey(new Preprocess.Pair<int, int>(v1.NodeID, v2.NodeID)) || this.EdgeNodes.ContainsKey(new Preprocess.Pair<int, int>(v2.NodeID, v1.NodeID)))
                    return true;
                return false;
            }

            public bool IsEdge(SingleNode v1, SingleNode v2, out Edge edge)
            {
                //foreach (var e in this.Edges)
                //{
                //    if ((e.StartID == v1.NodeID && e.EndID == v2.NodeID) || (e.EndID == v1.NodeID && e.StartID == v2.NodeID))
                //    {
                //        edge = new Edge(e.StartID, e.EndID, e.Weight, e.Type);
                //        return true;
                //    }
                //}
                if (this.EdgeNodes.ContainsKey(new Preprocess.Pair<int, int>(v1.NodeID, v2.NodeID)))
                {
                    Edge e = this.EdgeNodes[new Preprocess.Pair<int, int>(v1.NodeID, v2.NodeID)];
                    edge = new Edge(e.StartID, e.EndID, e.Weight, e.Type);
                    return true;
                }
                else if (this.EdgeNodes.ContainsKey(new Preprocess.Pair<int, int>(v2.NodeID, v1.NodeID)))
                {
                    Edge e = this.EdgeNodes[new Preprocess.Pair<int, int>(v2.NodeID, v1.NodeID)];
                    edge = new Edge(e.StartID, e.EndID, e.Weight, e.Type);
                    return true;
                }
                
                edge = null;
                return false;
            }

            public bool IsEdge(int v1, int v2, out Edge edge)
            {
                //foreach (var e in this.Edges)
                //{
                //    if ((e.StartID == v1 && e.EndID == v2) || (e.EndID == v1 && e.StartID == v2))
                //    {
                //        edge = new Edge(e.StartID, e.EndID, e.Weight, e.Type);
                //        return true;
                //    }
                //}
                if (this.EdgeNodes.ContainsKey(new Preprocess.Pair<int, int>(v1, v2)))
                {
                    Edge e = this.EdgeNodes[new Preprocess.Pair<int, int>(v1, v2)];
                    edge = new Edge(e.StartID, e.EndID, e.Weight, e.Type);
                    return true;
                }
                else if (this.EdgeNodes.ContainsKey(new Preprocess.Pair<int, int>(v2, v1)))
                {
                    Edge e = this.EdgeNodes[new Preprocess.Pair<int, int>(v2, v1)];
                    edge = new Edge(e.StartID, e.EndID, e.Weight, e.Type);
                    return true;
                }

                edge = null;
                return false;
            }

            public bool SetEdge(SingleNode v1, SingleNode v2, double weight, int type)
            {
                if(!IsNode(v1) || !IsNode(v2))
                {
                    Console.WriteLine("Node is not belonging to the graph.");
                    return false;
                }
                if (IsEdge(v1, v2))
                {
                    Console.WriteLine("The edge exists already.");
                    return false;
                }
                this.Edges.Add(new Edge(v1.NodeID, v2.NodeID, weight, type));
                this.EdgeNodes.Add(new Preprocess.Pair<int, int>(v1.NodeID, v2.NodeID), new Edge(v1.NodeID, v2.NodeID, weight, type));
                this.NumEdges++;
                return true;
            }

            public void DelEdge(SingleNode v1, SingleNode v2)
            {
                if (!IsNode(v1) || !IsNode(v2))
                {
                    Console.WriteLine("Node is not belonging to the graph.");
                }
                Edge edge; 
                if (IsEdge(v1, v2, out edge))
                {
                    this.Edges.Remove(edge);
                    this.NumEdges--;
                }
            }

            public void DelEdge(int v1, int v2)
            {
                int numNode = GetNumOfVertex();
                if (v1 < 0 || v2 < 0 || v1 >= numNode || v2 >= numNode)
                {
                    Console.WriteLine("Node is not belonging to the graph.");
                }
                Edge edge;
                if (IsEdge(v1, v2, out edge))
                {
                    this.Edges.Remove(edge);
                    this.NumEdges--;
                }
            }

            public void PrintGraph(string outFilename)
            {
                StreamWriter swGraph = new StreamWriter(outFilename);
                swGraph.WriteLine(this.NumModality);
                swGraph.WriteLine(this.GetNumOfVertex());
                for (int i = 0; i < this.NumModality; i++)
                {
                    List<Graph.SingleNode> mNodes = this.Nodes[i];
                    foreach(var node in mNodes)
                    {
                        swGraph.WriteLine(node.ToString());
                    }
                }
                swGraph.WriteLine(this.GetNumOfEdge());
                foreach(var edge in this.Edges)
                {
                    swGraph.WriteLine(edge.ToString());
                }
                swGraph.Close();
            }

            public static HyperGraph ReadGraph(string inFilename)
            {
                StreamReader srGraph = new StreamReader(inFilename);
                int modalityNum = Convert.ToInt32(srGraph.ReadLine());
                int nodeNum = Convert.ToInt32(srGraph.ReadLine());
                List<SingleNode>[] nodes = new List<SingleNode>[modalityNum];
                for (int i = 0; i < modalityNum; i++)
                {
                    nodes[i] = new List<SingleNode>();
                }
                List<Edge> edges = new List<Edge>();

                for (int i = 0; i < nodeNum; i++)
                {
                    string[] items = srGraph.ReadLine().Split(new char[] { ':', ' ' }, StringSplitOptions.RemoveEmptyEntries);
                    int modality = Convert.ToInt32(items[0]);
                    int index = Convert.ToInt32(items[1]);
                    int count = items.Count();
                    int day = Convert.ToInt32(items[count - 1]);
                    int label = Convert.ToInt32(items[count - 2]);
                    int gridIndex = Convert.ToInt32(items[count - 3]);
                    List<double> data = new List<double>();
                    for (int j = 2; j < count - 3; j++)
                    {
                        data.Add(Convert.ToDouble(items[j]));
                    }
                    nodes[modality - 1].Add(new SingleNode(index, data, label, modality, gridIndex, day));
                }

                int edgeNum = Convert.ToInt32(srGraph.ReadLine());
                while(srGraph.Peek()>0)
                {
                    string[] items = srGraph.ReadLine().Split(new char[] { ' ' });
                    int startID = Convert.ToInt32(items[0]);
                    int endID = Convert.ToInt32(items[1]);
                    double weight = Convert.ToDouble(items[2]);
                    int type = Convert.ToInt32(items[3]);
                    edges.Add(new Edge(startID, endID, weight, type));
                }
                srGraph.Close();

                HyperGraph graph = new HyperGraph(nodes, edges);
                return graph;
            }
        }

    }
}
