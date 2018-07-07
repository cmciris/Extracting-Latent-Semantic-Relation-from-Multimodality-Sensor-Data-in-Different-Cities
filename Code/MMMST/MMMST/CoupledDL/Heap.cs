using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MMMST.CoupledDL
{
    public class Heap
    {
        public class SubmodularHeap : MaxHeap<Graph.ClusteringEdge>
        {
            public SubmodularHeap(List<Graph.ClusteringEdge> data)
            {
                this.Data = new List<Graph.ClusteringEdge>();
                foreach (var d in data)
                {
                    Insert(d);
                } 
            }

            public void UpdateSubmodularHeap(List<Graph.CompositeNode> clusters, Graph.HyperGraph inGraph, double[] loop, int numLabelNodes, int numUnLabelNodes, double lambda, double gamma, double mu)
            {
                while (UpdateHeapValue(clusters, inGraph, loop, 0, numLabelNodes, numUnLabelNodes, lambda, gamma, mu) == 0) //only update the root node acoording to the diminishing return property of submodular functions
                {
                    if (this.IsEmpty())
                    {
                        break;
                    }

                    if (this.Data[0].Gain == 0)
                        this.ExtractMax();
                    else
                        this.MaxHeapify(0);
                    //this.CheckMaxHeap();
                }
            }

            public int UpdateHeapValue(List<Graph.CompositeNode> clusters, Graph.HyperGraph inGraph, double[] loop, int i, int numLabelNodes, int numUnLabelNodes, double lambda, double gamma, double mu)
            {
                if(this.IsEmpty())
                {
                    return 1;
                }
                double hGain = 0, pGain = 0, lGain = 0, mGain = 0; mu = 0 ;
                //store the old gain
                double oldGain = this.Data[i].Gain;
                //find the index of the clusters
                int clusterID1 = Function.FindIndex(inGraph, this.Data[i].StartID);
                int clusterID2 = Function.FindIndex(inGraph, this.Data[i].EndID);

                //if the edge forms a cycle, we make the gain zero.
                //later we will remove the the zero-gain edges from the heap
                if (clusterID1 == clusterID2)
                    this.Data[i].Gain = 0;
                else
                {
                    //recompute the entropy rate gain 
                    hGain = Function.CalculateHGain(this.Data[i].Weight, loop[this.Data[i].StartID] - this.Data[i].Weight, loop[this.Data[i].EndID] - this.Data[i].Weight);
                    //recopmute the pure function gain 
                    pGain = Function.CalculatePGain(clusters, numLabelNodes, clusterID1, clusterID2);
                    //recopmute the balance function gain
                    lGain = Function.CalculateLGain(clusters, numLabelNodes, numUnLabelNodes, clusterID1, clusterID2);
                    //recompute the modality function gain
                    mGain = Function.CalculateMGain(clusters, inGraph.NumEachModality, clusterID1, clusterID2);

                    this.Data[i].Gain = hGain + lambda * pGain + gamma * lGain + mu * mGain;
                }

                if (oldGain == this.Data[i].Gain)
                    return 1;
                return 0;
            }
        }


        public class MaxHeap<T> where T : IComparable
        {
            public List<T> Data = new List<T>();

            public MaxHeap() { }

            public MaxHeap(List<T> data)
            {
                this.Data = new List<T>();
                foreach(var d in data)
                {
                    Insert(d);
                }
            }

            public bool IsEmpty()
            {
                return Data.Count == 0;
            }

            public void CheckMaxHeap()
            {
                for (int i = 0; i < Data.Count; i++)
                {
                    int l = 2 * (i + 1) - 1;
                    int r = 2 * (i + 1) - 1 + 1;

                    //check if the left child is no greater than the parent
                    if (l < Data.Count)
                    {
                        if (Data[i].CompareTo(Data[l]) < 0)
                            Console.WriteLine("Error in heap!");
                    }
                    if (r < Data.Count)
                    {
                        if (Data[i].CompareTo(Data[r]) < 0)
                            Console.WriteLine("Error in heap!");
                    }
                }
            }

            public void Insert(T o)
            {
                Data.Add(o);

                int i = Data.Count - 1;
                while (i > 0)
                {
                    int j = (i + 1) / 2 - 1;

                    // Check if the invariant holds for the element in data[i]
                    T v = Data[j];
                    if (v.CompareTo(Data[i]) > 0 || v.CompareTo(Data[i]) == 0)
                    {
                        break;
                    }

                    // Swap the elements
                    T tmp = Data[i];
                    Data[i] = Data[j];
                    Data[j] = tmp;

                    i = j;
                }
            }

            public T ExtractMax()
            {
                if (Data.Count < 0)
                {
                    throw new ArgumentOutOfRangeException();
                }

                T max = Data[0];
                Data[0] = Data[Data.Count - 1];
                Data.RemoveAt(Data.Count - 1);
                this.MaxHeapify(0);
                return max;
            }

            public T Peek()
            {
                return Data[0];
            }

            public int Count
            {
                get { return Data.Count; }
            }

            public void MaxHeapify(int i)
            {
                int largest;
                int l = 2 * (i + 1) - 1;
                int r = 2 * (i + 1) - 1 + 1;

                if (l < Data.Count && (Data[l].CompareTo(Data[i]) > 0))
                {
                    largest = l;
                }
                else
                {
                    largest = i;
                }

                if (r < Data.Count && (Data[r].CompareTo(Data[largest]) > 0))
                {
                    largest = r;
                }

                if (largest != i)
                {
                    T tmp = Data[i];
                    Data[i] = Data[largest];
                    Data[largest] = tmp;
                    this.MaxHeapify(largest);
                }
            }

            public void PrintHeap(string filename)
            {
                System.IO.StreamWriter swHeap = new System.IO.StreamWriter(filename);

                foreach (var d in Data)
                {
                    swHeap.WriteLine(d);
                }

                swHeap.Close();
            }
        }

    }
}
