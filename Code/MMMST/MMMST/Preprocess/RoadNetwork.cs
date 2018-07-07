using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MMMST.Utils;

namespace MMMST.Preprocess
{
    public class RoadNetwork
    {
        public int NodeNum;
        public int EdgeNum;
        public List<Node> Nodes;
        public List<Edge> Edges;
        public LatLong MaxPoint;
        public LatLong MinPoint;


        public RoadNetwork(int nn, int en, List<Node> nodes, List<Edge> edges)
        {
            this.NodeNum = nn;
            this.EdgeNum = en;
            this.Nodes = new List<Node>();
            this.Edges = new List<Edge>();

            foreach(var nd in nodes)
            {
                this.Nodes.Add(nd);
            }
            foreach(var ed in edges)
            {
                this.Edges.Add(ed);
            }

            List<LatLong> latlongs = new List<LatLong>();
            foreach(var node in nodes)
            {
                latlongs.Add(node.LatLon);
            }
            Geo.CalculateBoundingBox(latlongs, out this.MaxPoint, out this.MinPoint);
        }

        public double[] CountHighwayLength()
        {
            double[] length = new double[2];
            foreach(var e in this.Edges)
            {
                if (e.RoadClass == 0)
                    length[0] += e.Length;
                else
                    length[1] += e.Length;
            }
            return length;
        }

        /// <summary>
        /// This function generates the features of the road network one region by one region.
        /// </summary>
        /// <param name="grid"></param>The grid partitions a specific space into different regions.
        /// <param name="roadFeatureFilename"></param>The parameter specifies where to output the feature matrix.
        public void RoadFeatures(Grid grid, string roadFeatureFilename)
        {
            int gridNum = grid.GridNumLat * grid.GridNumLon;
            List<double>[,] roadAttributes = new List<double>[gridNum, this.Edges[0].AttributeNum];// the array indicates the region. the first level of list encodes the different attributes. the second level of list encodes all the roads's an attribute of a region
            for (int i = 0; i < gridNum; i++)
            {
                for (int j = 0; j < this.Edges[0].AttributeNum; j++)
                {
                    roadAttributes[i, j] = new List<double>();
                }
            }

            //pre-save all the attributes of edges according to regions for future use.
            foreach (var edge in this.Edges)
            {
                List<int> indexes = new List<int>();
                foreach (var point in edge.Points)
                {
                    int pIndex = grid.RetrieveIndex(point);
                    if (!indexes.Contains(pIndex))
                    {
                        indexes.Add(pIndex);           //tranverse all the regions that the edge passes
                    }
                }

                int countAttribute = 0;
                roadAttributes[grid.RetrieveIndex(Nodes[edge.StartID].LatLon), countAttribute].Add(edge.StartID);
                roadAttributes[grid.RetrieveIndex(Nodes[edge.EndID].LatLon), countAttribute].Add(edge.EndID);
                countAttribute++;
                
                foreach (var id in indexes)
                {
                    countAttribute = 1;
                    //attribute 1: road segments
                    roadAttributes[id, countAttribute].Add(1.0);
                    countAttribute++;
                    //attribute 2: road length
                    roadAttributes[id, countAttribute].Add(edge.Length);
                    countAttribute++;
                    //attribute 3: road class
                    roadAttributes[id, countAttribute].Add(edge.RoadClass);
                    countAttribute++;
                    //attribute 4: road direction
                    roadAttributes[id, countAttribute].Add(edge.Direction);
                    countAttribute++;
                    //attribute 5: road formway
                    roadAttributes[id, countAttribute].Add(edge.FormWay);
                    countAttribute++;
                    //attribute 6: urban 
                    roadAttributes[id, countAttribute].Add(edge.Urban);
                    countAttribute++;
                    //attribute 7: max lanes
                    roadAttributes[id, countAttribute].Add(edge.MaxLanes);
                    countAttribute++;
                    //attribute 8: max speed
                    roadAttributes[id, countAttribute].Add(edge.MaxSpeed);
                    countAttribute++;
                    //attribute 9: road level
                    roadAttributes[id, countAttribute].Add(edge.Level);
                    countAttribute++;
                }
            }

            //start to design some road features 
            List<double>[] roadFeature = new List<double>[gridNum];
            for (int i = 0; i < gridNum; i++)
            {
                roadFeature[i] = new List<double>();
                if (roadAttributes[i, 1].Sum() > 0.0) //verify that there exist roads in the region
                {
                    //feature 1: the number of road segments
                    roadFeature[i].Add(roadAttributes[i, 1].Sum());
                    //feature 2: the number of intersections
                    double result = (from node in roadAttributes[i, 0]
                                     group node by node into newgroup
                                     where newgroup.Count() > 1 //the situation that two road segments share a node is also regarded as an intersection
                                     select newgroup.Key).Count();
                    roadFeature[i].Add(result);
                    //feature 3: total length of a specific road level (10 features)
                    double[] length = new double[Edge.MAXLEVEL + 1];
                    int edgeCount = roadAttributes[i, 9].Count();
                    for (int j = 0; j < edgeCount; j++)
                    {
                        length[(int)roadAttributes[i, 9][j]] += roadAttributes[i, 2][j];
                    }
                    foreach (var l in length)
                    {
                        roadFeature[i].Add(l);
                    }
                    //feature 4: total length
                    roadFeature[i].Add(roadAttributes[i, 2].Sum());
                    //feature 5: the majority level of the region
                    var majorLevel = (from level in roadAttributes[i, 9]
                                      group level by level into newgroup
                                      orderby newgroup.Count() descending
                                      select newgroup.Key).First();
                    roadFeature[i].Add(majorLevel);
                    //feature 6: the average number of maximum lanes 
                    roadFeature[i].Add(roadAttributes[i, 7].Average());
                    //feature 7: the average maximum speed 
                    roadFeature[i].Add(roadAttributes[i, 8].Average());
                    //feature 8: intersections with road > 2
                    result = (from node in roadAttributes[i, 0]
                                     group node by node into newgroup
                                     where newgroup.Count() > 2 //the situation that two road segments share a node is also regarded as an intersection
                                     select newgroup.Key).Count();
                    roadFeature[i].Add(result);
                    //feature 9: total length of a specific road class (12 features)
                    length = new double[12];
                    for (int j = 0; j < edgeCount; j++)
                    {
                        length[(int)roadAttributes[i, 3][j]] += roadAttributes[i, 2][j];
                    }
                    foreach (var l in length)
                    {
                        roadFeature[i].Add(l);
                    }
                    //feature 10: total length of highways
                    roadFeature[i].Add(length[0]);
                    //feature 11: total length of non-highways 
                    roadFeature[i].Add(length.Sum() - length[0]);
                    //feature 12: the majority maximum speed of the region
                    var majorSpeed = (from level in roadAttributes[i, 8]
                                      group level by level into newgroup
                                      orderby newgroup.Count() descending
                                      select newgroup.Key).First();
                    roadFeature[i].Add(majorSpeed);
                    //feature 13: the majority number of maximum lanes 
                    var majorLanes = (from level in roadAttributes[i, 7]
                                      group level by level into newgroup
                                      orderby newgroup.Count() descending
                                      select newgroup.Key).First();
                    roadFeature[i].Add(majorLanes);
                    //feature 14: the majority number of maximum road class
                    var majorClass = (from level in roadAttributes[i, 3]
                                      group level by level into newgroup
                                      orderby newgroup.Count() descending
                                      select newgroup.Key).First();
                    roadFeature[i].Add(majorClass);
                }
                else
                {
                    roadFeature[i].Add(Double.MaxValue);
                }
            }

            //output the road features into the file
            IO.WriteFeature(roadFeature, roadFeatureFilename);
        }
    }

    public class Node
    {
        public int NodeID;
        public LatLong LatLon;

        public Node(int nodeID, LatLong ll)
        {
            this.NodeID = nodeID;
            this.LatLon = new LatLong(ll.Latitude, ll.Longitude);
        }
    }

    public class Edge
    {
        public int EdgeID;
        public int StartID;
        public int EndID;
        public double Length;
        public int RoadClass;
        public int Direction;
        public int FormWay;
        public int Urban;
        public int MaxLanes;
        public int MaxSpeed;
        public int Level;
        public static int MAXLEVEL = 9;
        public int AttributeNum;
        public List<LatLong> Points;

        public Edge()
        {

        }

        public Edge(int edgeID, int sid, int eid, double l, int rc, int d, int fw, int ur, int ml, int ms, int level, int an, List<LatLong> points)
        {
            this.EdgeID = edgeID;
            this.StartID = sid;
            this.EndID = eid;
            this.Length = l;
            this.RoadClass = rc;
            this.Direction = d;
            this.FormWay = fw;
            this.Urban = ur;
            this.MaxLanes = ml;
            this.MaxSpeed = ms;
            this.Level = level;
            this.AttributeNum = an;
            this.Points = new List<LatLong>();
            
            foreach(var p in points)
            {
                Points.Add(p);
            }
        }
    }
}
