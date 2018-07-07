using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MMMST.Utils;
using System.IO;

namespace MMMST.Preprocess
{
    public class ToTalMobility:IDisposable
    {
        public List<Node> Nodes;
        public Mobility[] Mobilities;
        public int NumTimeSlots;

        bool disposed = false;

        public void Dispose()
        {
            Dispose(true);
            GC.SuppressFinalize(this);
        }

        protected virtual void Dispose(bool disposing)
        {
            if (disposed)
                return;
            if (disposing)
            {
                Nodes = null;
                Mobilities = null;
            }

            disposed = true;
        }

        ~ToTalMobility()
        {
            Dispose(false);
        }

        public ToTalMobility() { }

        public ToTalMobility(List<Node> nodes, Mobility[] m)
        {
            this.Nodes = new List<Node>();
            foreach (var node in nodes)
            {
                this.Nodes.Add(node);
            }
            this.NumTimeSlots = m.Count();
            this.Mobilities = new Mobility[this.NumTimeSlots];
            for (int i = 0; i < this.NumTimeSlots; i++)
            {
                Mobilities[i] = m[i];
            }
        }

        public void TotalMobilityFeature(Grid grid, string mobilityFeatureDirectory)
        {
            int gridNum = grid.GridNumLat * grid.GridNumLon;

            DirectoryInfo dir = new DirectoryInfo(mobilityFeatureDirectory);

            for (int i = 0; i < this.NumTimeSlots; i++)
            {
                Mobility mobility = this.Mobilities[i];

                if (mobility != null)
                {
                    //figure out which edges are contained in each region 
                    List<MobilityEdge>[] regionEdges = new List<MobilityEdge>[gridNum];
                    for (int j = 0; j < gridNum; j++)
                    {
                        regionEdges[j] = new List<MobilityEdge>();
                    }
                    foreach (var mobilityEdge in mobility.Edges)
                    {
                        List<int> indexRegion = new List<int>();
                        foreach (var latlong in mobilityEdge.Points)
                        {
                            int index = grid.RetrieveIndex(latlong);
                            if (!indexRegion.Contains(index))
                            {
                                indexRegion.Add(index);
                            }
                        }

                        foreach (var index in indexRegion)
                        {
                            regionEdges[index].Add(mobilityEdge);
                        }
                    }

                    //figure out which pickdroppoints are contained in each region 
                    List<PickDropPoint>[] regionPoints = new List<PickDropPoint>[gridNum];
                    for (int j = 0; j < gridNum; j++)
                    {
                        regionPoints[j] = new List<PickDropPoint>();
                    }
                    foreach (var point in mobility.Points)
                    {
                        int index = grid.RetrieveIndex(point.Location);
                        if (index < gridNum && index >= 0)
                        {
                            regionPoints[index].Add(point);
                        }
                    }

                    //start to generate features
                    List<double>[] mobilityFeature = new List<double>[gridNum];
                    for (int j = 0; j < gridNum; j++)
                    {
                        mobilityFeature[j] = new List<double>();
                        //feature 1: expectation of the average speed regarding to road lengths in the region
                        double avgAvgSpeed = 0.0;
                        //feature 2: deviation of the average speed regarding to road lengths in the region
                        double stdAvgSpeed = 0.0;
                        //feature 3: expectation of the deviation speed regarding to road lengths in the region
                        double avgStdSpeed = 0.0;
                        //feature 4: deviation of the deviation speed regarding to road lengths in the region 
                        double stdStdSpeed = 0.0;
                        //feature 5: number of vechicles that tranverse the region
                        int numVehicles = 0;
                        //feature 6: number of road segments with average speed in (0,6)
                        int numSpeedSix = 0;
                        //feature 7: number of road segments with average speed in (6-12)
                        int numSpeedTwelve = 0;
                        //feature 8: number of road segments with average speed in (12-18)
                        int numSpeedEighteen = 0;
                        //feature 9: number of road segments with average speed in (18-)
                        int numSpeedTop = 0;
                        //feature 12: number of road segments with average speed in (0,2)
                        int numSpeed1= 0;
                        //feature 13: number of road segments with average speed in (2,4)
                        int numSpeed2 = 0;
                        //feature 14: number of road segments with average speed in (4,6)
                        int numSpeed3 = 0;
                        //feature 15: number of road segments with average speed in (6,8)
                        int numSpeed4 = 0;
                        //feature 16: number of road segments with average speed in (8,10)
                        int numSpeed5 = 0;
                        //feature 17: number of road segments with average speed in (10,12)
                        int numSpeed6 = 0;
                        //feature 18: number of road segments with average speed in (12,14)
                        int numSpeed7 = 0;
                        //feature 19: number of road segments with average speed in (14,16)
                        int numSpeed8 = 0;
                        //feature 20: number of road segments with average speed in (16,18)
                        int numSpeed9 = 0;
                        //feature 21: number of road segments with average speed in (18-)
                        int numSpeed10 = 0;

                        double totalLength = 0.0;
                        if (regionEdges[j].Count() > 0)
                        {
                            foreach (var edge in regionEdges[j])
                            {
                                if (edge.AverageSpeed < 6)
                                {
                                    numSpeedSix++;
                                }
                                else if (edge.AverageSpeed >= 6 && edge.AverageSpeed < 12)
                                {
                                    numSpeedTwelve++;
                                }
                                else if (edge.AverageSpeed >= 12 && edge.AverageSpeed < 18)
                                {
                                    numSpeedEighteen++;
                                }
                                else
                                {
                                    numSpeedTop++;
                                }

                                if (edge.AverageSpeed < 2)
                                {
                                    numSpeed1++;
                                }
                                else if (edge.AverageSpeed >= 2 && edge.AverageSpeed < 4)
                                {
                                    numSpeed2++;
                                }
                                else if (edge.AverageSpeed >= 4 && edge.AverageSpeed < 6)
                                {
                                    numSpeed3++;
                                }
                                else if (edge.AverageSpeed >= 6 && edge.AverageSpeed < 8)
                                {
                                    numSpeed4++;
                                }
                                else if (edge.AverageSpeed >= 8 && edge.AverageSpeed < 10)
                                {
                                    numSpeed5++;
                                }
                                else if (edge.AverageSpeed >= 10 && edge.AverageSpeed < 12)
                                {
                                    numSpeed6++;
                                }
                                else if (edge.AverageSpeed >= 12 && edge.AverageSpeed < 14)
                                {
                                    numSpeed7++;
                                }
                                else if (edge.AverageSpeed >= 14 && edge.AverageSpeed < 16)
                                {
                                    numSpeed8++;
                                }
                                else if (edge.AverageSpeed >= 16 && edge.AverageSpeed < 18)
                                {
                                    numSpeed9++;
                                }
                                else
                                {
                                    numSpeed10++;
                                }

                                avgAvgSpeed += edge.AverageSpeed;// *edge.Length;
                                avgStdSpeed += edge.StdSpeed; //* edge.Length;
                                numVehicles += edge.NumVehicles;
                                totalLength += 1; // edge.Length;
                            }
                            avgAvgSpeed /= totalLength;
                            avgStdSpeed /= totalLength;
                            foreach (var edge in regionEdges[j])
                            {
                                stdAvgSpeed += Math.Pow((edge.AverageSpeed - avgAvgSpeed), 2);// *edge.Length;
                                stdStdSpeed += Math.Pow((edge.StdSpeed - avgStdSpeed), 2);// *edge.Length;
                            }
                            stdAvgSpeed = Math.Sqrt(stdAvgSpeed / totalLength);
                            stdStdSpeed = Math.Sqrt(stdStdSpeed / totalLength);
                        }
                        mobilityFeature[j].Add(avgAvgSpeed);
                        mobilityFeature[j].Add(stdAvgSpeed);
                        mobilityFeature[j].Add(avgStdSpeed);
                        mobilityFeature[j].Add(stdStdSpeed);
                        mobilityFeature[j].Add(numVehicles);
                        mobilityFeature[j].Add(numSpeedSix);
                        mobilityFeature[j].Add(numSpeedTwelve);
                        mobilityFeature[j].Add(numSpeedEighteen);
                        mobilityFeature[j].Add(numSpeedTop);


                        //feature 10: number of pick up points in the region: meaning number of people outgoing
                        int pickup = 0;
                        //feature 11: number of drop off points in the region: meaning number of people ingoing
                        int dropoff = 0;

                        foreach (var point in regionPoints[j])
                        {
                            if (point.Flag == 1)
                            {
                                pickup++;
                            }
                            else
                            {
                                dropoff++;
                            }
                        }
                        mobilityFeature[j].Add(pickup);
                        mobilityFeature[j].Add(dropoff);

                        //mobilityFeature[j].Add(numSpeed1);
                        //mobilityFeature[j].Add(numSpeed2);
                        //mobilityFeature[j].Add(numSpeed3);
                        //mobilityFeature[j].Add(numSpeed4);
                        //mobilityFeature[j].Add(numSpeed5);
                        //mobilityFeature[j].Add(numSpeed6);
                        //mobilityFeature[j].Add(numSpeed7);
                        //mobilityFeature[j].Add(numSpeed8);
                        //mobilityFeature[j].Add(numSpeed9);
                        //mobilityFeature[j].Add(numSpeed10);
                    }

                    //start to output
                    string dirName = Path.Combine(dir.FullName, mobility.TimeStamp.Year.ToString() + mobility.TimeStamp.Month.ToString("D2") + mobility.TimeStamp.Day.ToString("D2"));
                    if (!Directory.Exists(dirName))
                    {
                        Directory.CreateDirectory(dirName);
                    }

                    string path = Path.Combine(dir.FullName, Path.Combine(dirName, mobility.TimeStamp.TimeSlot.ToString("D2")));
                    IO.WriteFeature(mobilityFeature, path);

                }
            }
        }
    }

    public class Mobility
    {
        public List<MobilityEdge> Edges;
        public List<PickDropPoint> Points;
        public Time TimeStamp;


        public Mobility() { }

        public Mobility(List<MobilityEdge> edges, List<PickDropPoint> points, Time ts)
        {
            this.Edges = new List<MobilityEdge>();
            foreach(var edge in edges)
            {
                this.Edges.Add(edge);
            }
            this.Points = new List<PickDropPoint>();
            foreach(var point in points)
            {
                this.Points.Add(point);
            }
            this.TimeStamp = new Time(ts.Year, ts.Month, ts.Day, ts.TimeSlot);
        }

    }

    public class PickDropPoint
    {
        public LatLong Location;
        public int Flag; //flag = 1 represents "pick up" and flag = 2 represents "drop off"

        public PickDropPoint(LatLong l, int f)
        {
            this.Location = new LatLong(l.Latitude, l.Longitude);
            this.Flag = f;
        }
    }

    public class MobilityEdge:Edge
    {
        public double AverageSpeed;
        public int NumVehicles;
        public double StdSpeed;

        public MobilityEdge() { }

        public MobilityEdge(int edgeID, int sid, int eid, double l, int rc, int d, int fw, int ur, int ml, int ms, int level, int an, List<LatLong> points, double avgs, int nv, double stds)
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

            foreach (var p in points)
            {
                Points.Add(p);
            }

            this.AverageSpeed = avgs;
            this.NumVehicles = nv;
            this.StdSpeed = stds;
        }
    }

    public class Time:IComparable
    {
        public int Year;
        public int Month;
        public int Day;
        public int TimeSlot; //there are total 48 timeslots with each time slot as 30 minutes
        public static int NUMSLOTPERDAY = 24;

        public Time(int year, int month, int day, int ts)
        {
            this.Year = year;
            this.Month = month;
            this.Day = day;
            this.TimeSlot = ts;
        }

        public override int GetHashCode()
        {
            return this.Year ^ this.Month ^ this.Day ^ this.TimeSlot;
        }

        public override bool Equals(object obj)
        {
            if (obj == null)
                return false;
            Time time = obj as Time;
            if (time == null)
                return false;
            else
                return this.Year == time.Year && this.Month == time.Month && this.Day == time.Day && this.TimeSlot == time.TimeSlot;
        }

        public int CompareTo(object obj)
        {
            Time time = obj as Time;
            if (this.Year > time.Year)
            {
                return 1;
            }
            else if (this.Year < time.Year)
            {
                return -1;
            }
            else
            {
                if (this.Month > time.Month)
                {
                    return 1;
                }
                else if (this.Month < time.Month)
                {
                    return -1;
                }
                else
                {
                    if (this.Day > time.Day)
                    {
                        return 1;
                    }
                    else if (this.Day < time.Day)
                    {
                        return -1;
                    }
                    else
                    {
                        if (this.TimeSlot > time.TimeSlot)
                        {
                            return 1;
                        }
                        else if (this.TimeSlot < time.TimeSlot)
                        {
                            return -1;
                        }
                        else
                        {
                            return 0;
                        }
                    }
                }
            }
        }
    }
}
