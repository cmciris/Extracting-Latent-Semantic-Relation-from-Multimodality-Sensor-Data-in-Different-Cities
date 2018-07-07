using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace MMMST.Preprocess
{
    public class LabelSubset
    {
        public static void RetrieveLabelSubset(string cityDataDirectory, string cityFeatureDirectory, int timeSlot, int normType)
        {
            string roadFilename = Path.Combine(cityDataDirectory, @"roadnetwork\Road_Network_2014.txt");
            string roadFeatureFilename = Path.Combine(cityFeatureDirectory, @"roadNetworkFeature_0.015.txt");
            string poiFeatureFilename = Path.Combine(cityFeatureDirectory, @"cityPOIFeature_0.015.txt");
            string aqiDirectory = Path.Combine(cityFeatureDirectory, @"aqiFeature\0.015\");
            string mtlDirectory = Path.Combine(cityFeatureDirectory, @"meterologyFeature\mode0\0.015\");

            RoadNetwork roadNetwork = IO.ReadRoadNetwork(roadFilename);
            Grid grid = new Grid(roadNetwork.MinPoint, roadNetwork.MaxPoint, 0.015);

            int numDays1 = new DirectoryInfo(aqiDirectory).GetDirectories().Count();
            int numDays2 = new DirectoryInfo(mtlDirectory).GetDirectories().Count();
            int numDays = Math.Min(numDays1, numDays2);

            CoupledDL.GraphInput graphInput = new CoupledDL.GraphInput(timeSlot, numDays, 3, roadFeatureFilename, poiFeatureFilename, mtlDirectory, "", aqiDirectory, grid);
            CoupledDL.Graph.HyperGraph graph = graphInput.ReadGraphFromFiles2(normType);

            int numFeature = graph.Nodes[2].Count();
            List<double>[] meterologyFeature = new List<double>[numFeature];
            List<double>[] roadFeature = new List<double>[numFeature];
            List<double>[] poiFeature = new List<double>[numFeature];
            List<double>[] aqiLabel = new List<double>[numFeature];
            for (int n = 0; n < numFeature; n++)
            {
                meterologyFeature[n] = new List<double>();
                roadFeature[n] = new List<double>();
                poiFeature[n] = new List<double>();
                aqiLabel[n] = new List<double>();
            }

            int i=0;
            foreach (var node in graph.Nodes[2]) //meterology
            {
                int roadIndex = -1;
                foreach(var road in graph.Nodes[0])
                {
                    if (node.GridIndex == road.GridIndex)
                    {
                        roadIndex = graph.Nodes[0].IndexOf(road);
                        break;
                    }
                }
                int poiIndex = -1;
                foreach (var poi in graph.Nodes[1])
                {
                    if(poi.GridIndex == node.GridIndex)
                    {
                        poiIndex = graph.Nodes[1].IndexOf(poi);
                    }
                }
                if (roadIndex >= 0 && poiIndex >= 0)
                {
                    meterologyFeature[i] = node.Data;
                    roadFeature[i] = graph.Nodes[0][roadIndex].Data;
                    poiFeature[i] = graph.Nodes[1][poiIndex].Data;
                    aqiLabel[i].Add(node.Label);
                    i++;
                }
            }

            string outputDir = cityFeatureDirectory + @"\labeledSubset4\" + timeSlot.ToString("D2");
            if (!Directory.Exists(outputDir))
                Directory.CreateDirectory(outputDir);
            string filename0 = Path.Combine(outputDir, "road.txt");
            string filename1 = Path.Combine(outputDir, "poi.txt");
            string filename2 = Path.Combine(outputDir, "meterology.txt");
            string filename3 = Path.Combine(outputDir, "label.txt");
            IO.WriteFeature(roadFeature, filename0);
            IO.WriteFeature(poiFeature, filename1);
            IO.WriteFeature(meterologyFeature, filename2);
            IO.WriteFeature(aqiLabel, filename3);
        }

        public static void RetrieveLabelSubset4(string cityDataDirectory, string cityFeatureDirectory, int timeSlot, int normType)
        {
            string roadFilename = Path.Combine(cityDataDirectory, @"roadnetwork\Road_Network_2014.txt");
            string roadFeatureFilename = Path.Combine(cityFeatureDirectory, @"roadNetworkFeature_0.015.txt");
            string poiFeatureFilename = Path.Combine(cityFeatureDirectory, @"cityPOIFeature_0.015.txt");
            string aqiDirectory = Path.Combine(cityFeatureDirectory, @"aqiFeature\0.015\");
            string mtlDirectory = Path.Combine(cityFeatureDirectory, @"meterologyFeature\mode0\0.015\");
            string mobilityDirectory = Path.Combine(cityFeatureDirectory, @"mobilityFeature\original\timeSlot24\0.015\");

            RoadNetwork roadNetwork = IO.ReadRoadNetwork(roadFilename);
            Grid grid = new Grid(roadNetwork.MinPoint, roadNetwork.MaxPoint, 0.015);

            //int numDays1 = new DirectoryInfo(aqiDirectory).GetDirectories().Count();
            //int numDays2 = new DirectoryInfo(mtlDirectory).GetDirectories().Count();
            //int numDays3 = new DirectoryInfo(mobilityDirectory).GetDirectories().Count();
            //int numDays = Math.Min(numDays1, Math.Min(numDays2, numDays3));

            CoupledDL.GraphInput graphInput = new CoupledDL.GraphInput(timeSlot, 93, 4, roadFeatureFilename, poiFeatureFilename, mtlDirectory, mobilityDirectory, aqiDirectory, grid);
            CoupledDL.Graph.HyperGraph graph = graphInput.ReadGraphFromFiles3(normType);

            int numFeature = graph.Nodes[2].Count();
            List<double>[] meterologyFeature = new List<double>[numFeature];
            List<double>[] mobilityFeature = new List<double>[numFeature];
            List<double>[] roadFeature = new List<double>[numFeature];
            List<double>[] poiFeature = new List<double>[numFeature];
            List<double>[] aqiLabel = new List<double>[numFeature];
            List<double>[] metGridIndex = new List<double>[numFeature];
            for (int n = 0; n < numFeature; n++)
            {
                meterologyFeature[n] = new List<double>();
                roadFeature[n] = new List<double>();
                poiFeature[n] = new List<double>();
                aqiLabel[n] = new List<double>();
                mobilityFeature[n] = new List<double>();
                metGridIndex[n] = new List<double>();
            }

            int count = 0;
            int i = 0;
            foreach (var node in graph.Nodes[2]) //meterology
            {
                int roadIndex = -1;
                foreach (var road in graph.Nodes[0])
                {
                    if (node.GridIndex == road.GridIndex)
                    {
                        roadIndex = graph.Nodes[0].IndexOf(road);
                        break;
                    }
                }
                int poiIndex = -1;
                foreach (var poi in graph.Nodes[1])
                {
                    if (poi.GridIndex == node.GridIndex)
                    {
                        poiIndex = graph.Nodes[1].IndexOf(poi);
                        break;
                    }
                }
                if (roadIndex >= 0 && poiIndex >= 0)
                {
                    meterologyFeature[i] = node.Data;
                    mobilityFeature[i] = graph.Nodes[3][count].Data;
                    roadFeature[i] = graph.Nodes[0][roadIndex].Data;
                    poiFeature[i] = graph.Nodes[1][poiIndex].Data;
                    aqiLabel[i].Add(node.Label);
                    metGridIndex[i].Add(node.GridIndex);
                    i++;
                }
                count++;
                if (count != i)
                    Console.WriteLine();
            }

            string outputDir = cityFeatureDirectory + @"\labeledSubset4\" + timeSlot.ToString("D2");
            if (!Directory.Exists(outputDir))
                Directory.CreateDirectory(outputDir);
            string filename0 = Path.Combine(outputDir, "road.txt");
            string filename1 = Path.Combine(outputDir, "poi.txt");
            string filename2 = Path.Combine(outputDir, "meterology.txt");
            string filename3 = Path.Combine(outputDir, "label.txt");
            string filename4 = Path.Combine(outputDir, "mobility.txt");
            string filename5 = Path.Combine(outputDir, "metGridIndex.txt");
            IO.WriteFeature(roadFeature, filename0);
            IO.WriteFeature(poiFeature, filename1);
            IO.WriteFeature(meterologyFeature, filename2);
            IO.WriteFeature(mobilityFeature, filename4);
            IO.WriteFeature(aqiLabel, filename3);
            IO.WriteFeature(metGridIndex, filename5);
        }

    }
}
