using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MMMST.Preprocess;

namespace MMMST.Utils
{
    public class Geo
    {
        /// <summary>
        /// Initiate the bounding box of a specific dataset by calculating the max and min point
        /// </summary>
        /// <param name="latlongs"></param>the list of latitudes and longitudes to infer the max and the min
        /// <param name="maxPoint"></param>the point with maximum in both the latitude and longitute
        /// <param name="minPoint"></param>the point with minimum in both the latitude and longitude
        public static void CalculateBoundingBox(List<LatLong> latlongs, out LatLong maxPoint, out LatLong minPoint)
        {
            double minLat = Double.MaxValue;
            double minLon = Double.MaxValue;
            double maxLat = 0.0;
            double maxLon = 0.0;

            foreach (var single in latlongs)
            {
                if (single.Latitude < minLat)
                {
                    minLat = single.Latitude;
                }
                if (single.Latitude > maxLat)
                {
                    maxLat = single.Latitude;
                }
                if (single.Longitude < minLon)
                {
                    minLon = single.Longitude;
                }
                if (single.Longitude > maxLon)
                {
                    maxLon = single.Longitude;
                }
            }
            maxPoint = new LatLong(maxLat, maxLon);
            minPoint = new LatLong(minLat, minLon);
        }

        /// <summary>
        /// Find the top k nearest neighbors to the target point
        /// </summary>
        /// <param name="target"></param>specifies the target point
        /// <param name="neighbours"></param>specifies all the neighbours pool to be possibly retrieved 
        /// <returns></returns>
        public static List<LatLong> FindKNN(LatLong target, List<LatLong> neighbours, int k)
        {
            Dictionary<LatLong, double> distances = new Dictionary<LatLong, double>();
            foreach (var n in neighbours)
            {
                double distance = GeoDistance(target, n);
                if (!distances.ContainsKey(n))
                {
                    distances.Add(n, distance);
                }
            }

            var distancesSort = (from d in distances
                                 orderby d.Value
                                 ascending
                                 select d.Key).ToList();
            return distancesSort.Take(k).ToList();
        }

        /// <summary>
        /// Calculate the geographical distance between point a and point b
        /// </summary>
        /// <param name="a"></param>geopoint a
        /// <param name="b"></param>geopoint b
        /// <returns></returns>
        public static double GeoDistance(LatLong x1, LatLong x2)
        {
            double radLat1 = rad(x1.Latitude);
            double radLat2 = rad(x2.Latitude);
            double a = radLat1 - radLat2;
            double b = rad(x1.Longitude) - rad(x2.Longitude);
            double s = 2 * Math.Asin(Math.Sqrt(Math.Pow(Math.Sin(a / 2), 2) +
            Math.Cos(radLat1) * Math.Cos(radLat2) * Math.Pow(Math.Sin(b / 2), 2)));
            s = s * EARTH_RADIUS;
            s = Math.Round(s * 10000) / 10000;
            return s;  
        }

        private const double EARTH_RADIUS = 6378.137;
        private static double rad(double d)
        {
            return d * Math.PI / 180.0;
        } 
    }
}
