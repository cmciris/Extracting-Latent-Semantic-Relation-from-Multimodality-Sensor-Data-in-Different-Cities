using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MMMST.Preprocess
{
    public class Grid
    {
        public LatLong Min;
        public LatLong Max;
        public double GridSize;
        public int GridNumLon;
        public int GridNumLat;

        public Grid(LatLong min, LatLong max, double gs)
        {
            this.Min = new LatLong(min.Latitude, min.Longitude);
            this.Max = new LatLong(max.Latitude, max.Longitude);
            this.GridSize = gs;

            GridNumLon = (int)(Math.Ceiling((this.Max.Longitude - this.Min.Longitude) / GridSize));
            GridNumLat = (int)(Math.Ceiling((this.Max.Latitude - this.Min.Latitude) / GridSize));
        }

        /// <summary>
        /// given a target location point, retrive the index of it in this grid
        /// </summary>
        /// <param name="target"></param>specifies the location to be retrieved
        /// <returns></returns>
        public int RetrieveIndex(LatLong target)
        {
            int index = 0;
            int indexLon = (int)Math.Floor((target.Longitude - this.Min.Longitude) / this.GridSize);
            int indexLat = (int)Math.Floor((target.Latitude - this.Min.Latitude) / this.GridSize);
            
            index = indexLat * GridNumLon + indexLon;
            return index;
        }

        /// <summary>
        /// given an index, retrieve the max and min point of its boundary in this grid
        /// </summary>
        /// <param name="index"></param>specifies the index to be retrieved
        /// <param name="maxPoint"></param>output, the maxpoint of the boundary
        /// <param name="minPoint"></param>output, the minpoint of the boundary
        public void RetrieveBoundary(int index, out LatLong maxPoint, out LatLong minPoint)
        {
            int indexLon = index % GridNumLon;
            int indexLat = index / GridNumLon;

            double minLat = this.Min.Latitude + this.GridSize * indexLat;
            double minLon = this.Min.Longitude + this.GridSize * indexLon;
            double maxLat = this.Min.Latitude + this.GridSize * (indexLat + 1);
            double maxLon = this.Min.Longitude + this.GridSize * (indexLon + 1);

            maxPoint = new LatLong(maxLat, maxLon);
            minPoint = new LatLong(minLat, minLon);
        }

        /// <summary>
        /// judge wether two indexed regions are adjacent or not given the specified radius.
        /// </summary>
        /// <param name="index1"></param>the first region's index
        /// <param name="index2"></param>the second region's index
        /// <param name="radius"></param>specifies the radius of adjacency
        /// <returns></returns>
        public bool IsAdjacent(int index1, int index2, int radius)
        {
            int minIndex = 0, maxIndex = 0;
            if (index1 < this.GridNumLon)
                minIndex = index1 < radius ? index1 : index1 - radius;
            else
                minIndex = index1 - this.GridNumLon * radius - radius;
            if (GridNumLat * GridNumLon - index1 <= this.GridNumLon)
                maxIndex = GridNumLat * GridNumLon - index1 <= radius ? index2 : index1 + radius;
            else
                maxIndex = index1 + this.GridNumLon * radius + radius;


            LatLong minPoint1, maxPoint1, minPoint2, maxPoint2, minPoint, maxPoint;
            RetrieveBoundary(minIndex, out maxPoint1, out minPoint1);
            RetrieveBoundary(maxIndex, out maxPoint2, out minPoint2);
            RetrieveBoundary(index2, out maxPoint, out minPoint);

            if (minPoint.Latitude >= minPoint1.Latitude && minPoint.Longitude >= minPoint1.Longitude && maxPoint.Latitude <= maxPoint2.Latitude && maxPoint.Longitude <= maxPoint2.Longitude)
                return true;
            else
                return false;
        }
    }

    public class LatLong : IEquatable<LatLong>
    {
        public double Latitude;
        public double Longitude;

        public LatLong(double lat, double lon)
        {
            this.Latitude = lat;
            this.Longitude = lon;
        }

        public bool Equals(LatLong other)
        {
            if (other == null)
                return false;

            if (this.Latitude == other.Latitude && this.Longitude == other.Longitude)
                return true;
            else
                return false;
        }

        public override bool Equals(Object obj)
        {
            if (obj == null)
                return false;

            LatLong latlong = obj as LatLong;
            if (latlong == null)
                return false;
            else
                return Equals(latlong);
        }

        public override int GetHashCode()
        {
            return this.Latitude.GetHashCode() ^ this.Longitude.GetHashCode();
        }

        public override string ToString()
        {
            return string.Format(this.Latitude + ", " + this.Longitude);
        }
    }
}
