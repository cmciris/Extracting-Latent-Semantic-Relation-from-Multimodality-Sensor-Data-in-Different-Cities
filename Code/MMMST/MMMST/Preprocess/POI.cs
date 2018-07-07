using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MMMST.Utils;

namespace MMMST.Preprocess
{
    public class CityPOI
    {
        public int POINum;
        public List<POI> POIs;

        public CityPOI(int pn, List<POI> pois)
        {
            this.POINum = pn;
            List<LatLong> latlongs = new List<LatLong>();
            this.POIs = new List<POI>();
            foreach(var poi in pois)
            {
                latlongs.Add(poi.Location);
                this.POIs.Add(poi);
            }
        }

        /// <summary>
        /// generate the features of POI one region by one region given the grid 
        /// </summary>
        /// <param name="grid"></param>specifies the way that partitions the space into different regions
        /// <param name="cityPOIFeatureFilename"></param>specifies the path where the feature matrix to be saved
        public void CityPOIFeature(Grid grid, string cityPOIFeatureFilename)
        {
            int gridNum = grid.GridNumLat * grid.GridNumLon;

            //figure out which pois are contained in per region
            List<POI>[] regionPOIs = new List<POI>[gridNum];
            for (int i = 0; i < gridNum; i++)
            {
                regionPOIs[i] = new List<POI>();
            }
            foreach (var poi in POIs)
            {
                if (grid.RetrieveIndex(poi.Location) < gridNum && grid.RetrieveIndex(poi.Location) >= 0)
                {
                    regionPOIs[grid.RetrieveIndex(poi.Location)].Add(poi);
                }

            }

            //start to generate features of each region
            List<double>[] cityPOIFeature = new List<double>[gridNum];
            for (int i = 0; i < gridNum; i++)
            {
                cityPOIFeature[i] = new List<double>();
                if (regionPOIs[i].Count() > 0)
                {
                    //feature 1: the total number of pois
                    cityPOIFeature[i].Add(regionPOIs[i].Count());
                    //feature 2: the density of pois 
                    LatLong minPoint, maxPoint;
                    grid.RetrieveBoundary(i, out maxPoint, out minPoint);
                    Grid smallGrid = new Grid(minPoint, maxPoint, grid.GridSize / 10);
                    List<int> occupiedIndex = new List<int>();
                    //feature 3: the categorical count
                    int[] catCount = new int[POI.CATEGORYNUM];
                    //feature 4: average rank
                    List<double> rankScore = new List<double>();

                    foreach (var regionPOI in regionPOIs[i])
                    {
                        int index = smallGrid.RetrieveIndex(regionPOI.Location);
                        if (!occupiedIndex.Contains(index))
                        {
                            occupiedIndex.Add(index);
                        }
                        catCount[regionPOI.Category]++;
                        rankScore.Add(regionPOI.Rank);
                    }
                    cityPOIFeature[i].Add(((double)occupiedIndex.Count()) / (smallGrid.GridNumLat * smallGrid.GridNumLon));
                    foreach (var cat in catCount)
                    {
                        cityPOIFeature[i].Add(cat);
                    }
                    cityPOIFeature[i].Add(rankScore.Average());
                }
                else
                {
                    cityPOIFeature[i].Add(Double.MaxValue);
                }
            }

            //output the feature matrix of the city-level pois
            IO.WriteFeature(cityPOIFeature, cityPOIFeatureFilename);
        }
    }



    public class POI
    {
        public string ID;
        public LatLong Location;
        public int Category; //01. 汽车服务 02. 汽车销售 03. 汽车维修 04. 摩托车服务 05. 餐饮服务 06. 购物服务 07. 生活服务 08. 体育休闲服务 09. 医疗保健服务 10. 住宿服务
                             //11. 风景名胜 12. 商务住宅 13. 政府机构及社会团体 14. 科教文化服务 15. 交通设施服务 16. 金融保险服务 17. 公司企业 18. 道路附属设施 19. 地名地址信息 20. 公共设施
        public int Rank;
        public static int CATEGORYNUM = 20;

        public POI(string id, LatLong loc, int cat, int r)
        {
            this.ID = id;
            this.Category = cat;
            this.Rank = r;
            this.Location = new LatLong(loc.Latitude, loc.Longitude);
        }


    }
}
