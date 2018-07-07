using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace MMMST.Preprocess
{
    public class TotalAirQuality
    {
        public List<AirQuality>[] AirQualities;
        public int NumTimeSlots;
        public List<string> Days;
        public Dictionary<int, LatLong> AQIStations;

        public TotalAirQuality() { }

        public TotalAirQuality(List<AirQuality>[] airQualities, Dictionary<int, LatLong> aqiStations, List<string> days)
        {
            this.NumTimeSlots = airQualities.Count();
            this.AirQualities = new List<AirQuality>[this.NumTimeSlots]; 

            for (int i = 0; i < this.NumTimeSlots; i++)
            {
                AirQualities[i] = new List<AirQuality>();
                foreach(var aq in airQualities[i])
                {
                    this.AirQualities[i].Add(aq);
                }
            }
            this.AQIStations = new Dictionary<int, LatLong>();
            this.AQIStations = aqiStations;
            this.Days = new List<string>();
            foreach(var d in days)
            {
                this.Days.Add(d);
            }
        }

        public void TotalAirQualityFeature(Grid grid, string aqiFeatureDirectory)
        {
            int gridNum = grid.GridNumLat * grid.GridNumLon;
            int countDay = this.NumTimeSlots / Time.NUMSLOTPERDAY;

            for (int i = 0; i < this.NumTimeSlots; i++)
            {
                List<double>[] airQualityFeature = new List<double>[gridNum];
                List<AirQuality> airQualities = this.AirQualities[i];

                List<AirQuality>[] regionAQIs = new List<AirQuality>[gridNum];
                for (int j = 0; j < gridNum; j++)
                {
                    airQualityFeature[j] = new List<double>();
                    regionAQIs[j] = new List<AirQuality>();
                }
                foreach (var aq in airQualities)
                {
                    int regionIndex = grid.RetrieveIndex(this.AQIStations[aq.StationID]);
                    if (regionIndex < gridNum && regionIndex >= 0)
                        regionAQIs[regionIndex].Add(aq);
                }

                for (int j = 0; j < gridNum; j++)
                {
                    //label 1: aqi
                    double aqi = 0.0;
                    //label 2: pm25
                    double pm25 = 0.0;
                    //label 3: pm10
                    double pm10 = 0.0;
                    //label 4: no2
                    double no2 = 0.0;
                    //label 5: co
                    double co = 0.0;
                    //label 6: o3
                    double o3 = 0.0;
                    //label 7: so2
                    double so2 = 0.0;

                    int count = regionAQIs[j].Count();
                    if (count > 0)
                    {
                        foreach (var aq in regionAQIs[j])
                        {
                            aqi += aq.AQI;
                            pm25 += aq.PM25;
                            pm10 += aq.PM10;
                            no2 += aq.NO2;
                            co += aq.CO;
                            o3 += aq.O3;
                            so2 += aq.SO2;
                        }

                        aqi /= count; airQualityFeature[j].Add(aqi);
                        pm25 /= count; airQualityFeature[j].Add(pm25);
                        pm10 /= count; airQualityFeature[j].Add(pm10);
                        no2 /= count; airQualityFeature[j].Add(no2);
                        co /= count; airQualityFeature[j].Add(co);
                        o3 /= count; airQualityFeature[j].Add(o3);
                        so2 /= count; airQualityFeature[j].Add(so2);
                    }
                    else
                    {
                        airQualityFeature[j].Add(double.MaxValue);
                    }
                }

                string directory = Path.Combine(aqiFeatureDirectory, this.Days[i / Time.NUMSLOTPERDAY]);
                if (!Directory.Exists(directory))
                {
                    Directory.CreateDirectory(directory);
                }

                string filename = Path.Combine(directory, (i % Time.NUMSLOTPERDAY).ToString("D2"));
                IO.WriteFeature(airQualityFeature, filename);
            }
        }
    }

    public class AirQuality
    {
        public int StationID;
        public Time Timestamp;
        public double PM25;
        public double PM10;
        public double NO2;
        public double CO;
        public double O3;
        public double SO2;
        public int AQI;

        public AirQuality() { }

        public AirQuality(int sid, Time ts, double pm25, double pm10, double no2, double co, double o3, double so2, int aqi)
        {
            this.StationID = sid;
            this.PM25 = pm25;
            this.PM10 = pm10;
            this.NO2 = no2;
            this.CO = co;
            this.O3 = o3;
            this.SO2 = so2;
            this.AQI = aqi;

            this.Timestamp = new Time(ts.Year, ts.Month, ts.Day, ts.TimeSlot);
        }
    }
}
