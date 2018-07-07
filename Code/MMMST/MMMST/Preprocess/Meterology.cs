using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MMMST.Utils;
using System.IO;

namespace MMMST.Preprocess
{
    public class TotalMeterology
    {
        public List<Meterology>[] Meterologies;
        public int NumTimeSlots;
        public Dictionary<LatLong, int> StationLocations;
        public List<string> Days;

        public TotalMeterology()
        { }

        public TotalMeterology(List<Meterology>[] m, Dictionary<LatLong, int> sl, List<string> days)
        {
            this.NumTimeSlots = m.Count();
            this.Meterologies = new List<Meterology>[this.NumTimeSlots];
            for (int i = 0; i < this.NumTimeSlots; i++)
            {
                this.Meterologies[i] = new List<Meterology>();
                foreach (var e in m[i])
                {
                    this.Meterologies[i].Add(e);
                }
            }
            this.StationLocations = sl;
            this.Days = new List<string>();
            foreach(var d in days)
            {
                this.Days.Add(d);
            }
        }

        public void TotalMeterologyFeature(Grid grid, string meterologyFeatureDirectory, int mode)
        {
            int gridNum = grid.GridNumLat * grid.GridNumLon;
            int countDay = this.NumTimeSlots / Time.NUMSLOTPERDAY;

            List<LatLong> stations = new List<LatLong>();
            foreach(var s in this.StationLocations)
            {
                stations.Add(s.Key);
            }
            List<double>[] meterologyFeatureOld = new List<double>[gridNum];


            for (int i = 0; i < this.NumTimeSlots; i++)
            {
                List<Meterology> meterology = this.Meterologies[i];
                List<double>[] meterologyFeature = new List<double>[gridNum];
                for (int g = 0; g < gridNum;g++ )
                {
                    meterologyFeatureOld[g] = new List<double>();
                }

                if (meterology.Count() > 0)
                {
                    for (int j = 0; j < gridNum; j++)
                    {
                        LatLong maxPoint, minPoint;
                        grid.RetrieveBoundary(j, out maxPoint, out minPoint);
                        LatLong central = new LatLong((maxPoint.Latitude - minPoint.Latitude) / 2 + minPoint.Latitude, (maxPoint.Longitude - minPoint.Longitude) / 2 + minPoint.Longitude);
                        LatLong nearest = Geo.FindKNN(central, stations, 1)[0];

                        int districtID = this.StationLocations[nearest];
                        meterologyFeature[j] = new List<double>();
                        int indexStation = grid.RetrieveIndex(nearest);

                        if (mode == 0 || (mode == 1 && indexStation == j))
                        {
                            int count = 0;
                            //feature 1: average rainfall 
                            double avgRainFall = 0.0;
                            //feature 2: average temperature
                            int avgTemperature = 0;
                            //feature 3: average pressure
                            int avgPressure = 0;
                            //feature 4: average humidity
                            int avgHumidity = 0;
                            //feature 5: average windspeed
                            double avgWindSpeed = 0.0;
                            //feautre 6: average wind direction
                            int avgWindDirection = 0;
                            //feature 7: the majority of weather flag
                            Dictionary<int, int> weatherFlag = new Dictionary<int, int>();

                            foreach (var mtrl in meterology)
                            {
                                if (mtrl.DistrictID == districtID) //infer the meterology record that the region has
                                {
                                    avgRainFall += mtrl.RainFall;
                                    avgTemperature += mtrl.Temperature;
                                    avgPressure += mtrl.Pressure;
                                    avgHumidity += mtrl.Humidity;
                                    avgWindSpeed += mtrl.WindSpeed;
                                    avgWindDirection += mtrl.WindDirection;
                                    if (!weatherFlag.ContainsKey(mtrl.WeatherFlag))
                                    {
                                        weatherFlag.Add(mtrl.WeatherFlag, 1);
                                    }
                                    else
                                    {
                                        weatherFlag[mtrl.WeatherFlag]++;
                                    }
                                    count++;
                                }
                            }

                            if (count > 0)
                            {
                                meterologyFeature[j].Add(avgRainFall / count);
                                meterologyFeature[j].Add(((double)avgTemperature) / count);
                                meterologyFeature[j].Add(((double)avgPressure) / count);
                                meterologyFeature[j].Add(((double)avgHumidity) / count);
                                meterologyFeature[j].Add(avgWindSpeed / count);
                                meterologyFeature[j].Add(((double)avgWindDirection) / count);

                                var flagSort = (from w in weatherFlag
                                                orderby w.Value
                                                descending
                                                select w.Key).ToList();
                                meterologyFeature[j].Add((double)flagSort[0]);
                            }
                            else
                            {
                                meterologyFeature[j].Add(Double.MaxValue);
                            }
                        }
                        else
                        {
                            meterologyFeature[j].Add(Double.MaxValue);
                        }

                    }
                    for (int j = 0; j < gridNum; j++)
                    {
                        meterologyFeatureOld[j] = new List<double>();
                        foreach (var mf in meterologyFeature[j])
                        {
                            meterologyFeatureOld[j].Add(mf);
                        }
                    }
                }
                else
                {
                    meterologyFeature = meterologyFeatureOld;
                }

                
                string directory = Path.Combine(meterologyFeatureDirectory, this.Days[i / Time.NUMSLOTPERDAY]);
                if (!Directory.Exists(directory))
                {
                    Directory.CreateDirectory(directory);
                }

                string filename = Path.Combine(directory, (i % Time.NUMSLOTPERDAY).ToString("D2"));
                IO.WriteFeature(meterologyFeature, filename);

            }
        }
    }

    public class Meterology
    {
        public int DistrictID;
        public Time Timestamp;
        public double RainFall;
        public int Temperature;
        public int Pressure;
        public int Humidity;
        public double WindSpeed;
        public int WindDirection;
        public int WeatherFlag;

        public Meterology() { }

        public Meterology(int did, Time ts, double rf, int temp, int prs, int hmd, double ws, int wd, int wf)
        {
            this.DistrictID = did;
            this.Timestamp = new Time(ts.Year, ts.Month, ts.Day, ts.TimeSlot);
            this.RainFall = rf;
            this.Temperature = temp;
            this.Pressure = prs;
            this.Humidity = hmd;
            this.WindSpeed = ws;
            this.WindDirection = wd;
            this.WeatherFlag = wf;
        }
        
    }
}
