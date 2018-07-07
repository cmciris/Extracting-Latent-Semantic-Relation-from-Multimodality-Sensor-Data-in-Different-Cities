using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace MMMST.Preprocess
{
    public class IO
    {
        /// <summary>
        /// read in the road network from the file
        /// </summary>
        /// <param name="roadnetworkFilename"></param>specifies the path to the filename of road network
        /// <returns></returns>
        public static RoadNetwork ReadRoadNetwork(string roadnetworkFilename)
        {
            Console.WriteLine("start to read the road network.");
            FileInfo fRead = new FileInfo(roadnetworkFilename);
            var totalSize = fRead.Length;
            var fs = fRead.OpenRead();
            StreamReader sr = new StreamReader(fs);
            long completeOld = 0;


            int nodeNum = Convert.ToInt32(sr.ReadLine());

            // add the node one by one
            List<Node> nodes = new List<Node>();
            int countNode = 0;
            while (countNode < nodeNum)
            {
                string[] items = sr.ReadLine().Split(new char[] { '\t' });
                Node node = new Node(Convert.ToInt32(items[0]), new LatLong(Convert.ToDouble(items[1]), Convert.ToDouble(items[2])));
                nodes.Add(node);
                long complete = fs.Position * 10 / totalSize;
                if (complete != completeOld)
                {
                    Console.WriteLine(complete * 10 + "% completed!");
                }
                completeOld = complete;
                countNode++;
            }

            //add the edge one by one
            int edgeNum = Convert.ToInt32(sr.ReadLine());
            List<Edge> edges = new List<Edge>();
            int countEdge = 0;
            while (sr.Peek() > 0)
            {
                string[] items = sr.ReadLine().Split(new char[] { '\t' });
                string[] points = sr.ReadLine().Split(new char[] { ',' });
                int startID = Convert.ToInt32(items[0]);
                int endID = Convert.ToInt32(items[1]);
                double length = Convert.ToDouble(items[2]);
                int roadClass = Convert.ToInt32(items[3]);
                switch(roadClass)
                {
                    case 41000:
                        roadClass = 0;
                        break;
                    case 42000:
                        roadClass = 1;
                        break;
                    case 51000:
                        roadClass = 2;
                        break;
                    case 52000:
                        roadClass = 3;
                        break;
                    case 53000:
                        roadClass = 4;
                        break;
                    case 54000:
                        roadClass = 5;
                        break;
                    case 43000:
                        roadClass = 6;
                        break;
                    case 44000:
                        roadClass = 7;
                        break;
                    case 45000:
                        roadClass = 8;
                        break;
                    case 47000:
                        roadClass = 9;
                        break;
                    case 49:
                        roadClass = 10;
                        break;
                    case 100:
                        roadClass = 11;
                        break;
                    default:
                        break;
                }
                int level = Convert.ToInt32(items[9]);
                int direction = Convert.ToInt32(items[4]);
                int formway = Convert.ToInt32(items[5]);
                int urban = Convert.ToInt32(items[6]);
                int maxLanes = Convert.ToInt32(items[7]);
                int maxSpeed = Convert.ToInt32(items[8]);

                List<LatLong> pointLatlong = new List<LatLong>();
                foreach(var p in points)
                {
                    string[] subs = p.Split(new char[] { ' ' });
                    pointLatlong.Add(new LatLong(Convert.ToDouble(subs[0]), Convert.ToDouble(subs[1])));
                }

                edges.Add(new Edge(countEdge, startID, endID, length, roadClass, direction, formway, urban, maxLanes, maxSpeed, level, 10, pointLatlong));
                countEdge++;

                long complete = fs.Position * 10 / totalSize;
                if (complete != completeOld)
                {
                    Console.WriteLine(complete * 10 + "% completed!");
                }
                completeOld = complete;
            }


            sr.Close();

            RoadNetwork roadNetwork = new RoadNetwork(nodeNum, edgeNum, nodes, edges);

            return roadNetwork;
        }

        /// <summary>
        /// read in all pois in the city level from the file
        /// </summary>
        /// <param name="cityPOIFilename"></param>specifies the path to the city poi filename
        /// <returns></returns>
        public static CityPOI ReadCityPOI(string cityPOIFilename)
        {
            Console.WriteLine("start to read the poi file.");
            FileInfo fRead = new FileInfo(cityPOIFilename);
            var totalSize = fRead.Length;
            var fs = fRead.OpenRead();
            StreamReader sr = new StreamReader(fs);
            long completeOld = 0;

            List<POI> pois = new List<POI>();
            int countPOI = 0;

            string[] attributes = sr.ReadLine().Split(new char[] { ',' });

            while(sr.Peek()>0)
            {
                string[] items = sr.ReadLine().Split(new char[] { ',' });
                //find the ID 
                int indexID = Array.IndexOf(attributes, "ID");
                //find the Category
                int indexCat = Array.IndexOf(attributes, "CATEGORY");
                //find the latitude and longitude
                int indexLat = Array.IndexOf(attributes, "Y_COORD");
                int indexLon = Array.IndexOf(attributes, "X_COORD");
                //find the rank
                int indexRank = Array.IndexOf(attributes, "RANK");

                POI poi = new POI(items[indexID], new LatLong(Convert.ToDouble(items[indexLat]), Convert.ToDouble(items[indexLon])), Convert.ToInt32(items[indexCat].Substring(0, 2)) - 1, Convert.ToInt32(items[indexRank]));
                pois.Add(poi);
                countPOI++;

                long complete = fs.Position * 10 / totalSize;
                if (complete != completeOld)
                {
                    Console.WriteLine(complete * 10 + "% completed!");
                }
                completeOld = complete;
            }

            sr.Close();
            CityPOI cityPOI = new CityPOI(countPOI, pois);
            return cityPOI;
        }

        /// <summary>
        /// read all files regarding to the mobility of the city in all time.
        /// </summary>
        /// <param name="roadnetworkFilename"></param>specifies the path of road network to read.
        /// <param name="roadMapFilename"></param>specifies the path to the map of mobility edge to the edge of road network.
        /// <param name="speedDirectory"></param>specifies the directory to all the speed information.
        /// <param name="pickDropPointsDirectory"></param>specifies the directory to the pickup and dropoff points.
        /// <returns></returns>
        public static ToTalMobility ReadMobility(string roadNetworkFilename, string roadMapFilename, string speedFilename, string pickDropPointsFilename, string dayName)
        {
            RoadNetwork roadNetwork = ReadRoadNetwork(roadNetworkFilename);

            Console.WriteLine("start to read the road map file.");
            StreamReader srMap = new StreamReader(roadMapFilename);
            Dictionary<int, List<int>> edgeRecords = new Dictionary<int, List<int>>();

            while (srMap.Peek() > 0)
            {
                string[] items = srMap.ReadLine().Split(new char[] { '\t' });
                if (Convert.ToInt32(items[1]) >= 0)
                {
                    if (!edgeRecords.ContainsKey(Convert.ToInt32(items[1])))
                    {
                        List<int> oriEdges = new List<int>();
                        oriEdges.Add(Convert.ToInt32(items[0]));
                        edgeRecords.Add(Convert.ToInt32(items[1]), oriEdges);
                    }
                    else
                    {
                        edgeRecords[Convert.ToInt32(items[1])].Add(Convert.ToInt32(items[0]));
                    }
                }
            }
            srMap.Close();

            Console.WriteLine("start to read the speed file.");


            Time[] timeStamps = new Time[Time.NUMSLOTPERDAY];
            List<PickDropPoint>[] points = new List<PickDropPoint>[Time.NUMSLOTPERDAY];
            List<MobilityEdge>[] edges = new List<MobilityEdge>[Time.NUMSLOTPERDAY];
            for (int i = 0; i < Time.NUMSLOTPERDAY; i++)
            {
                points[i] = new List<PickDropPoint>();
                edges[i] = new List<MobilityEdge>();
            }

            FileInfo fRead = new FileInfo(speedFilename);
            var totalSize = fRead.Length;
            var fs = fRead.OpenRead();
            StreamReader srSpeed = new StreamReader(fs);
            long completeOld = 0;

            while (srSpeed.Peek() > 0)
            {
                string[] items = srSpeed.ReadLine().Split(new char[] { ' ' });
                int edgeID = Convert.ToInt32(items[0]);
                //int timeSlot = Convert.ToInt32(items[1]);
                int timeSlot = Convert.ToInt32(items[1]) / 2;
                double avgs = Convert.ToDouble(items[2]);
                int numVehicles = Convert.ToInt32(items[3]);
                if (numVehicles < 0)
                {
                    numVehicles = 0;
                }
                double stds = Convert.ToDouble(items[4]);
                int year = Convert.ToInt32(dayName.Substring(0, 4));
                int month = Convert.ToInt32(dayName.Substring(4, 2));
                int day = Convert.ToInt32(dayName.Substring(6, 2));

                timeStamps[timeSlot] = new Time(year, month, day, timeSlot);
                foreach (var ind in edgeRecords[edgeID])
                {
                    Edge oriEdge = roadNetwork.Edges[ind];
                    edges[timeSlot].Add(new MobilityEdge(oriEdge.EdgeID, oriEdge.StartID, oriEdge.EndID, oriEdge.Length, oriEdge.RoadClass,
                                                oriEdge.Direction, oriEdge.FormWay, oriEdge.Urban, oriEdge.MaxLanes, oriEdge.MaxSpeed, oriEdge.Level,
                                                oriEdge.AttributeNum, oriEdge.Points, avgs, numVehicles, stds));
                }

                long complete = fs.Position * 10 / totalSize;
                if (complete != completeOld)
                {
                    Console.WriteLine(complete * 10 + "% completed!");
                }
                completeOld = complete;
            }
            srSpeed.Close();


            Console.WriteLine("start to read the pickdroppoint directory.");
            DirectoryInfo dirPoint = new DirectoryInfo(pickDropPointsFilename);
            FileInfo[] files2 = dirPoint.GetFiles();

            StreamReader srPoint = new StreamReader(pickDropPointsFilename);
            while (srPoint.Peek() > 0)
            {
                string[] items = srPoint.ReadLine().Split(new char[] { ' ', ':', '/' });
                int hour = Convert.ToInt32(items[3]);
                int min = Convert.ToInt32(items[4]);
                //int timeSlot = (min <= 30) ? (hour * 2) : (hour * 2 + 1);
                int timeSlot = hour;
                LatLong location = new LatLong(Convert.ToDouble(items[6]), Convert.ToDouble(items[7]));
                points[timeSlot].Add(new PickDropPoint(location, Convert.ToInt32(items[8])));
            }
            srPoint.Close();


            //start to built the totalModality
            Mobility[] mobilities = new Mobility[Time.NUMSLOTPERDAY];
            for (int i = 0; i < Time.NUMSLOTPERDAY; i++)
            {
                mobilities[i] = new Mobility(edges[i], points[i], timeStamps[i]);
            }
            ToTalMobility totalModality = new ToTalMobility(roadNetwork.Nodes, mobilities);
            return totalModality;
        }

        /// <summary>
        /// read all files regarding to the mobility of the city in all time.
        /// </summary>
        /// <param name="roadNetwork"></param>the road network that use.
        /// <param name="roadMapFilename"></param>specifies the path to the map of mobility edge to the edge of road network.
        /// <param name="speedDirectory"></param>specifies the directory to all the speed information.
        /// <param name="pickDropPointsDirectory"></param>specifies the directory to the pickup and dropoff points.
        /// <returns></returns>
        public static ToTalMobility ReadMobility(RoadNetwork roadNetwork, string roadMapFilename, string speedFilename, string pickDropPointsFilename, string dayName)
        {
            Console.WriteLine("start to read the road map file.");
            StreamReader srMap = new StreamReader(roadMapFilename);
            Dictionary<int, List<int>> edgeRecords = new Dictionary<int, List<int>>();

            while (srMap.Peek() > 0)
            {
                string[] items = srMap.ReadLine().Split(new char[] { '\t' });
                if (Convert.ToInt32(items[1]) >= 0)
                {
                    if (!edgeRecords.ContainsKey(Convert.ToInt32(items[1])))
                    {
                        List<int> oriEdges = new List<int>();
                        oriEdges.Add(Convert.ToInt32(items[0]));
                        edgeRecords.Add(Convert.ToInt32(items[1]), oriEdges);
                    }
                    else
                    {
                        edgeRecords[Convert.ToInt32(items[1])].Add(Convert.ToInt32(items[0]));
                    }
                }
            }
            srMap.Close();

            Console.WriteLine("start to read the speed file.");


            Time[] timeStamps = new Time[Time.NUMSLOTPERDAY];
            List<PickDropPoint>[] points = new List<PickDropPoint>[Time.NUMSLOTPERDAY];
            List<MobilityEdge>[] edges = new List<MobilityEdge>[Time.NUMSLOTPERDAY];
            for (int i = 0; i < Time.NUMSLOTPERDAY; i++)
            {
                points[i] = new List<PickDropPoint>();
                edges[i] = new List<MobilityEdge>();
            }

            FileInfo fRead = new FileInfo(speedFilename);
            var totalSize = fRead.Length;
            var fs = fRead.OpenRead();
            StreamReader srSpeed = new StreamReader(fs);
            long completeOld = 0;

            StreamWriter swSpeed = new StreamWriter(@"D:\Result\Temporary\SpeedDistribution\speed.txt");

            int count = 0;
            while (srSpeed.Peek() > 0)
            {
                string[] items = srSpeed.ReadLine().Split(new char[] { ' ' });
                int edgeID = Convert.ToInt32(items[0]);
                int timeSlot = Convert.ToInt32(items[1]) / 2;
                if(Convert.ToInt32(items[1]) == 1)
                {
                    ;
                }
                double avgs = Convert.ToDouble(items[2]);
                int numVehicles = Convert.ToInt32(items[4]);
                if (numVehicles < 0)
                {
                    numVehicles = 0;
                }
                double stds = Convert.ToDouble(items[3]);
                int year = Convert.ToInt32(dayName.Substring(0, 4));
                int month = Convert.ToInt32(dayName.Substring(4, 2));
                int day = Convert.ToInt32(dayName.Substring(6, 2));

                swSpeed.WriteLine(avgs);

                timeStamps[timeSlot] = new Time(year, month, day, timeSlot);
                foreach (var ind in edgeRecords[edgeID])
                {
                    Edge oriEdge = roadNetwork.Edges[ind];
                    edges[timeSlot].Add(new MobilityEdge(oriEdge.EdgeID, oriEdge.StartID, oriEdge.EndID, oriEdge.Length, oriEdge.RoadClass,
                                                oriEdge.Direction, oriEdge.FormWay, oriEdge.Urban, oriEdge.MaxLanes, oriEdge.MaxSpeed, oriEdge.Level,
                                                oriEdge.AttributeNum, oriEdge.Points, avgs, numVehicles, stds));
                }

                long complete = fs.Position * 10 / totalSize;
                if (complete != completeOld)
                {
                    Console.WriteLine(complete * 10 + "% completed!");
                }
                completeOld = complete;
                count++;
            }
            srSpeed.Close();
            swSpeed.Close();

            Console.WriteLine("start to read the pickdroppoint directory.");
            StreamReader srPoint = new StreamReader(pickDropPointsFilename);
            while (srPoint.Peek() > 0)
            {
                string[] items = srPoint.ReadLine().Split(new char[] { ' ', ':', '/' });
                int hour = Convert.ToInt32(items[3]);
                int min = Convert.ToInt32(items[4]);
                //int timeSlot = (min <= 30) ? (hour * 2) : (hour * 2 + 1);
                int timeSlot = hour;
                LatLong location = new LatLong(Convert.ToDouble(items[6]), Convert.ToDouble(items[7]));
                points[timeSlot].Add(new PickDropPoint(location, Convert.ToInt32(items[8])));
            }
            srPoint.Close();


            //start to built the totalModality
            Mobility[] mobilities = new Mobility[Time.NUMSLOTPERDAY];
            for (int i = 0; i < Time.NUMSLOTPERDAY; i++)
            {
                if (timeStamps[i] != null)
                {
                    mobilities[i] = new Mobility(edges[i], points[i], timeStamps[i]);
                }
            }
            ToTalMobility totalModality = new ToTalMobility(roadNetwork.Nodes, mobilities);
            return totalModality;
        }

        /// <summary>
        /// read the meterology data from the file.
        /// </summary>
        /// <param name="meterologyFilename"></param>specifies the path to the meterology data file.
        /// <param name="stationLocationFilename"></param>specifies the path to the station-location correspondence file.
        /// <returns></returns>
        public static TotalMeterology ReadMeterology(string meterologyFilename, string stationLocationFilename)
        {
            //decide the number of days
            StreamReader srPre = new StreamReader(meterologyFilename);
            int month_first = 0, month_old = 0;
            int day_first = 0, day_old = 0;
            int district_first = 0;

            int countDay = 0;
            List<string> timeValue = new List<string>();
            srPre.ReadLine();
            while (srPre.Peek() > 0)
            {
                string[] items = srPre.ReadLine().Split(new char[] { ',', '-', ' ' });
                int month = Convert.ToInt32(items[2]);
                int day = Convert.ToInt32(items[3]);
                int district = Convert.ToInt32(items[0]);

                if (month_first == 0)
                {
                    month_first = month;
                    day_first = day;
                    district_first = district;
                }
                if (month == month_first && day == day_first && district != district_first)
                {
                    break;
                }
                if (day != day_old || month != month_old)
                {
                    countDay++;
                    timeValue.Add(items[1] + items[2] + items[3]);
                }

                month_old = month;
                day_old = day;
            }
            srPre.Close();

            Console.WriteLine("start to read the meterology data.");
            FileInfo fRead = new FileInfo(meterologyFilename);
            var totalSize = fRead.Length;
            var fs = fRead.OpenRead();
            StreamReader sr = new StreamReader(fs);
            long completeOld = 0;

            double rainFall, rainFallOld = 0.0;
            int temperature, temperatureOld = 0;
            int pressure, pressureOld = 0;
            int humidity, humidityOld = 0;
            double windSpeed, windSpeedOld = 0;
            int windDirection, windDirectionOld = 0;
            int weatherFlag, weatherFlagOld = 0;

            List<Meterology>[] meterologies = new List<Meterology>[countDay * Time.NUMSLOTPERDAY];
            for (int i = 0; i < countDay * Time.NUMSLOTPERDAY; i++)
            {
                meterologies[i] = new List<Meterology>();
            }
            sr.ReadLine(); //read the attributes line
            while (sr.Peek() > 0)
            {
                string[] items = sr.ReadLine().Split(new char[] { ',', ':', ' ' });

                int districtID = Convert.ToInt32(items[0]);
                string[] time = items[1].Split(new char[] { '-' });
                int month = Convert.ToInt32(time[1]);
                int day = Convert.ToInt32(time[2]);
                int year = Convert.ToInt32(time[0]);
                int hour = Convert.ToInt32(items[2]);
                int min = Convert.ToInt32(items[3]);
                //int timeSlot = (min < 30) ? (hour * 2) : (hour * 2 + 1);
                int timeSlot = hour;
    
                rainFall = items[5] != "NULL" ? Convert.ToDouble(items[5]) : rainFallOld;
                temperature = items[6] != "NULL" ? (int)Convert.ToDouble(items[6]) : temperatureOld;
                pressure = items[7] != "NULL" ? (int)Convert.ToDouble(items[7]) : pressureOld;
                humidity = items[8] != "NULL" ? (int)Convert.ToDouble(items[8]) : humidityOld;
                windSpeed = items[9] != "NULL" ? Convert.ToDouble(items[9]) : windSpeedOld;
                windDirection = items[10] != "-1" && items[10] != "NULL" ? Convert.ToInt32(items[10]) : windDirectionOld;
                weatherFlag = items[11] != "NULL" ? Convert.ToInt32(items[11]) : weatherFlagOld;

                rainFallOld = rainFall;
                temperatureOld = temperature;
                pressureOld = pressure;
                humidityOld = humidity;
                windSpeedOld = windSpeed;
                windDirectionOld = windDirection;
                weatherFlagOld = weatherFlag;

                Meterology meterology = new Meterology(districtID, new Time(year, month, day, timeSlot), rainFall, temperature, pressure, humidity, windSpeed, windDirection, weatherFlag);
                //decide the index of the time
                int index = timeValue.IndexOf(time[0] + time[1] + time[2]);
                meterologies[index * Time.NUMSLOTPERDAY + timeSlot].Add(meterology);


                long complete = fs.Position * 10 / totalSize;
                if (complete != completeOld)
                {
                    Console.WriteLine(complete * 10 + "% completed!");
                }
                completeOld = complete;

            }
            sr.Close();

            Dictionary<LatLong, int> stationLocations = new Dictionary<LatLong, int>();
            StreamReader srStation = new StreamReader(stationLocationFilename);
            srStation.ReadLine();
            while (srStation.Peek() > 0)
            {
                string[] items = srStation.ReadLine().Split(new char[] { ',' });
                LatLong latlong = new LatLong(Convert.ToDouble(items[2]), Convert.ToDouble(items[3]));
                if (!stationLocations.ContainsKey(latlong))
                    stationLocations.Add(latlong, Convert.ToInt32(items[4]));
            }
            srStation.Close();

            TotalMeterology totalMeterology = new TotalMeterology(meterologies, stationLocations, timeValue);
            return totalMeterology;
        }

        
        public static TotalAirQuality ReadAirQuality(string airQualityFilename, string aqiFilename, string aqiStationFilename)
        {
            Console.WriteLine("start to read the air quality file.");
            StreamReader srAQI = new StreamReader(aqiFilename);
            int month_first = 0, month_old = 0;
            int day_first = 0, day_old = 0;
            int district_first = 0;
            bool stop = false;

            int countDay = 0;
            List<string> timeValue = new List<string>();

            Dictionary<Pair<int, Time>, int> aqis = new Dictionary<Pair<int, Time>, int>();
            srAQI.ReadLine();

            while(srAQI.Peek()>0)
            {
                string[] items = srAQI.ReadLine().Split(new char[] { ',', '/', ':', ' ' });
                int sid = Convert.ToInt32(items[0]);
                int month = Convert.ToInt32(items[1]);
                int day = Convert.ToInt32(items[2]);
                int year = Convert.ToInt32(items[3]);
                int hour = Convert.ToInt32(items[4]);
                if (hour == 12 && items[7] == "AM")
                    hour = 0;
                else if (items[7] == "PM" && hour != 12)
                    hour = hour + 12;


                Time time = new Time(year, month, day, hour);
                int aqi = Convert.ToInt32(items[8]);
                var pair = new Pair<int, Time>(sid, time);

                if(!aqis.ContainsKey(pair))
                    aqis.Add(pair, aqi);

                //if (!stop)
                //{
                //    if (month_first == 0)
                //    {
                //        month_first = month;
                //        day_first = day;
                //        district_first = sid;
                //    }
                //    if (month == month_first && day == day_first && sid != district_first)
                //    {
                //        stop = true;
                //    }
                //    if (day != day_old || month != month_old)
                //    {
                        timeValue.Add(Convert.ToInt32(items[3]).ToString("D2") + Convert.ToInt32(items[1]).ToString("D2") + Convert.ToInt32(items[2]).ToString("D2"));
                    //}

                    //month_old = month;
                    //day_old = day;
                //}
            }
            srAQI.Close();

            timeValue = timeValue.Distinct().ToList();
            countDay = timeValue.Count();

            FileInfo fRead = new FileInfo(airQualityFilename);
            var totalSize = fRead.Length;
            var fs = fRead.OpenRead();
            StreamReader sr = new StreamReader(fs);
            long completeOld = 0;

            List<AirQuality>[] airQualities = new List<AirQuality>[countDay * Time.NUMSLOTPERDAY];
            for (int i = 0; i < countDay * Time.NUMSLOTPERDAY; i++)
            {
                airQualities[i] = new List<AirQuality>();
            }
            sr.ReadLine();
            while (sr.Peek() > 0)
            {
                string[] items = sr.ReadLine().Split(new char[] { ',', '-', ':', ' ' });
                int sid = Convert.ToInt32(items[0]);
                int year = Convert.ToInt32(items[1]);
                int month = Convert.ToInt32(items[2]);
                int day = Convert.ToInt32(items[3]);
                int hour = Convert.ToInt32(items[4]);


                double pm25 = items[7] == "NULL" ? -1 : Convert.ToDouble(items[7]);
                double pm10 = items[8] == "NULL" ? -1 : Convert.ToDouble(items[8]);
                double no2 = items[9] == "NULL" ? -1 : Convert.ToDouble(items[9]);
                double co = items[10] == "NULL" ? -1 : Convert.ToDouble(items[10]);
                double o3 = items[11] == "NULL" ? -1 : Convert.ToDouble(items[11]);
                double so2 = items[12] == "NULL" ? -1 : Convert.ToDouble(items[12]);

                var pair = new Pair<int, Time>(sid, new Time(year, month, day, hour));
                int aqi = 0;
                try
                {
                    aqi = aqis[pair];
                }
                catch (Exception)
                {
                    aqi = int.MaxValue;
                }
                AirQuality airQuality = new AirQuality(sid, new Time(year, month, day, hour), pm25, pm10, no2, co, o3, so2, aqi);
                int index = timeValue.IndexOf(items[1] + items[2] + items[3]);
                airQualities[index * Time.NUMSLOTPERDAY + hour].Add(airQuality);

                long complete = fs.Position * 10 / totalSize;
                if (complete != completeOld)
                {
                    Console.WriteLine(complete * 10 + "% completed!");
                }
                completeOld = complete;
            }

            sr.Close();

            StreamReader srStation = new StreamReader(aqiStationFilename);
            Dictionary<int, LatLong> stations = new Dictionary<int, LatLong>();
            srStation.ReadLine();
            while(srStation.Peek()>0)
            {
                string[] items = srStation.ReadLine().Split(new char[] { '\t' });
                int sid = Convert.ToInt32(items[0]);
                LatLong latlong = new LatLong(Convert.ToDouble(items[2]), Convert.ToDouble(items[3]));
                stations.Add(sid, latlong);
            }
            srStation.Close();

            TotalAirQuality totalAirQuality = new TotalAirQuality(airQualities, stations, timeValue);
            return totalAirQuality;
        }
       

        /// <summary>
        /// write the feature matrix into the file specified
        /// </summary>
        /// <param name="featureMatrix"></param>the array encodes the instances, the list encodes the features per instance
        /// <param name="featureFilename"></param>specifies the path to the feature file to be saved
        public static void WriteFeature(List<double>[] featureMatrix, string featureFilename)
        {
            Console.WriteLine("start to output the feature matrix.");
            StreamWriter sw = new StreamWriter(featureFilename);
            int instNum = featureMatrix.Count();

            for (int i = 0; i < instNum; i++)
            {
                bool first = true;
                foreach(var feature in featureMatrix[i])
                {
                    if (first)
                    {
                        sw.Write(feature);
                        first = false;
                    }
                    else
                    {
                        sw.Write(" " + feature);
                    }
                }
                sw.WriteLine();
                //display the progress of writing
                if (instNum >= 10)
                {
                    if (i % (instNum / 10) == 0)
                    {
                        Console.WriteLine((i * 10 / (instNum / 10)) + "% writing has been completed.");
                    }
                }
            }
            sw.Close();
        }
    }

    public class Pair<T1, T2> : IEquatable<Pair<T1, T2>>
    {
        public T1 Value1;
        public T2 Value2;

        public Pair() { }

        public Pair(T1 v1, T2 v2)
        {
            this.Value1 = v1;
            this.Value2 = v2;
        }

        public bool Equals(Pair<T1, T2> other)
        {
            if (other == null)
                return false;
            if (this.Value1.Equals(other.Value1) && this.Value2.Equals(other.Value2))
                return true;
            else
                return false;
        }

        public override int GetHashCode()
        {
            return this.Value1.GetHashCode() ^ this.Value2.GetHashCode();
        }

        public override bool Equals(object obj)
        {
            Console.WriteLine("is executed!");
            if (obj == null)
                return false;
            Pair<T1, T2> pair = obj as Pair<T1, T2>;
            if (pair == null)
                return false;
            else
                return Equals(obj);
        }
    }
}
