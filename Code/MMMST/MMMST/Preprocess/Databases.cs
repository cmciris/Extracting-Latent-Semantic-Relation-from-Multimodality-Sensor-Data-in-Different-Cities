using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Data.SqlClient;
using System.IO;

namespace MMMST.Preprocess
{
    public class Databases
    {
        public void MeterologyStation()
        {
            SqlConnection conn = new SqlConnection(@"Data Source=urbcomp01;Initial Catalog=ruiyuan_test;Persist Security Info=True;User ID=sa;Password=abcd1234!");
            conn.Open();
            SqlCommand cmd = new SqlCommand("SELECT [weather_station_id],[name_chinese],[latitude],[longitude],[district_id] FROM [ruiyuan_test].[dbo].[WeatherStation]", conn);
            SqlDataReader reader = cmd.ExecuteReader();
            string cityID = "017"; //004 shenzhen 006 tianjin 022 langang 017 baoding

            StreamWriter swStation = new StreamWriter(@"D:\Data\Baoding\new meterolgy\station.csv");
            swStation.WriteLine("weather_station_id,name_chinese,latitude,longitude,district_id");
            while(reader.Read())
            {
                if (reader[0].ToString().Substring(0, 3) == cityID)
                {
                    swStation.WriteLine(reader[0] + "," + reader[1] + "," + reader[2] + "," + reader[3] + "," + reader[4]);
                }
            }
            swStation.Close();

            reader.Close();
            conn.Close();
        }

        public void Meterology()
        {
            SqlConnection conn = new SqlConnection(@"Data Source=urbcomp01;Initial Catalog=ruiyuan_test;Persist Security Info=True;User ID=sa;Password=abcd1234!");
            conn.Open();
            //read in the district id
            SqlCommand cmd1 = new SqlCommand("SELECT [district_id] FROM [ruiyuan_test].[dbo].[WeatherStation] WHERE district_id like '006%'", conn);
            SqlDataReader reader1 = cmd1.ExecuteReader();
            List<string> districtID = new List<string>();
            while(reader1.Read())
            {
                districtID.Add(reader1[0].ToString());
            }
            reader1.Close();

            districtID = districtID.Distinct().ToList();
            //districtID.Clear();
            //districtID.Add("017");

            StreamWriter sw = new StreamWriter(@"D:\Data\Tianjin\new meterolgy\meterology.csv");
            sw.WriteLine("id,time,RainFall,Temperature,Pressure,Humidity,WindSpeed,WindDirection,Weather");

            foreach (var did in districtID)
            {
                //find all the district ID's meterology
                SqlCommand cmd = new SqlCommand("SELECT [id],[time],[RainFall],[Temperature],[Pressure],[Humidity],[WindSpeed],[WindDirection],[Weather],[update_time],[Tag] FROM [ruiyuan_test].[dbo].[Meteorology]  WHERE id =" + did + " AND len(id)=5 AND time>='2014-08-01' AND time < '2014-12-01' AND update_time>='2014-08-01' AND update_time < '2014-12-01'", conn);
                SqlDataReader reader = cmd.ExecuteReader();

                int year_old = 0, month_old = 0, day_old = 0, hour_old = 0, tag_old = 0;
                string value = "", value_old = "";
                bool flag_old = true, flag = false;
                while (reader.Read())
                {
                    string[] nReader = new string[11];
                    for (int i = 0; i < 11; i++)
                    {
                        if (reader[i].ToString() == "")
                            nReader[i] = "NULL";
                        else
                            nReader[i] = reader[i].ToString();
                    }
                    string[] items;
                    if(reader[9].ToString() != "")
                    {
                        items = reader[9].ToString().Split(new char[] { '/', ' ', ':' });
                    }
                    else
                    {
                        items = reader[1].ToString().Split(new char[] { '/', ' ', ':' });
                    }
                    int month = Convert.ToInt32(items[0]);
                    int day = Convert.ToInt32(items[1]);
                    int year = Convert.ToInt32(items[2]);
                    int hour = Convert.ToInt32(items[3]);
                    if (hour == 12 && items[6] == "AM")
                        hour = 0;
                    else if (items[6] == "PM" && hour != 12)
                        hour = hour + 12;
                    int min = Convert.ToInt32(items[4]);
                    int sec = Convert.ToInt32(items[5]);
                    string time = year.ToString("D2") + "-" + month.ToString("D2") + "-" + day.ToString("D2") + " " + hour.ToString("D2") + ":" + min.ToString("D2") + ":" + sec.ToString("D2");
                    int tag = Convert.ToInt32(reader[10]);
                    switch (tag)
                    {
                        case 6:
                        case 7:
                            tag = 5;
                            break;
                        case 1:
                            tag = 4;
                            break;
                        case 0:
                            tag = 3;
                            break;
                        case 2:
                            tag = 2;
                            break;
                        case 4:
                        case 5:
                            tag = 1;
                            break;
                        case 3:
                            tag = 0;
                            break;
                        default:
                            break;
                    }

                    value = nReader[0] + "," + time + "," + nReader[2] + "," + nReader[3] + "," + nReader[4] + "," + nReader[5] + "," + nReader[6] + "," + nReader[7] + "," + nReader[8];
                    string toWrite = value;

                    flag = false;
                    if(year_old == year && month_old == month && day_old == day && hour_old == hour)
                    {
                        flag = true;
                        toWrite = tag >= tag_old ? value : value_old;
                    }

                    if (flag == true)
                        sw.WriteLine(toWrite);
                    if (flag == false && flag_old == false)
                    {
                        sw.WriteLine(value_old);
                    }

                    flag_old = flag;

                    value_old = value;
                    year_old = year;
                    month_old = month;
                    day_old = day;
                    hour_old = hour;
                    tag_old = tag;
                }
                sw.WriteLine(value);
                reader.Close();
            }
            sw.Close();
            conn.Close();
        }

        public void AQIStation()
        {
            SqlConnection conn = new SqlConnection(@"Data Source=urbcomp01;Initial Catalog=ruiyuan_test;Persist Security Info=True;User ID=sa;Password=abcd1234!");
            conn.Open();
            SqlCommand cmd = new SqlCommand("SELECT TOP 1000 [station_id],[name_Chinese],[latitude],[longitude],[district_id] FROM [ruiyuan_test].[dbo].[Station]", conn);
            SqlDataReader reader = cmd.ExecuteReader();
            string cityID = "006"; //004 shenzhen 006 tianjin 022 langang 017 baoding

            StreamWriter swStation = new StreamWriter(@"D:\Data\Tianjin\aqi\aqi station.txt");
            swStation.WriteLine("station_id\tname_Chinese\tlatitude\tlongitude\tdistrict_id");
            while (reader.Read())
            {
                if (reader[0].ToString().Substring(0, 3) == cityID)
                {
                    swStation.WriteLine(reader[0] + "\t" + reader[1] + "\t" + reader[2] + "\t" + reader[3] + "\t" + reader[4]);
                }
            }
            swStation.Close();

            reader.Close();
            conn.Close();
        }

        public void AQI()
        {
            SqlConnection conn = new SqlConnection(@"Data Source=urbcomp01;Initial Catalog=ruiyuan_test;Persist Security Info=True;User ID=sa;Password=abcd1234!");
            conn.Open();
            SqlCommand cmd = new SqlCommand("SELECT [station_id],[time],[PM25_Concentration],[PM10_Concentration],[NO2_Concentration],[CO_Concentration],[O3_Concentration],[SO2_Concentration],[update_time],[Tag] FROM [ruiyuan_test].[dbo].[AirQuality] WHERE station_id like '006%' AND len(station_id) = 6 AND time >= '2014-08-01' AND time < '2014-12-01'", conn);
            SqlDataReader reader = cmd.ExecuteReader();

            StreamWriter swStation = new StreamWriter(@"D:\Data\Tianjin\aqi\aqi.txt");
            StreamWriter sw = new StreamWriter(@"D:\Data\Tianjin\aqi\airquality.csv");
            swStation.WriteLine("station_id,time,AQI");
            sw.WriteLine("station_id,time,PM25_Concentration,PM10_Concentration,NO2_Concentration,CO_Concentration,O3_Concentration,SO2_Concentration");

            int year_old = 0, month_old = 0, day_old = 0, hour_old = 0, tag_old = 0;
            string value = "", value_old = "";
            string value2 = "", value2_old = "";
            bool flag_old = true, flag = false;

            while (reader.Read())
            {
                string[] nReader = new string[10];
                for (int i = 0; i < 10; i++)
                {
                    if (reader[i].ToString() == "")
                        nReader[i] = "NULL";
                    else
                        nReader[i] = reader[i].ToString();
                }

                double pm25 = 0.0;
                double pm10 = 0.0;
                double no2 = 0.0;
                double so2 = 0.0;
                if (reader[2].ToString() != "")//&& reader[3].ToString() != "" && reader[4].ToString() != "" && reader[5].ToString() != "")
                {
                    pm25 = Convert.ToDouble(reader[2]);
                    //pm10 = Convert.ToDouble(reader[3]);
                   // no2 = Convert.ToDouble(reader[4]);
                   // so2 = Convert.ToDouble(reader[5]);
                }
                else
                {
                    continue;
                }

                int aqi = CalculateAllAQI(pm25, pm10, no2, so2);

                string time = "";
                if (reader[8].ToString() != "")
                    time = reader[8].ToString();
                else
                    time = reader[1].ToString();

                string[] items = time.Split(new char[] { '/', ' ', ':' });
                int month = Convert.ToInt32(items[0]);
                int day = Convert.ToInt32(items[1]);
                int year = Convert.ToInt32(items[2]);
                int hour = Convert.ToInt32(items[3]);
                if (hour == 12 && items[6] == "AM")
                    hour = 0;
                else if (items[6] == "PM" && hour != 12)
                    hour = hour + 12;
                int min = Convert.ToInt32(items[4]);
                int sec = Convert.ToInt32(items[5]);

                string time2 = year.ToString("D2") + "-" + month.ToString("D2") + "-" + day.ToString("D2") + " " + hour.ToString("D2") + ":" + min.ToString("D2") + ":" + sec.ToString("D2");

                int tag = Convert.ToInt32(reader[9]);
                switch (tag)
                {
                    case 3:
                        tag = 4;
                        break;
                    case 4:
                        tag = 3;
                        break;
                    case 5:
                    case 0:
                        tag = 2;
                        break;
                    case 2:
                    case 1:
                        tag = 1;
                        break;
                    default:
                        break;
                }

                value = reader[0] + "," + time + "," + aqi;
                value2 = reader[0] + "," + time2 + "," + nReader[2] + "," + nReader[3] + "," + nReader[4] + "," + nReader[5] + "," + nReader[6] + "," + nReader[7];
                string toWrite = value;
                string toWrite2 = value2;

                flag = false;
                if (year_old == year && month_old == month && day_old == day && hour_old == hour)
                {
                    flag = true;
                    toWrite = tag >= tag_old ? value : value_old;
                    toWrite2 = tag >= tag_old ? value2 : value2_old;
                }

                if (flag == true)
                {
                    swStation.WriteLine(toWrite);
                    sw.WriteLine(toWrite2);
                }
                if (flag == false && flag_old == false)
                {
                    swStation.WriteLine(value_old);
                    sw.WriteLine(value2_old);
                }

                flag_old = flag;

                value_old = value;
                value2_old = value2;

                year_old = year;
                month_old = month;
                day_old = day;
                hour_old = hour;
                tag_old = tag;
            }
            swStation.WriteLine(value);
            sw.WriteLine(value2);

            swStation.Close();
            sw.Close();

            reader.Close();
            conn.Close();
        }

        private int CalculateAllAQI(double pm25, double pm10, double no2, double so2)
        {
            //return Math.Max(CalculateAQI(pm25), Math.Max(CalculateAQI(pm10), Math.Max(CalculateAQI(no2), CalculateAQI(so2))));
            return CalculateAQI(pm25);
        }

        private int CalculateAQI(double pollutant)
        {
            if (pollutant < 35)
                return (int)((pollutant - 0) * (50 - 0) / (35 - 0)) + 0;
            else if (pollutant >= 35 && pollutant < 75)
                return (int)((pollutant - 35) * (100 - 50) / (75 - 35)) + 50;
            else if (pollutant >= 75 && pollutant < 115)
                return (int)((pollutant - 75) * (150 - 100) / (115 - 75)) + 100;
            else if (pollutant >= 115 && pollutant < 150)
                return (int)((pollutant - 115) * (200 - 150) / (150 - 115)) + 150;
            else if (pollutant >= 150 && pollutant < 250)
                return (int)((pollutant - 150) * (300 - 200) / (250 - 150)) + 200;
            else if (pollutant >= 250 && pollutant < 350)
                return (int)((pollutant - 250) * (400 - 300) / (350 - 250)) + 300;
            else
                return (int)((pollutant - 350) * (500 - 400) / (500 - 350)) + 400;
        }
    }
}
