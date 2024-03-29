#ifndef STRUCT_FILE
#define STRUCT_FILE
#include <vector>
using namespace std;
const int MAX_OBS_SATE_NUM =40;
const int MAX_OBS_TYPE = 4;
//时间结构体
typedef struct utctime{
	int year,month,day,hour,minute;
	double second;
}UTC;
typedef UTC * ptc;

typedef struct gpstime{
	int weeknum;
	long secondnum;
	double secondfrac;
}GPST;
typedef GPST * pgpst;

typedef struct bdtime{
	int week;
	double second;
}BDT;

typedef BDT * pbdt;
typedef struct juliantime{
	long daynum;
	int secondnum;
	double secondfrac;
}JULIANDAY;
typedef JULIANDAY * pjulian;

typedef struct doytime{
	int yearnum;
	int daynum;
	int secondnum;
	double secondfrac;
}DOY;
typedef DOY * pdoy;

//坐标结构体
typedef struct xyz_coordinate{
	double x;
	double y;
	double z;
}XYZ;
typedef XYZ * pxyz;

typedef struct enu_coordinate{
	double northing;
	double easting;
	double upping;
}ENU;
typedef ENU * penu;

typedef struct enupolar_coordinate{
	double range;
	double azimuth;
	double elevation;;
}ENUPOLAR;
typedef ENUPOLAR * penupolar;

typedef struct blh_coordinate{
	double longitude;
	double latitude;
	double height;
}BLH;
typedef BLH * pblh;
//sp3文件结构体
typedef struct single_sp3_ephemeris{//一颗卫星一个时间的xyz坐标
	UTC utime_n;
	double x, y, z;
}s_sp3_ephe;
typedef struct one_sate_ephemris{//一颗卫星所有时间的xyz坐标
	string prn;
	vector<s_sp3_ephe> sate_ephem;
}one_sate_ephem;
typedef struct all_sate_ephemeris{//所有卫星所有时间的xyz坐标
	one_sate_ephem all_ephem[32];//32颗卫星
}all_sate_ephem;
typedef all_sate_ephem* psp3;

typedef struct obs_head_data{
	XYZ approx_coordinate;
	ENU antena_height;
	string obsdata_type,station;
}ohd;

//对应vtec计算的o文件读取结构体
typedef struct observe_value_vtec{//一颗卫星一个时间的观测值
	UTC utime_o;
	double p1, p2, l1, l2;
}obs_value;
typedef struct observe_value{//一颗卫星所有时间的观测值
	string prn;
	vector<obs_value> obs_data;
}obs_one_sate;
typedef struct all_sate_obs{//所有卫星所有时间的观测值
	ohd obsheaddata;
	obs_one_sate all_sate[32];//32颗卫星
}obs;
typedef obs* pobs;

//vtec结果结构体
typedef struct one_vtec_result{
	//double sta_lat, sta_lon;//测站纬度、经度
	UTC rtime;//时间
	double lat, lon;//纬度、经度
	double vtec,vtec1, ionex_vtec;//vtec大小
	double mf;//投影函数
	string prn;
	//double p1, p2, p1s, p2s;
}one_result;
typedef struct result_out{
	vector<one_result> allresult;
}result;
typedef result * pvt;
//ionex文件结构体
//卫星dcb
typedef struct satellitedcb{
	string prn;
	double bias, rms;
}satedcb;
//测站dcb
typedef struct stationdcb{
	string station;
	double bias, rms;
}stadcb;
//vtec格网
typedef struct tecmaps{
	UTC vtime;
	double tec_values[71][73];
}vtecmap;
typedef struct ionex_file{
	vector<satedcb> sate_dcb;
	vector<stadcb> sta_dcb;
	vector<vtecmap> vtec_map;
	vector<vtecmap> vtec_rms;
}ionex;
typedef  ionex* pio;
//解dcb的方程组的系数和右端项
/*typedef struct dcb_equations{
	double E0, E1, E2, E3,dp;
}dcb;*/
//dcb文件结构体
typedef struct dcb_file{
	string prn,station,obs1,obs2;
	double bias, rms;
}dcb_f;
typedef struct all_dcb_file{
	vector<dcb_f> dcb_val;
}dcb;
typedef  dcb* pdcb;
#endif