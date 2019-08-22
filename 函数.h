#ifndef FUNCTION_H
#define FUNCTION_H
#include<string>
#include"结构体.h"
#include"常量.h"
//#include"matrix.h"
using namespace std;
//读sp3文件
void readsp3file(string strn, psp3 sp3file);
//按照vtec计算的方式读取o文件
void readofile_vtec(string stro, pobs obsfile);
//读取ionex文件
void read_ionex(string stri, pio ionfile);
//读取dcb文件
void read_dcb(string strd, pdcb dcbfile);
//GPST到UTC
void gpsttoutc(pgpst gt, ptc ut);
//JULIANDAY到BDT
void juliandaytobdt(pjulian ju, pbdt bt);
//UTC到BDT
void utctobdt(ptc ut, pbdt bt);
//UTC到GPST
void utctogpst(ptc ut,pgpst gt);

//UTC到JULIANDAY
void utctojulianday(ptc ut,pjulian ju);

//JULIANDAY到UTC
void juliandaytoutc(pjulian ju,ptc ut);

//GPST到JULIANDAY
void gpsttojulianday(pgpst gt, pjulian ju);

//JULIANDAY到GPST
void juliandaytogpst(pjulian ju, pgpst gt);
//求儒略日差值
double deltjulianday(ptc u1, ptc u2);
//儒略日变换
void transjulian(pjulian ju1, double * dt, pjulian ju2);
//

//XYZ到BLH
void xyztoblh(pxyz px,pblh pb);
//BLH到XYZ
void blhtoxyz(pblh pb, pxyz px);
//XYZ到ENU
void xyztoenu(pxyz pxcenter, pxyz px, penu pe);
//ENU到XZY
void enutoxyz(pxyz pxcenter,penu pe, pxyz px);
//ENU到ENUPOLAR
void enutoenupolar(penu pe, penupolar pep);
//求卫星方位角和高度角
void sate_azi_ele(pxyz pxcenter, pxyz px, penupolar pep);
//对流层延迟改正
void tropo(double mjd,pxyz px1, pxyz px2,double * dx);
//求气象数据
void get_nominal_metedata(double mjd, double lat, double lon,double dhgt, double* pres, double* temp, double* rhumi, double* undu);
//求较小值
int min(int a, int b);
//接收机天线高改正
void antena_height_correction(pxyz px1, pxyz px2, penu pe, double*p1, double*p2, double*pp1, double*pp2);
//精密星历计算卫星坐标
void cal_sate_coor(string prn, ptc ut, psp3 sp3file, bool &flag, pxyz coor);
//由平近点角计算偏近点角
double calE(double* M, double* e);
//由偏近点角求真近点角
double calf(double* E, double* e);
//旋转矩阵
//void rotamatrix(string str, double seta,Matrix r);
//输出结果
void putresult(pvt pt);
//vtec函数
void vtec(pobs obsfile, psp3 sp3file, pio ionfile, pdcb dcbfile, pvt pt);
//由ionex文件插值计算任意点的vtec
double ionex_vtec(pio ionfile, UTC vtime, double lat, double lon);//lat和lon是大地经纬度
//从dcb文件中找到对应的dcb值
bool find_dcb(bool b, string prn, string station, pdcb dcbfile, double* dcb_val);
//bool find_dcb(bool b, string prn, string station, pio ionfile, double* dcb_val);
#endif