#include"函数.h"
#include"stdafx.h"
#include<iostream>
#include<string>
#include"stdlib.h"
int _tmain(int argc, _TCHAR* argv[])
{
	//XYZ sta_coor;
	//BLH sta_coor_blh;//测站和穿刺点的大地坐标
	
	string strn = "E:\\GNSS\\数据\\igs19332.sp3";
	string stro = "E:\\GNSS\\数据\\cebr0240.17o";//cas1
	string stri = "E:\\GNSS\\数据\\codg0240.17i";
	string strd = "E:\\GNSS\\数据\\CAS0MGXRAP_20170240000_01D_01D_DCB.BSX";
	pobs po = new obs;
	psp3 pn = new all_sate_ephem;
	pio ion = new ionex;
	pdcb pd = new dcb;
	pvt pt = new result;
	readofile_vtec(stro, po);
	//sta_coor =po->obsheaddata.approx_coordinate;
	//xyztoblh(&sta_coor, &sta_coor_blh);
	//cout << "测站地理坐标："<<po->obsheaddata.station<<"（" << sta_coor_blh.latitude*180.0 / PI << "," << sta_coor_blh.longitude*180.0 / PI << ")" <<  endl;
	//return 0;
	readsp3file(strn, pn);
	read_ionex(stri, ion);
	read_dcb(strd, pd);
	vtec(po, pn,ion,pd, pt);
	putresult(pt);
	delete po;
	delete pn;
	delete ion;
	delete pd;
	delete pt;
	return 0;
}