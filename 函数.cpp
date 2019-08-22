#include"函数.h"
#include<math.h>
#include<iostream>
#include<fstream>
#include<vector>
//#include"matrix.h"
#include<iomanip>
#include<string>
#include <Eigen/Dense>
using Eigen::MatrixXd;
using namespace std;
/*#ifndef _NO_NAMESPACE
using namespace std;
using namespace math;
#define STD std
#else
#define STD
#endif

#ifndef _NO_TEMPLATE
typedef matrix<double> Matrix;
#else
typedef matrix Matrix;
#endif

#ifndef _NO_EXCEPTION
#  define TRYBEGIN()	try {
#  define CATCHERROR()	} catch (const STD::exception& e) { \
						cerr << "Error: " << e.what() << endl; }
#else
#  define TRYBEGIN()
#  define CATCHERROR()
#endif*/
//读取sp3文件
void readsp3file(string strn, psp3 sp3file)
{
	s_sp3_ephe onesp3;
	ifstream nfile(strn, ios::in);
	if (!nfile){
		cout << "sp3文件打开错误" << endl;
		exit(0);
	}
	cout << "开始读sp3文件" << endl;
	for (int i = 0; i < 32; i++){//prn初始化
		if (i < 9){
			sp3file->all_ephem[i].prn = "G0" + to_string(i + 1);
		}
		else{sp3file->all_ephem[i].prn = "G" + to_string(i + 1);}
	}
	string str1;
	getline(nfile, str1);
	while (str1.substr(0, 1) != "*")
	{
		getline(nfile, str1);
	}
	while (!nfile.eof())
	{
		if (str1.substr(0, 1) == "*"){
			onesp3.utime_n.year = atoi(str1.substr(3, 4).c_str());
			onesp3.utime_n.month = atoi(str1.substr(8, 2).c_str());
			onesp3.utime_n.day = atoi(str1.substr(11, 2).c_str());
			onesp3.utime_n.hour = atoi(str1.substr(14, 2).c_str());
			onesp3.utime_n.minute = atoi(str1.substr(17, 2).c_str());
			onesp3.utime_n.second = atof(str1.substr(20, 10).c_str());
		}
		else if(str1.substr(0, 2) == "PG"){
			onesp3.x = atof(str1.substr(5, 13).c_str())*1000.0;
			onesp3.y = atof(str1.substr(19, 13).c_str())*1000.0;
			onesp3.z = atof(str1.substr(33, 13).c_str())*1000.0;
			for (int i = 0; i < 32; i++){
				if (sp3file->all_ephem[i].prn == str1.substr(1, 3).c_str()){
					sp3file->all_ephem[i].sate_ephem.push_back(onesp3);
				}
			}
		}
		getline(nfile, str1);
	}
	cout << "sp3文件读取完毕" << endl;
}
/*
void readofile_vtec(string stro, pobs obsfile)//按照相位平滑伪距计算方式读取，即同一颗卫星放一起
{
	for (int i = 0; i < 32; i++){//prn初始化
		if (i < 9){
			obsfile->all_sate[i].prn = "G0" + to_string(i + 1);
		}
		else{ obsfile->all_sate[i].prn = "G" + to_string(i + 1); }
	}
	obs_value val;
	//GPS
	int pos_p1 = -1;//L1C在所给数据中的位置
	int pos_p2 = -1;//L2W在所给数据中的位置
	int pos_l1 = -1;//C1I在所给数据中的位置
	int pos_l2 = -1;//C6I在所给数据中的位置
	//BDS
	int pos_c1 = -1;//C1I在所给数据中的位置
	int pos_c2 = -1;//C6I在所给数据中的位置
	ifstream ofile(stro, ios::in);
	if (!ofile){
		cout << "o文件打开错误" << endl;
		exit(0);
	}
	//memset(&obsfile->obsdata, 0, sizeof(obsfile->obsdata));
	//cout << "开始读o文件" << endl;
	string str1, str2;
	getline(ofile, str1);
	//cout << str1 << endl;
	while (str1.substr(60, 13) != "END OF HEADER")
	{
		//cout << "开始读文件头"<< endl;
		str2 = str1.substr(60, 19);
		//cout <<str1<< endl;
		//GPS
		if (str1.substr(0, 1) == "G" && str1.substr(60, 19) == "SYS / # / OBS TYPES")
		{
			int a = atoi(str1.substr(4, 2).c_str());
			int b = 0;//防止一行数据写不完所有的测距码类型
			//cout << "a的值：" <<a<< endl;
			for (int i = 0; i < a; i++)
			{
				if (str1.substr(7 + 4 * (i-b), 3) == "C1C") pos_p1 = i;
				if (str1.substr(7 + 4 * (i-b), 3) == "C2W") pos_p2 = i;
				if (str1.substr(7 + 4 * (i-b), 3) == "L1C") pos_l1 = i;
				if (str1.substr(7 + 4 * (i-b), 3) == "L2W") pos_l2 = i;
				if (i == 12) getline(ofile, str1),b=12;
			}
			//cout << "p1位置的值：" << pos_p1 << endl;
			//cout << "p2位置的值：" << pos_p2 << endl;
			if (pos_p1 == -1)cout << "没有C1C数据" << endl;
			if (pos_p2 == -1)cout << "没有C2W数据" << endl;
			if (pos_l1 == -1)cout << "没有L1C数据" << endl;
			if (pos_l2 == -1)cout << "没有L2W数据" << endl;
		}
		//北斗
		if (str1.substr(0, 1) == "C" && str1.substr(60, 19) == "SYS / # / OBS TYPES")
		{
			int a = atoi(str1.substr(4, 2).c_str());
			int b = 0;
			//cout << "a的值：" <<a<< endl;
			for (int i = 0; i < a; i++)
			{
				if (str1.substr(7 + 4 * (i - b), 3) == "C1I" || str1.substr(7 + 4 * (i - b), 3) == "C2I") pos_c1 = i;
				if (str1.substr(7 + 4 * (i - b), 3) == "C7I" || str1.substr(7 + 4 * (i - b), 3) == "C6I") pos_c2 = i;
				if (i == 12) getline(ofile, str1), b = 12;
			}
			//cout << "p1位置的值：" << pos_p1 << endl;
			//cout << "p2位置的值：" << pos_p2 << endl;
			if (pos_c1 == -1)cout << "没有C1I和C2I数据" << endl;//C1I对应B1伪距
			if (pos_c2 == -1)cout << "没有C7I和C6I数据" << endl;//C6I对应B3伪距
			//C7I对应B2伪距
		}

		if (str1.substr(60, 20) == "RINEX VERSION / TYPE"){
			obsfile->obsheaddata.obsdata_type = str1.substr(40, 1);
			//cout << str1.substr(60, 19) << endl;
		}
		if (str1.substr(60, 19) == "APPROX POSITION XYZ"){

			obsfile->obsheaddata.approx_coordinate.x = atof(str1.substr(1, 14).c_str());
			obsfile->obsheaddata.approx_coordinate.y = atof(str1.substr(15, 14).c_str());
			obsfile->obsheaddata.approx_coordinate.z = atof(str1.substr(29, 14).c_str());
			//cout << obsfile->obsheaddata.approx_coordinate.x << endl;
		}
		if (str1.substr(60, 20) == "ANTENNA: DELTA H/E/N"){
			obsfile->obsheaddata.antena_height.upping = atof(str1.substr(7, 7).c_str());
			obsfile->obsheaddata.antena_height.easting = atof(str1.substr(21, 7).c_str());
			obsfile->obsheaddata.antena_height.northing = atof(str1.substr(35, 7).c_str());
			//cout << obsfile->obsheaddata.antena_height.upping << endl;
		}
		if (str1.substr(60, 11) == "MARKER NAME"){
			obsfile->obsheaddata.station = str1.substr(0, 4);
			//cout << obsfile->obsheaddata.antena_height.upping << endl;
		}
		getline(ofile, str1);
	}
	cout << "读取o文件头成功" << endl;
	while (!ofile.eof())
	{
		getline(ofile, str1);
		//if (str1.empty() == true) break;
		if (str1.substr(0, 1) == ">")
		{
			int sat = 0;
			val.utime_o.year = atoi(str1.substr(2, 4).c_str());
			val.utime_o.month = atoi(str1.substr(7, 2).c_str());
			val.utime_o.day = atoi(str1.substr(10, 2).c_str());
			val.utime_o.hour = atoi(str1.substr(13, 2).c_str());
			val.utime_o.minute = atoi(str1.substr(16, 2).c_str());
			val.utime_o.second = atof(str1.substr(19, 10).c_str());
			sat = atoi(str1.substr(32, 3).c_str());
			//cout << obsd.obs_n << endl;
			//cout << "观测时间：" << obsd.utime_o.hour << ":" << obsd.utime_o.minute << ":" << obsd.utime_o.second << obsd.obsvalue[sat][0] << endl;
			for (int i = 0; i < sat; i++)
			{
				//cout <<"i的值" <<i<< endl;
				getline(ofile, str1);
				//GPS			
				if (str1.substr(0, 1) == "G" && str1.length()>abs(16 * (1 + pos_p1)) && str1.length() > abs(16 * (1 + pos_p2))
					&& str1.length() > abs(16 * (1 + pos_l1)) && str1.length() > abs(16 * (1 + pos_l2)))//保证读取的数据长度够长，因为atof函数的代入值不能为空
				{
					val.p1 = atof(str1.substr(4 + 16 * pos_p1, 14).c_str());
					val.p2 = atof(str1.substr(4 + 16 * pos_p2, 14).c_str());
					val.l1 = atof(str1.substr(4 + 16 * pos_l1, 14).c_str());
					val.l2 = atof(str1.substr(4 + 16 * pos_l2, 14).c_str());
					if (val.p1 == 0 || val.p1< 20000000.0 || val.p2 == 0 || val.p2 < 20000000.0 || val.l1 == 0 || val.l1*c / p1 < 20000000.0 || val.l2 == 0 || val.l2*c / p2 < 20000000.0)continue;//==0时字符串为空，防止为空字符串，< 20000000.0证明相位值有问题，卫星轨道离地面高度是大于两万米的
					for (int j = 0; j < 32; j++)
					{
						if (obsfile->all_sate[j].prn == str1.substr(0, 3))
						{
							(obsfile->all_sate[j].obs_data.push_back(val));
						}
					}
				}
			}
		}
	}
	ofile.close();
	cout << "读取o文件成功 " << endl;
}*/
void readofile_vtec(string stro, pobs obsfile)//按照相位平滑伪距计算方式读取，即同一颗卫星放一起
{
	for (int i = 0; i < 32; i++){//prn初始化
		if (i < 9){
			obsfile->all_sate[i].prn = "G0" + to_string(i + 1);
		}
		else{ obsfile->all_sate[i].prn = "G" + to_string(i + 1); }
	}
	obs_value val;
	//GPS
	int pos[4] = { -1,-1,-1,-1};//测距码在所给数据中的位置
	double value;//对应的数值
	/*
	int pos_p1 = -1;//L1C在所给数据中的位置
	int pos_p2 = -1;//L2W在所给数据中的位置
	int pos_l1 = -1;//C1I在所给数据中的位置
	int pos_l2 = -1;//C6I在所给数据中的位置
	*/
	ifstream ofile(stro, ios::in);
	if (!ofile){
		cout << "o文件打开错误" << endl;
		exit(0);
	}
	//memset(&obsfile->obsdata, 0, sizeof(obsfile->obsdata));
	//cout << "开始读o文件" << endl;
	string str1,str_ins;//str_ins临时数据，用于str2和str3的pushback函数使用
	vector<string> str2, str3;//str2存储卫星prn，str3存储一颗卫星的所有数据
	getline(ofile, str1);
	//cout << str1 << endl;
	int a = 0;//所有的测距码类型
	while (str1.substr(60, 13) != "END OF HEADER")
	{
		//cout << "开始读文件头"<< endl;
		//cout <<str1<< endl;
		//GPS
		if (str1.substr(60, 19) == "# / TYPES OF OBSERV")
		{
			a = atoi(str1.substr(4, 2).c_str());
			int b = 0;//防止一行数据写不完所有的测距码类型
			//cout << "a的值：" <<a<< endl;
			for (int i = 0; i < a; i++)
			{
				if (str1.substr(10 + 6 * (i - b), 2) == "C1") pos[0] = i;
				if (str1.substr(10 + 6 * (i - b), 2) == "P2") pos[1] = i;
				if (str1.substr(10 + 6 * (i - b), 2) == "L1") pos[2] = i;
				if (str1.substr(10 + 6 * (i - b), 2) == "L2") pos[3] = i;
				if (i == 8) getline(ofile, str1), b = 9;//一行最多写9种类型的数据
			}
			/*cout << "c1位置的值：" << pos[0] << endl;
			cout << "p2位置的值：" << pos[1] << endl;
			cout << "l1位置的值：" << pos[2] << endl;
			cout << "l2位置的值：" << pos[3] << endl;*/
			if (pos[0] == -1)cout << "没有P1数据" << endl;
			if (pos[1] == -1)cout << "没有P2数据" << endl;
			if (pos[2] == -1)cout << "没有L1数据" << endl;
			if (pos[3] == -1)cout << "没有L2数据" << endl;
		}

		if (str1.substr(60, 20) == "RINEX VERSION / TYPE"){
			obsfile->obsheaddata.obsdata_type = str1.substr(40, 1);
			//cout << str1.substr(60, 19) << endl;
		}
		if (str1.substr(60, 19) == "APPROX POSITION XYZ"){

			obsfile->obsheaddata.approx_coordinate.x = atof(str1.substr(1, 14).c_str());
			obsfile->obsheaddata.approx_coordinate.y = atof(str1.substr(15, 14).c_str());
			obsfile->obsheaddata.approx_coordinate.z = atof(str1.substr(29, 14).c_str());
			//cout << obsfile->obsheaddata.approx_coordinate.x << endl;
		}
		if (str1.substr(60, 20) == "ANTENNA: DELTA H/E/N"){
			obsfile->obsheaddata.antena_height.upping = atof(str1.substr(7, 7).c_str());
			obsfile->obsheaddata.antena_height.easting = atof(str1.substr(21, 7).c_str());
			obsfile->obsheaddata.antena_height.northing = atof(str1.substr(35, 7).c_str());
			//cout << obsfile->obsheaddata.antena_height.upping << endl;
		}
		if (str1.substr(60, 11) == "MARKER NAME"){
			obsfile->obsheaddata.station = str1.substr(0, 4);
			//cout << obsfile->obsheaddata.antena_height.upping << endl;
		}
		getline(ofile, str1);
	}
	cout << "读取o文件头成功" << endl;
	//getline(ofile, str1);
	while (!ofile.eof())
	{
		//getline(ofile, str1);
		// << "str1:" << str1 << endl;
		//cout << "str1.length():" << str1.length() << endl;
		//if (str1.empty() == true && ofile.eof()) break;//最后一行为空行则为true，字符串为空时substr函数会报错
		do{
			getline(ofile, str1);
			//cout << "str1:" << str1 << endl;
		} while (str1.length() == 0 && !ofile.eof());//防止连续出现一个或者多个空行
		if (str1.empty() == true && ofile.eof()) break;//最后一行为空行则为true，字符串为空时substr函数会报错
		//cout << "str1:" << str1 << endl;
		if (str1.substr(1, 2) == "17"&&str1.length()>=58)// == "17"意味着读取的是17年的数据，待修改
		{
			//先清空vector变量
			str2.clear();
			int sat = 0;
			int b = 0;//防止一行数据写不完所有的测距码类型
			val.utime_o.year = atoi(str1.substr(1, 2).c_str())+2000;
			val.utime_o.month = atoi(str1.substr(4, 2).c_str());
			val.utime_o.day = atoi(str1.substr(7, 2).c_str());
			val.utime_o.hour = atoi(str1.substr(10, 2).c_str());
			val.utime_o.minute = atoi(str1.substr(13, 2).c_str());
			val.utime_o.second = atof(str1.substr(16, 10).c_str());
			sat = atoi(str1.substr(30, 2).c_str());
			for (int k = 0; k < sat / 12; k++){//一共sat/12行数据，一行12个
				getline(ofile, str_ins);
				str2.push_back(str_ins);
			}
			/*for (int iter = 0; iter < str2.size(); iter++){
				cout << "str2[" << iter << "]" << str2[iter] << endl;
			}*/
			//cout << sat / 12 + 1 << endl;
			//if (sat > 12) getline(ofile, str2);//一行数据写不下所有的观测到的卫星
			//cout << obsd.obs_n << endl;
			//cout << "观测时间：" << val.utime_o.year << ":" << val.utime_o.month << ":" << val.utime_o.day << ":" << val.utime_o.hour << ":" << val.utime_o.minute << ":" << val.utime_o.second << endl;
			for (int i = 0; i < sat; i++)
			{
				//先清空vector变量
				str3.clear();
				val.p1 = 0.0; val.p2 = 0.0; val.l1 = 0.0; val.l2 = 0.0;
				if (i % 12 == 0 && i != 0){ b = i; str1 = str2[i / 12 - 1];};
				//cout <<"i的值" <<i<< endl;
				for (int k = 0; k < a / 5+1; k++){//一共a/5行数据，一行5种类型数据
					getline(ofile, str_ins);
					str3.push_back(str_ins);
				}
				/*for (int iter = 0; iter < str3.size(); iter++){
					cout << "str3["<<iter<<"]" << str3[iter] << endl;
				}*/
				/*cout << "str3:" << str3 << endl;
				cout << "str4:" << str4 << endl;
				cout << "str3.length():" << str3.length() << endl;
				cout << "str4.length():" << str4.length() << endl;
				//GPS
				cout << "(i+1)%12:" << (i + 1) % 12 << endl;
				cout << "b:" << b<< endl;
				cout << "32 + (i - b) * 3:" << 32 + (i - b) * 3 << endl;
				cout << "str1:" << str1 << endl;
				cout << "str1.substr(32 + (i - b) * 3, 1):" << str1.substr(32 + (i - b) * 3, 1) << endl;*/
				//cout << "str1:" << str1 << endl;
				if (str1.substr(32 + (i - b) * 3, 1) != "G" ) continue;//保证读取的数据为GPS
				//cout << "prn:" << str1.substr(32 + (i - b) * 3, 3) << endl;
				for (int p = 0; p < 4; p++)
				{
					if (str3[pos[p]/5].length() < 14 + 16 * (pos[p]%5))break;//保证读取的数据长度够长，因为atof函数的代入值不能为空
					value = atof(str3[pos[p] / 5].substr(1 + 16 * (pos[p]%5), 14).c_str());//pos[p]种类型数据所在行数为pos[p]/5，列数为pos[p]%5
					switch (p)
					{
					case 0:val.p1 = value; break;
					case 1:val.p2 = value; break;
					case 2:val.l1 = value; break;
					case 3:val.l2 = value; break;
					}
				}
				//cout << "hhh" << endl;
				/*if (val.utime_o.hour == 9 && val.utime_o.minute == 26 && str1.substr(32 + (i - b) * 3, 3) == "G07"){
					cout.setf(ios::fixed);
					cout << "val.p1:" << val.p1 << endl;
					cout << "val.p2:" << val.p2 << endl;
					cout << "val.l1:" << val.l1 << endl;
					cout << "val.l2:" << val.l2 << endl;
				}*/
				if (val.p1 == 0 || val.p1 < 20000000.0 || val.p2 == 0 || val.p2 < 20000000.0 || val.l1 == 0 || val.l2 == 0) continue;
				//if (val.p1 == 0 || val.p1 < 20000000.0 || val.p2 == 0 || val.p2 < 20000000.0 || val.l1 == 0 || val.l1*c / f1 < 20000000.0 || val.l2 == 0 || val.l2*c / f2 < 20000000.0) continue;//==0时字符串为空，防止为空字符串，< 20000000.0证明相位值有问题，卫星轨道离地面高度是大于两万米的
				//if (val.p1 < 20000000.0 || val.p2 < 20000000.0 ||  val.l1*c / f1 < 20000000.0  || val.l2*c / f2 < 20000000.0)cout << "数据太小" << endl;
				for (int j = 0; j < 32; j++)
				{
					//cout << "prn:" << str1.substr(32 + (i - b) * 3, 3) << endl;
					//cout << "obsfile->all_sate[j].prn:" << obsfile->all_sate[j].prn << endl;
					if (atoi(obsfile->all_sate[j].prn.substr(1, 2).c_str()) == atoi(str1.substr(32 + (i - b) * 3+1, 2).c_str()))//有的测站给出的卫星prn号是G 1的形式
					{
						(obsfile->all_sate[j].obs_data.push_back(val));
					}
				}
			}
		}
		//getline(ofile, str1);//放在最后一行防止出现读取空行用substr函数判断出错
	}
	ofile.close();
	cout << "读取o文件成功 " << endl;
}
void read_ionex(string stri, pio ionfile)//按照相位平滑伪距计算方式读取，即同一颗卫星放一起
{
	cout << "开始读ionex文件 " << endl;
	satedcb dcb_sate;
	stadcb dcb_sta;
	vtecmap vtec_val;
	//vtecmap val_rms;
	ifstream ofile(stri, ios::in);
	int lat = 0;//纬度个数
	int num;//一行数据包含的点数，第五行9个，其余16个，共73个
	if (!ofile){
		cout << "ionex文件打开错误" << endl;


		exit(0);
	}
	//memset(&obsfile->obsdata, 0, sizeof(obsfile->obsdata));
	//cout << "开始读o文件" << endl;
	string str1, str2;
	getline(ofile, str1);
	//cout << str1 << endl;
	while (str1.substr(60, 13) != "END OF HEADER")
	{
		if (str1.substr(60, 16) == "PRN / BIAS / RMS")
		{
			//cout << "卫星dcb：" << str1 << endl;
			dcb_sate.prn = str1.substr(3, 3);
			dcb_sate.bias = atof(str1.substr(10, 6).c_str());
			dcb_sate.rms = atof(str1.substr(20, 6).c_str());
			ionfile->sate_dcb.push_back(dcb_sate);
		}
		if (str1.substr(60, 20) == "STATION / BIAS / RMS")
		{
			dcb_sta.station = str1.substr(6, 4);
			dcb_sta.bias = atof(str1.substr(30, 6).c_str());
			dcb_sta.rms = atof(str1.substr(40, 6).c_str());
			ionfile->sta_dcb.push_back(dcb_sta);
		}
		getline(ofile, str1);
	}
	cout << "读取ionex文件头成功" << endl;
	//cout <<str1<< endl;
	while (str1.substr(60, 11) != "END OF FILE")//getline函数遇到空行不再往下读
	{
		getline(ofile, str1);//跳过END OF TEC MAP或者END OF RMS MAP
		lat = 0;
		if (str1.substr(60, 16) == "START OF TEC MAP" || str1.substr(60, 16) == "START OF RMS MAP")//新的时间点对应的全天vtec
		{
			str2 = str1;//用于判断是START OF TEC MAP还是START OF RMS MAP
			getline(ofile, str1);
			if (str1.substr(60, 20) == "EPOCH OF CURRENT MAP"){
				//cout << str1 << endl;
				vtec_val.vtime.year = atoi(str1.substr(2, 4).c_str());
				vtec_val.vtime.month = atoi(str1.substr(10, 2).c_str());
				vtec_val.vtime.day = atoi(str1.substr(16, 2).c_str());
				vtec_val.vtime.hour = atoi(str1.substr(22, 2).c_str());
				vtec_val.vtime.minute = atoi(str1.substr(28, 2).c_str());
				vtec_val.vtime.second = atof(str1.substr(34, 10).c_str());
			}
			getline(ofile, str1);
			while (str1.substr(60, 20) == "LAT/LON1/LON2/DLON/H")//新的纬度对应的vtec
			{
				num = 16;
				for (int i = 0; i < 5; i++)
				{
					if (i == 4) num = 9;
					getline(ofile, str1);
					for (int j = 0; j < num; j++)
					{
						vtec_val.tec_values[lat][i * 16 + j] = atof(str1.substr(j * 5, 5).c_str());
					}
				}
				lat++;
				getline(ofile, str1);
			}
			if (str2.substr(60, 16) == "START OF TEC MAP"){
				ionfile->vtec_map.push_back(vtec_val);//vtec值
			}
			else{
				ionfile->vtec_rms.push_back(vtec_val);//vtec误差值
			}
			//getline(ofile, str1);//跳过END OF TEC MAP或者END OF RMS MAP
		}
	}
	ofile.close();
	cout << "读取ionex文件成功 " << endl;
}
void read_dcb(string strd, pdcb dcbfile){
	cout << "开始读dcb文件 " << endl;
	dcb_f one_dcb;//一行数据
	ifstream ofile(strd, ios::in);
	if (!ofile){
		cout << "dcb文件打开错误" << endl;
		exit(0);
	}
	string str1;
	int n = 0;
	//cout << str1 << endl;
	while (!ofile.eof())
	{
		getline(ofile, str1);
		//cout << "dcb文件列数：" << str1.length() << endl;
		//cout << "str1.substr(1, 3)：" << str1.substr(1, 3) << endl;
		//cout << str1 << endl;
		if ((str1.substr(1, 3) == "DSB" || str1.substr(1, 3) == "DCB") && str1.length()>100)
		{
			//cout << "dcb文件行数：" << ++n << endl;
			one_dcb.prn = str1.substr(11, 3);
			one_dcb.station = str1.substr(15, 4);
			//cout << "str1.substr(11, 3):" << one_dcb.prn << "  str1.substr(15, 4):" << one_dcb.station << endl;
			one_dcb.obs1 = str1.substr(25, 3);
			one_dcb.obs2 = str1.substr(30, 3);
			one_dcb.bias = atof(str1.substr(82, 9).c_str());
			one_dcb.rms = atof(str1.substr(94, 9).c_str());
			dcbfile->dcb_val.push_back(one_dcb);
		}
		if (str1.substr(0, 8) == "%=ENDBIA")break;
	}
	ofile.close();
	cout << "读取dcb文件成功 " << endl;
	cout << "dcb文件结果个数：" << dcbfile->dcb_val.size()<< endl;
	//for (int i = 0; i < dcbfile->dcb_val.size(); i++){
	//	cout << "prn" << dcbfile->dcb_val[i].prn << "station" << dcbfile->dcb_val[i].station << "obs1" << dcbfile->dcb_val[i].obs1 << "obs2" << dcbfile->dcb_val[i].obs2 << "dcb" << dcbfile->dcb_val[i].bias << endl;
//	}
}
void gpsttoutc(pgpst gt, ptc ut)
{
	pjulian ju;
	ju = (pjulian)malloc(sizeof(JULIANDAY));
	gpsttojulianday(gt,ju);
	juliandaytoutc(ju, ut);
	free(ju);
}
void juliandaytobdt(pjulian ju, pbdt bt)
{
	double dt = ju->daynum + (ju->secondfrac + ju->secondnum) / _DAY_IN_SECOND - 2453736.500000;
	bt->week = int(dt / 7);
	bt->second = (dt - bt->week*7)*_DAY_IN_SECOND;
}
void utctobdt(ptc ut, pbdt bt)
{
	if (ut->year < 2006 || ut->month>12 || ut->month < 0 || ut->day>31 || ut->day < 0 || ut->hour>24 || ut->hour < 0 || ut->minute>60 || ut->minute < 0 || ut->second>60 || ut->second < 0)
	{
		cout << "时间大小有误"<<endl;
	}
	pjulian ju;
	ju = (pjulian)malloc(sizeof(JULIANDAY));
	utctojulianday(ut,ju);
	juliandaytobdt(ju, bt);
	free(ju);
}

void utctogpst(ptc ut, pgpst gt)
{
	pjulian ju;
	ju = (pjulian)malloc(sizeof(JULIANDAY));
	utctojulianday(ut, ju);
	juliandaytogpst(ju, gt);
	free(ju);
}

void utctojulianday(ptc ut, pjulian ju)
{
	int		m;
	int		y;
	double	dhour;

	dhour = ut->hour + ut->minute / (double)_HOUR_IN_MINUTE
		+ ut->second / (double)_HOUR_IN_SECOND;

	if (ut->month <= 2) {
		y = ut->year - 1;
		m = ut->month + 12;
	}
	else {
		y = ut->year;
		m = ut->month;
	}

	ju->daynum = (long)(365.25*y) + (long)(30.6001*(m + 1))
		+ ut ->day + (long)(dhour / 24 + 1720981.5);
	ju->secondnum = ((ut->hour + 12) % _DAY_IN_HOUR)*_HOUR_IN_SECOND
		+ ut->minute*_MINUTE_IN_SECOND + (long)ut->second;
	ju->secondfrac = ut->second-(long)ut->second;
}

void juliandaytoutc(pjulian ju, ptc ut)
{
	int a, b, c, d, e;
	double JD;
	JD = ju->daynum + (ju->secondnum + ju->secondfrac) / _DAY_IN_SECOND;

	a = static_cast<int>(JD + 0.5);
	b = a + 1537;
	c = static_cast<int>((b - 122.1) / 365.25);
	d = static_cast<int>(365.25*c);
	e = static_cast<int>((b - d) / 30.6001);
	
	double day = b - d-(long)(30.6001*e) + JD + 0.5 - a;
	ut->day = int(day);
	ut->month = e - 1 - 12 * (int)(e / 14);
	ut->year = c - 4715 - (int)((7 + ut->month) / 10);

	ut->hour = int((day - ut->day)*24.0);
	ut->minute = (int)(((day - ut->day)*24.0 - ut->hour)*60.0);
	ut->second = ju->secondnum + ju->secondfrac - (int((ju->secondnum + ju->secondfrac)/60))*60.0;

}

void gpsttojulianday(pgpst gt, pjulian ju)
{
	double JD;
	JD = gt->weeknum * 7 + (gt->secondnum + gt->secondfrac) / _DAY_IN_SECOND + 2444244.5;
	ju->daynum= long(JD);
	
	ju->secondnum = long(gt->secondnum + (gt->weeknum * 7 + 2444244.5 - ju->daynum)*_DAY_IN_SECOND);
	ju->secondfrac = gt->secondfrac;
}


void juliandaytogpst(pjulian ju, pgpst gt)
{
	double JD;
	JD = ju->daynum +( ju->secondnum + ju->secondfrac) / _DAY_IN_SECOND;
	gt->weeknum = int((JD - 2444244.5) / 7);

	gt->secondnum = long((JD - 2444244.5 - gt->weeknum * 7)*_DAY_IN_SECOND);
	gt->secondfrac = ju->secondfrac;

}

//求儒略日差值
double deltjulianday(ptc u1, ptc u2)
{
	JULIANDAY j1, j2;
	utctojulianday(u1, &j1);
	utctojulianday(u2, &j2);
	double delt,d1,d2;
	d1 = j1.daynum + (j1.secondnum + j1.secondfrac) / _DAY_IN_SECOND;
	d2 = j2.daynum + (j2.secondnum + j2.secondfrac) / _DAY_IN_SECOND;
	//delt = (ju1->daynum - ju2->daynum)*_DAY_IN_SECOND + (ju1->secondnum - ju2->secondnum) + (ju1->secondfrac - ju2->secondfrac);
	delt = (d1 - d2)*_DAY_IN_SECOND;

	/*if (delt>302400)
		delt -= 604800;
	else if (delt<-302400)
		delt += 604800;
	else
		delt = delt;*/
	return delt;
}

//儒略日变换
void transjulian(pjulian ju1, double * dt, pjulian ju2)
{
	double JDold, JDnew;
	JDold = ju1->daynum + (ju1->secondnum + ju1->secondfrac) / _DAY_IN_SECOND;
	JDnew = JDold-*dt / _DAY_IN_SECOND;

	ju2->daynum = long(JDnew);
	ju2->secondnum = long((JDnew - long(JDnew))*_DAY_IN_SECOND);
	ju2->secondfrac = (JDnew - long(JDnew))*_DAY_IN_SECOND
		- long((JDnew - long(JDnew))*_DAY_IN_SECOND);
}






//坐标函数

//空间直角坐标系到大地坐标系
void xyztoblh(pxyz px, pblh pb)
{
	double pi = 4.0*atan(1.0);
	double E2 = 2.0*flattening - flattening * flattening;
	double E4 = E2*E2;
	double ALFA = (px->x*px->x + px->y*px->y + (1.0 - E2)*px->z*px->z) / (a*a);
	double BATA = (px->x*px->x + px->y*px->y - (1.0 - E2)*px->z*px->z) / (a*a);
	double Q = 1.0 + 13.50*E4*(ALFA*ALFA - BATA*BATA) / pow(ALFA - E4, 3);
	double A1 = -Q + sqrt(Q*Q - 1.0);
	double AL = (1.0 / 3.0)*log(-A1);
	AL = -exp(AL);
	double A2 = AL + 1.0 / AL;
	double A3 = AL - 1.0 / AL;
	double T23 = (ALFA + E4 / 2.0) / 3.0 - (ALFA - E4)*A2 / 12.0;
	double T32 = sqrt(T23 *T23 + ((ALFA - E4)*A3)*((ALFA - E4)*A3) / 48.0);
	double T1 = -E2*BATA / (4.0*T32);
	double DK = T1 + sqrt(T23 + T32) - (1.0 - E2 / 2.0);
	double EK = (1.0 + DK) / (1.0 + DK - E2);
	/*cout << "T1的值：" << T1 << endl;
	cout << "T32的值：" << T32 << endl;
	cout << "T23的值：" << T23 << endl;
	cout << "E2的值：" << E2 << endl;
	cout << "ALFA的值：" << ALFA << endl;
	cout << "AL的值：" << AL << endl;
	cout << "DK的值：" << DK << endl;
	cout << "EK的值：" << EK << endl;*/
	pb->height = (DK / (1.0 + DK))*sqrt(px->x*px->x + px->y*px->y + (EK*px->z)*(EK*px->z));
	double P = sqrt(px->x*px->x + px->y*px->y);
	pb->latitude = atan(EK*px->z / P);
	double COSFL = px->x / P;
	double SINFL = px->y / P;
	pb->longitude = asin(SINFL);
	if (SINFL > 0.0&&COSFL<0.0) pb->longitude = pi - pb->longitude;
	if (SINFL < 0.0&&COSFL>0.0) pb->longitude = 2.0*pi + pb->longitude;
	if (SINFL <0.0&&COSFL<0.0) pb->longitude = pi - pb->longitude;

	/*double e2;//第一偏心率的平方
	e2 = 2 * flattening - flattening*flattening;

	pb->longitude = atan(px->y / px->x);
	double W, N, N1 = 0, B, B1;
	B1 = atan(px->z / sqrt(px->x*px->x + px->y*px->y));
	while (1)
	{
		W = sqrt(1 - e2*sin(B1)*sin(B1));
		N1 = a / W;
		B = atan((px->z + N1*e2*sin(B1)) / sqrt(px->x*px->x + px->y*px->y));

		if (fabs(B - B1)<delta)
			break;
		else
			B1 = B;
	}

	pb->latitude = B;
	N = a / sqrt(1 - e2*sin(pb->latitude)*sin(pb->latitude));
	pb->height = sqrt(px->x*px->x + px->y*px->y) / cos(B) - N;*/
}

//大地坐标系到空间直角坐标系
void blhtoxyz(pblh pb, pxyz px)
{
	double e2;//第一偏心率的平方
	double N;//卯酉圈半径
	e2 = 2 * flattening - flattening*flattening;
	N = a / sqrt(1 - e2*sin(pb->latitude)*sin(pb->latitude));

	px->x = (N + pb->height)*cos(pb->latitude)*cos(pb->longitude);
	px->y = (N + pb->height)*cos(pb->latitude)*sin(pb->longitude);
	px->z = (N*(1 - e2) + pb->height)*sin(pb->latitude);
}

//笛卡尔坐标系到站心空间直角坐标系
void xyztoenu(pxyz pxcenter, pxyz px, penu pe)
{
	double dx, dy, dz;
	dx = px->x - pxcenter->x;
	dy = px->y - pxcenter->y;
	dz = px->z - pxcenter->z;

	pblh pd;
	pd = (pblh)malloc(sizeof(BLH));

	xyztoblh(pxcenter, pd);

	pe->northing = -sin(pd->latitude)*cos(pd->longitude)*dx
		- sin(pd->latitude)*sin(pd->longitude)*dy
		+ cos(pd->latitude)*dz;
	pe->easting = -sin(pd->longitude)*dx
		+ cos(pd->longitude)*dy;
	pe->upping = cos(pd->latitude)*cos(pd->longitude)*dx
		+ cos(pd->latitude)*sin(pd->longitude)*dy
		+ sin(pd->latitude)*dz;
	free(pd);
}

//站心空间直角坐标系到笛卡尔坐标系
 void enutoxyz(pxyz pxcenter, penu pe, pxyz px)
{
	pblh pd;
	pd = (pblh)malloc(sizeof(BLH));
	xyztoblh(pxcenter, pd);
	MatrixXd H(3, 3), DB(3, 1), DX(3, 1);
	DB(0, 0) = pe->northing;
	DB(1, 0) = pe->easting;
	DB(2, 0) = pe->upping;
	H(0, 0) = -sin(pd->latitude)*cos(pd->longitude);
	H(0, 1) = -sin(pd->latitude)*sin(pd->longitude);
	H(0, 2) = cos(pd->latitude);
	H(1, 0) = -sin(pd->longitude);
	H(1, 1) = cos(pd->longitude);
	H(1, 2) = 0;
	H(2, 0) = cos(pd->latitude)*cos(pd->longitude);
	H(2, 1) = cos(pd->latitude)*sin(pd->longitude);
	H(2, 2) = sin(pd->latitude);
	DX = (H.inverse())*DB;
	double dx, dy, dz;
	dx = DX(0,0 );
	dy = DX(1, 0);
	dz = DX(2, 0);
	px->x = pxcenter->x + dx;
	px->y = pxcenter->y + dy;
	px->z = pxcenter->z + dz;
	free(pd);
}

 //站心直角坐标系到站心极坐标系
 void enutoenupolar(penu pe, penupolar pep)
 {
	 pep->range = sqrt(pe->northing*pe->northing + pe->easting*pe->easting + pe->upping*pe->upping);
	 pep->azimuth = atan(pe->easting / pe->northing);
	 //atan2函数返回的范围为(-PI,PI]，当返回值大于零时，在1,2象限，小于零时在3,4象限
	 pep->azimuth = atan2(pe->easting / sqrt(pe->northing*pe->northing + pe->easting*pe->easting), pe->northing / sqrt(pe->northing*pe->northing + pe->easting*pe->easting));
	 if (pep->azimuth < 0.0)pep->azimuth += PI*2.0;
	 pep->elevation = atan(pe->upping / sqrt(pe->northing*pe->northing + pe->easting*pe->easting));
 }
 //求卫星高度角和方位角
 void sate_azi_ele(pxyz pxcenter, pxyz px, penupolar pep){
	 penu pe;
	 pe = (penu)malloc(sizeof(ENU));
	 xyztoenu(pxcenter,  px, pe);
	 enutoenupolar(pe, pep);
	 //cout << "卫星方位角：" << pep->azimuth*180.0/PI << endl;
	 //cout << "卫星高度角：" << pep->elevation*180.0 / PI << endl;
	 //cout << "卫星距离：" << pep->range<< endl;
	 free(pe);
 }
 //求穿刺点地理坐标以及投影角
 void ipp_pos(pblh pb, penupolar pep,pblh pb1,double* MF){
	 double psi;//张角
	 psi = PI / 2.0 - pep->elevation - asin(ave_a / (ave_a + hion)*cos(pep->elevation));
	 pb1->latitude = asin(sin(pb->latitude)*cos(psi) + cos(pb->latitude)*sin(psi)*cos(pep->azimuth));
	 pb1->longitude = pb->longitude + atan(cos(pb->latitude)*sin(psi)*sin(pep->azimuth) / (cos(psi) - sin(pb->latitude)*sin(pb1->latitude)));
	 if (pb1->longitude < 0.0)pb1->longitude += PI*2.0;
	 if (pb1->longitude > PI*2.0)pb1->longitude -= PI*2.0;
	 pb1->height = hion;
	 *MF = 1.0 / sqrt(1.0 - (ave_a / (ave_a + hion)*cos(pep->elevation)*ave_a / (ave_a + hion)*cos(pep->elevation)));
	 //cout << "穿刺点纬度：" << pb1->latitude*180.0 / PI << endl;
	 //cout << "穿刺点经度：" << pb1->longitude*180.0 / PI << endl;
	 //cout << "倾斜角：" << *MF << endl;
 }
 //由ionex文件插值计算任意点的vtec
 double ionex_vtec(pio ionfile, UTC vtime, double lat, double lon){//lat和lon是大地经纬度
	 double p, q;//权
	 double v,v1, v2;//前后两个时间点的插值大小
	 //大地经纬度转为地心经纬度
	 double lat_geo, lon_geo;
	 lon_geo = lon;
	 if (lon_geo > 180.0)lon_geo -= 360.0;
	 lat_geo = atan((1 - flattening)*(1 - flattening)*tan(lat*PI/180.0))*180.0/PI;
	 if (fabs(lon_geo) > 180.0 || fabs(lat_geo) > 87.5 || deltjulianday(&vtime, &ionfile->vtec_map[0].vtime)<0){
		 cout << "穿刺点不在计算范围内" << endl;
		 return 0.0;
	 }
	 for (int i = 0; i < ionfile->vtec_map.size(); i++){
		 if (deltjulianday(&vtime, &ionfile->vtec_map[i].vtime) <= 0.0)//找到插值时间点后面的历元
		 {
			 for (int j = 0; j < 71; j++)
			 {
				 if (lat_geo > 87.5 - j*2.5)//找到插值点纬度后面的纬度节点
				 {
					 for (int k = 0; k < 73; k++)
					 {
						 if (lon_geo < -180.0 + k *5.0)//找到插值点经度后面的经度节点
						 {
							 p = fabs((lat_geo - (87.5-j*2.5))/2.5);
							 q = fabs((lon_geo -(-180.0+(k-1)*5.0)) / 5.0);
							//空间插值
							 //前一时间点
							 v1 = (1 - p)*(1 - q)*ionfile->vtec_map[i-1].tec_values[j][k - 1] +
								  (1 - p)*q*ionfile->vtec_map[i-1].tec_values[j][k] +
								   p*(1 - q)*ionfile->vtec_map[i-1].tec_values[j-1][k - 1] +
								   p*q*ionfile->vtec_map[i-1].tec_values[j-1][k];
							 //后一时间点
							 v2 = (1 - p)*(1 - q)*ionfile->vtec_map[i].tec_values[j][k - 1] +
								 (1 - p)*q*ionfile->vtec_map[i].tec_values[j][k] +
								 p*(1 - q)*ionfile->vtec_map[i].tec_values[j - 1][k - 1] +
								 p*q*ionfile->vtec_map[i].tec_values[j - 1][k];
							/* if (vtime.hour == 1 && vtime.minute == 59 && vtime.second > 0.0){
								 cout << "1:" << "i:" << i << "j:" << j << "k:" << k << "lat:" << lat_geo << "lon:" << lon_geo<< endl;
								 cout << ionfile->vtec_map[i - 1].tec_values[j][k - 1] <<" "<< ionfile->vtec_map[i - 1].tec_values[j][k]
									 << " " << ionfile->vtec_map[i - 1].tec_values[j - 1][k - 1] << " " << ionfile->vtec_map[i - 1].tec_values[j - 1][k] << " " << v1 << " " << v2<<endl;
								 cout << ionfile->vtec_map[i ].tec_values[j][k - 1] << " " << ionfile->vtec_map[i ].tec_values[j][k]
									 << " " << ionfile->vtec_map[i ].tec_values[j - 1][k - 1] << " " << ionfile->vtec_map[i ].tec_values[j - 1][k] << " " << v1 << " " << v2 << endl;
							 }
							 if (vtime.hour == 2 && vtime.minute == 0 && vtime.second <30.0){
								 cout << "2:" << "i:" << i << "j:" << j << "k:" << k<< "lat:" << lat_geo << "lon:" << lon_geo << endl;
								 cout << ionfile->vtec_map[i - 1].tec_values[j][k - 1] << " " << ionfile->vtec_map[i - 1].tec_values[j][k]
									 << " " << ionfile->vtec_map[i - 1].tec_values[j - 1][k - 1] << " " << ionfile->vtec_map[i - 1].tec_values[j - 1][k] << " " << v1 << " " << v2 << endl;
								 cout << ionfile->vtec_map[i].tec_values[j][k - 1] << " " << ionfile->vtec_map[i].tec_values[j][k]
									 << " " << ionfile->vtec_map[i].tec_values[j - 1][k - 1] << " " << ionfile->vtec_map[i].tec_values[j - 1][k] << " " << v1 << " " << v2 << endl;
							 }*/
							 //时间插值
							 v = (fabs(deltjulianday(&vtime, &ionfile->vtec_map[i - 1].vtime))*v2 +
								 fabs(deltjulianday(&vtime, &ionfile->vtec_map[i].vtime))*v1) / fabs(deltjulianday(&ionfile->vtec_map[i].vtime, &ionfile->vtec_map[i-1].vtime));
							 return v;
						 }
					 }
				 }
			 }
		 }
	 }
 
 }
 //寻找chazh星历
 bool find_sp3_ephem(string prn, ptc ut, psp3 sp3file, s_sp3_ephe sp3[])//n：拉格朗日插值法阶数，sp3：返回的差值节点值
 {
	 int num;//差值点前面的节点数
	 if (NN % 2 == 0)num = NN / 2;//差值阶数为偶数
	 else num = (NN + 1) / 2;
	 JULIANDAY jld1, jld2;//jlld1卫星信号发射时间，jlld2卫星信星历时间
	 /*cout << "卫星信号发射时间：" << ut->year<< ":" << ut->month << ":" << ut->day;
	 cout << ":" << ut->hour << ":" << ut->minute << ":" << ut->second << endl; */
	 int kk = 32;
	 for (int i = 0; i < kk; i++)
	 {
		 if (prn.substr(0, 3) == sp3file->all_ephem[i].prn.substr(0, 3))
		 {
			 //cout << "找到对应的卫星号" << endl;
			 // satn = true;
			 int w = sp3file->all_ephem[i].sate_ephem.size();
			 for (int j = 0; j < w; j++)
			 {
				 double delta;
				 delta = deltjulianday(ut, &sp3file->all_ephem[i].sate_ephem[j].utime_n);
				 if (delta > 0.0) continue;//找到差值点前面的差指点时间
				 else{
					 //cout << "找到历元" << endl;
					 if (j == 0) return false;//插指点在所有插值节点之外（之前）
					 else if (j  < num){//j刚好是插值时间点前面的可用插值点数量，从第一个插值点开始取前n+1个节点用来插值
						 for (int k = 0; k < NN + 1; k++)
						 {
							 sp3[k] = sp3file->all_ephem[i].sate_ephem[k];
							// cout << "前" << "x:" << sp3[k].x << " y:" << sp3[k].y << " z:" << sp3[k].z << endl;
						 }
					 }
					 else if (w - j < NN + 1 - num){//插值点在最后面，后面节点数不足够
						 for (int k = 0; k < NN + 1; k++)
						 {
							 sp3[k] = sp3file->all_ephem[i].sate_ephem[w-1 - NN + k];
							 //cout << "后" << "x:" << sp3[k].x << " y:" << sp3[k].y << " z:" << sp3[k].z << endl;
						 }
					 }
					 else{
						 for (int k = 0; k < NN + 1; k++)
						 {
							 sp3[k] = sp3file->all_ephem[i].sate_ephem[j - num + k];//j是第num+1个
							// cout << "中" << "x:" << sp3[k].x << " y:" << sp3[k].y << " z:" << sp3[k].z << endl;
						 }
					 }
					 return true;
				 }
			 }
		 }
	 }
	 // cout
	 return false;
 }
 
 //精密星历计算卫星坐标
 void cal_sp3_sate_coor(string prn, ptc ut, psp3 sp3file, bool &flag, pxyz coor){
	 coor->x = 0.0;
	 coor->y = 0.0;
	 coor->z = 0.0;
	 s_sp3_ephe sp3[NN + 1];
	 double s = 1.0;
	 if(find_sp3_ephem(prn, ut, sp3file, sp3))
	 {
		 //cout <<"找到星历"<< endl;
		 flag = true;
		 for (int i = 0; i < NN + 1; i++)
		 {//拉格朗日法求解卫星坐标
			 s = 1.0;
			 for (int j = 0; j < NN + 1; j++)
			 {
				 if (i == j)continue;
				 //cout << "x-xj:" << deltjulianday(ut, &sp3[j].utime_n) << "xi-xj:" << deltjulianday(&sp3[i].utime_n, &sp3[j].utime_n) << endl;
				 s = s*deltjulianday(ut, &sp3[j].utime_n) / deltjulianday(&sp3[i].utime_n, &sp3[j].utime_n);
			 }
			 coor->x += s*sp3[i].x;
			 coor->y += s*sp3[i].y;
			 coor->z += s*sp3[i].z;
		 }
	 }
	 else{ 
		 flag = false; 
	 }

 }
 //从dcb文件中找到对应的dcb值
 bool find_dcb(bool b,string prn, string station, pdcb dcbfile, double* dcb_val){//b==0计算卫星dcb，b!=0计算测站dcb
	 //cout << "dcb文件大小：" << dcbfile->dcb_val.size() << endl;
	 for (int i = 0; i < dcbfile->dcb_val.size(); i++)
	 {
		 //cout << "卫星prn:" << dcbfile->dcb_val[i].prn << "  测站:" << dcbfile->dcb_val[i].station << "  obs1:" << dcbfile->dcb_val[i].obs1 << "  obs2:" << dcbfile->dcb_val[i].obs2 << endl;
		 //cout << "prn:" << prn << "  station:" << station <<endl;
		 if (((b==0&&dcbfile->dcb_val[i].prn == prn)||(b!=0&&dcbfile->dcb_val[i].station == station))&&dcbfile->dcb_val[i].obs1 == "C1C"&&dcbfile->dcb_val[i].obs2 == "C2W")
		 {
			 *dcb_val = dcbfile->dcb_val[i].bias*1e-9*c;//返回距离值
			 return true;
		 }
	 }
	 return false;
 }
 /*bool find_dcb(bool b, string prn, string station, pio ionfile, double* dcb_val){//b==0计算卫星dcb，b!=0计算测站dcb
	 //cout << "卫星dcb文件大小：" << ionfile->sate_dcb.size() << endl;
	 //cout << "测站dcb文件大小：" << ionfile->sta_dcb.size() << endl;
	 //cout << "prn:" << prn << endl;
	 if (b == 0)
	 {
		 for (int i = 0; i < ionfile->sate_dcb.size(); i++)
		 {
			 //cout << "prn:" << prn << endl;
			// cout << "卫星prn:" << ionfile->sate_dcb[i].prn << endl;
			 if (ionfile->sate_dcb[i].prn == prn)
			 {
				 //cout << "卫星dcb(ns):" << ionfile->sate_dcb[i].bias << endl;
				 *dcb_val = ionfile->sate_dcb[i].bias*1e-9*c;//返回距离值
				// cout << "卫星dcb(m):" << *dcb_val << endl;
				 return true;
			 }
		 }
	 }
	 else
	 {
		 for (int i = 0; i < ionfile->sta_dcb.size(); i++)
		 {
			 //cout << "prn:" << prn << endl;
			 //cout << "测站:" << ionfile->sta_dcb[i].station << endl;
			 if (ionfile->sta_dcb[i].station == station)
			 {
				 //cout << "测站dcb(ns):" << ionfile->sta_dcb[i].bias << endl;
				 *dcb_val = ionfile->sta_dcb[i].bias*1e-9*c;//返回距离值
				 //cout << "测站dcb(m):" << *dcb_val << endl;
				 return true;
			 }
		 }
	 }
	 return false;
 }*/
 //结果输出
 void putresult(pvt pt)
 {
	 ofstream outfile("E:\\GNSS\\结果\\vtec.txt", ios::out);
	 if (!outfile)
	 {
		 cerr << "文件创建失败" << endl;
		 abort();
	 }
	 cout << "开始写文件" << endl;
	 outfile << "  年 ,月,日,时,分,    秒,      穿刺点纬度     穿刺点经度        vtec1           VTEC         ionex_vtec       d_ion      mf         prn"<<endl;
	 for (unsigned int i = 0; i < pt->allresult.size(); i++)
	 {
		 outfile.setf(ios_base::fixed, ios_base::floatfield);
		 outfile.precision(4);
		 outfile << setw(4) << pt->allresult[i].rtime.year << setw(3) << pt->allresult[i].rtime.month << setw(3) << pt->allresult[i].rtime.day;
		 outfile << setw(3) << pt->allresult[i].rtime.hour << setw(3) << pt->allresult[i].rtime.minute << setw(9) << pt->allresult[i].rtime.second;
		 outfile << std::right << setw(15) << pt->allresult[i].lat*180.0 / PI << std::right << setw(15) << pt->allresult[i].lon*180.0 / PI << std::right << setw(15) << pt->allresult[i].vtec1 << std::right << setw(15) << pt->allresult[i].vtec
			 << std::right << setw(15) << pt->allresult[i].ionex_vtec / 10.0 << std::right << setw(15) << pt->allresult[i].vtec - pt->allresult[i].ionex_vtec / 10.0 << setw(10) << pt->allresult[i].mf << setw(10) << pt->allresult[i].prn << endl;
	 }
	 outfile.close();
	 cout << "文件写入完毕" << endl;
 }
 //计算VTEC
 void vtec(pobs obsfile, psp3 sp3file, pio ionfile, pdcb dcbfile,pvt pt){
	// cout <<"开始计算" << endl;
	 int num=0,qq=0;//num:穿刺点数量,qq:连续观测一颗卫星的观测历元数,cycle_slip:周跳个数
	 XYZ sta_coor;//测站wgs84坐标
	 BLH sta_coor_blh,ipp;//测站和穿刺点的大地坐标
	 pxyz sate_coor=new XYZ;//卫星坐标
	 penupolar pep=new ENUPOLAR;//卫星方位角和高度角
	 bool flag;
	 double* mf=new double;//投影函数
	 JULIANDAY ju1, ju2;
	 UTC ut;
	 one_result re;
	 double dt1,dt2=0.0;//dt1：卫星信号传播时间,dt2：同一颗卫星相邻历元的时间差
	 double w=1.0;//权
	 vector<dcb> dcb_efficient;//方程组的所有系数和右端项
	 dcb one_equ;//一个方程的系数和右端项
	 double dphi, ds;//纬度差和时角差
	 double mf0=1.0;//上一步的mf
	 double dcb_sat = 0.0, dcb_sta = 0.0;//测站和卫星dcb
	 MatrixXd  Q(5, 5), V(5, 1);//5个未知数
	 //周跳探测
	 double n_mw1 = 0.0, n_mw2 = 0.0;//MW组合
	 double n0 = 0.0, n1 = 0.0;//均值
	 double sigma0 = 0.0, sigma1 = 0.0;//标准差
	 sta_coor = obsfile->obsheaddata.approx_coordinate;
	 xyztoblh(&sta_coor, &sta_coor_blh);
	 cout << "测站地理坐标：（" << sta_coor_blh.latitude*180.0 / PI << "," << sta_coor_blh.longitude*180.0 / PI << ")" << endl;
	 //cout << "测站纬度;" << sta_coor_blh.latitude*180.0/PI << "测站经度：" << sta_coor_blh.longitude*180.0/PI << endl;
	 for ( int i = 0; i < 32; i++)//所有卫星依次使用相位平滑伪距的方法求vtec
	 {
		 re.prn = obsfile->all_sate[i].prn;
		 int cycle_slip = 0;
		 //cout << "i:" << i << endl;
		 cout << i + 1 << "卫星历元数:" << obsfile->all_sate[i].obs_data.size()<<endl;
		 if (!find_dcb(0, obsfile->all_sate[i].prn, obsfile->obsheaddata.station, dcbfile, &dcb_sat) || !find_dcb(1, obsfile->all_sate[i].prn, obsfile->obsheaddata.station, dcbfile, &dcb_sta))
		 {
			 cout << "没有找到卫星" << obsfile->all_sate[i].prn << "或者测站dcb" << endl;
			 continue;
		 }
		  // cout << "卫星dcb：" << dcb_sat*1.0e9 / c << "  测站dcb：" << dcb_sta*1.0e9 / c << endl;
		 //cout << "卫星dcb：" << dcb_sat << "  测站dcb：" << dcb_sta << endl;
		//初始化
		 qq = 0;
		 n_mw1 = 0.0, n_mw2 = 0.0;//MW组合
		 n0=0.0, n1=0.0;//均值
		 sigma0=0.0, sigma1=0.0;//标准差
		 for (int j = 0; j < obsfile->all_sate[i].obs_data.size(); j++)//一颗卫星所有历元按照相位平滑伪距的方法求vtec
		 {
			// cout << "prn:" << obsfile->all_sate[i].prn << " 时间：" << ut.year << " " << ut.month << " " << ut.day << " " << ut.hour << " " << ut.minute << " " << ut.second << endl;
			 cal_sp3_sate_coor(obsfile->all_sate[i].prn, &obsfile->all_sate[i].obs_data[j].utime_o, sp3file, flag, sate_coor);//求卫星坐标
			 if (!flag){
				 continue;
			 }//卫星坐标计算失败
			 sate_azi_ele(&sta_coor, sate_coor, pep);//求卫星方位角和高度角
			 if (pep->elevation*180.0 / PI < 15.0)continue;//卫星高度角要大于15度
			 ipp_pos(&sta_coor_blh, pep, &ipp, mf);//求穿刺点纬度经度
			 re.rtime = obsfile->all_sate[i].obs_data[j].utime_o;
			/* if (obsfile->all_sate[i].obs_data[j].utime_o.hour == 8 && obsfile->all_sate[i].obs_data[j].utime_o.minute == 1){
				 double psi = PI / 2.0 - pep->elevation - asin(ave_a / (ave_a + hion)*cos(pep->elevation)), la = asin(sin(sta_coor_blh.latitude)*cos(psi) + cos(sta_coor_blh.latitude)*sin(psi)*cos(pep->azimuth));
				 penu pe;
				 pe = (penu)malloc(sizeof(ENU));
				 xyztoenu(&sta_coor, sate_coor, pe);
				 enutoenupolar(pe, pep);
				 xyztoblh(&sta_coor, &sta_coor_blh);
				 double dx, dy, dz;
				 dx = sate_coor->x - sta_coor.x;
				 dy = sate_coor->y - sta_coor.y;
				 dz = sate_coor->z - sta_coor.z;
				 cout.setf(ios::fixed);
				 cout << "大地纬度：" << sta_coor_blh.latitude << endl;
				 cout << "大地经度：" << sta_coor_blh.longitude << endl;
				 cout << "大地北卫星x：" << -sin(sta_coor_blh.latitude)*cos(sta_coor_blh.longitude)*dx << endl;
				 cout << "大地北卫星y：" << -sin(sta_coor_blh.latitude)*sin(sta_coor_blh.longitude)*dy << endl;
				 cout << "大地北卫星z：" << cos(sta_coor_blh.latitude)*dz << endl;
				 cout << "卫星x：" << sate_coor->x << endl;
				 cout << "卫星y：" << sate_coor->y << endl;
				 cout << "卫星z：" << sate_coor->z << endl;
				 cout << "dx：" << sate_coor->x - sta_coor.x << endl;
				 cout << "dy：" << sate_coor->y - sta_coor.y << endl;
				 cout << "dz：" << sate_coor->z - sta_coor.z << endl;
				 cout << "卫星北：" << pe->northing  << endl;
				 cout << "卫星东：" << pe->easting << endl;
				 cout << "卫星天：" << pe->upping  << endl;
				 cout << "比值：" << pe->easting / pe->northing <<"方位角"<< atan(pe->easting / pe->northing) << endl;
				 free(pe);
				 cout << "卫星方位角：" << pep->azimuth*180.0 / PI << "卫星高度角：" << pep->elevation*180.0 / PI << endl;
				 cout << "张角：" << psi*180.0 / PI << endl;
				 cout << "赤纬：" << la*180.0 / PI << endl;
				 cout << "赤经：" << sta_coor_blh.longitude*180.0 / PI + atan(cos(sta_coor_blh.latitude)*sin(psi)*sin(pep->azimuth) / (cos(psi) - sin(sta_coor_blh.latitude)*sin(la)))*180.0 / PI << endl;
			 }*/
			 re.lat = ipp.latitude;
			 re.lon = ipp.longitude;
			 re.mf = *mf;
			 re.vtec1 = 9.52437*(obsfile->all_sate[i].obs_data[j].p2 - obsfile->all_sate[i].obs_data[j].p1 + dcb_sta + dcb_sat) / (*mf);
			 re.ionex_vtec = ionex_vtec(ionfile, re.rtime, re.lat*180.0 / PI, re.lon*180.0 / PI);
			 if (j != 0&&num>0){//第一步不计算dt2
				 dt2 = fabs(deltjulianday(&obsfile->all_sate[i].obs_data[j].utime_o, &pt->allresult[num - 1].rtime));
			 }
			 if (qq==0||num==0||dt2 > 300.0)//qq=0意味着上一步发现周跳或者第一步，相邻历元时间间隔大于300s，证明卫星落下又升起，需要选为时间起点计算vtec
			 {
				 /*cout << "卫星：" << re.prn << endl;
				 cout << " 时间：" << obsfile->all_sate[i].obs_data[j].utime_o.hour << " " << obsfile->all_sate[i].obs_data[j].utime_o.minute << " " << obsfile->all_sate[i].obs_data[j].utime_o.second << endl;
				 cout << "重新开始计算"<<endl;
				 */
				 re.vtec = 9.52437*(obsfile->all_sate[i].obs_data[j].p2 - obsfile->all_sate[i].obs_data[j].p1+dcb_sta+dcb_sat)/ (*mf);
				 w = 1.0;
				 n_mw1 = obsfile->all_sate[i].obs_data[j].l1 - obsfile->all_sate[i].obs_data[j].l2 - (f1 - f2) / (f1 + f2)*(obsfile->all_sate[i].obs_data[j].p1*f1 / c + obsfile->all_sate[i].obs_data[j].p2*f2 / c);
				 qq = 1;
				 n1 = n_mw1;
				 sigma1 = 0.0;
				// qq = 0;  
			 }
			 else{
				 n0 = n1;
				 sigma0 = sigma1;
				 n_mw1 = obsfile->all_sate[i].obs_data[j].l1 - obsfile->all_sate[i].obs_data[j].l2 - (f1 - f2) / (f1 + f2)*(obsfile->all_sate[i].obs_data[j].p1*f1 / c + obsfile->all_sate[i].obs_data[j].p2*f2 / c);
				 if (j<obsfile->all_sate[i].obs_data.size() - 1){ n_mw2 = obsfile->all_sate[i].obs_data[j+1].l1 - obsfile->all_sate[i].obs_data[j+1].l2 - (f1 - f2) / (f1 + f2)*(obsfile->all_sate[i].obs_data[j+1].p1*f1 / c + obsfile->all_sate[i].obs_data[j+1].p2*f2 / c); }
				 n1 = n0 + 1.0 / qq*(n_mw1-n0);
				 sigma1 = sqrt(sigma0*sigma0 + 1.0 / qq*((n_mw1 - n0)*(n_mw1 - n0) - sigma0*sigma0));
				 if (qq>2&&fabs(n_mw1 - n0) >= 4.0*sigma0&&fabs(n_mw2 - n_mw1) < 1.0){//发现周跳，最后一个历元n_mw2 == n_mw1，因为没有重新赋值给n_mw2
					/* cout << "qq:" << qq << endl;
					 cout << "sigma0:" << sigma0 << endl;
					 cout << "n_mw1:" << n_mw1<< endl;
					 cout << "n0:" << n0 << endl;
					 cout << "fabs(n_mw1 - n0):" << fabs(n_mw1 - n0) << endl;
					 cout << "n_mw2:" << n_mw2 << endl;
					 cout << "发现周跳" << endl;
					 cout << " 时间：" << obsfile->all_sate[i].obs_data[j].utime_o.hour << " " << obsfile->all_sate[i].obs_data[j].utime_o.minute << " " << obsfile->all_sate[i].obs_data[j].utime_o.second << endl;
					 */
					 /*
					 cout << "卫星：" << re.prn << endl;
					 cout << " 时间：" << obsfile->all_sate[i].obs_data[j].utime_o.hour << " " << obsfile->all_sate[i].obs_data[j].utime_o.minute << " " << obsfile->all_sate[i].obs_data[j].utime_o.second << endl;
					 cout << "发现周跳" << endl;
					*/
					 cycle_slip++;
					 qq = 0;
					 continue;
				 }
				 w = 1.0 / qq;
				 //w = 1.0 - 0.01*qq;
				 //if (w < 0.01)w = 0.01;
				/* if (re.rtime.hour == 9 && re.rtime.minute == 26 && re.prn == "G07"){
					 cout.setf(ios::fixed);
					 cout << "观测值：" << 9.52437*(obsfile->all_sate[i].obs_data[j].p2 - obsfile->all_sate[i].obs_data[j].p1 + dcb_sta + dcb_sat)*w << endl;
					 cout << "平滑值：" << 9.52437*(1.0 - w)*(re.vtec / 9.52437*mf0 +
						 ((obsfile->all_sate[i].obs_data[j].l1*c / f1 - obsfile->all_sate[i].obs_data[j].l2*c / f2) - (obsfile->all_sate[i].obs_data[j - 1].l1*c / f1 - obsfile->all_sate[i].obs_data[j - 1].l2*c / f2))) << endl;
					 cout << "平滑值后一步：" << 9.52437*( obsfile->all_sate[i].obs_data[j].l1*c / f1 - obsfile->all_sate[i].obs_data[j].l2*c / f2 )<< endl;
					 cout << "平滑值前一步：" << 9.52437*(obsfile->all_sate[i].obs_data[j - 1].l1*c / f1 - obsfile->all_sate[i].obs_data[j - 1].l2*c / f2) << endl;
					 cout << "权w：" << w << endl;
					 cout << "相位：" << obsfile->all_sate[i].obs_data[j].l1 << " " << obsfile->all_sate[i].obs_data[j].l2 << " " << obsfile->all_sate[i].obs_data[j - 1].l1 << " " << obsfile->all_sate[i].obs_data[j - 1].l2 << endl;
				 }*/
				 re.vtec = 9.52437*((obsfile->all_sate[i].obs_data[j].p2 - obsfile->all_sate[i].obs_data[j].p1 + dcb_sta + dcb_sat)*w + (1.0 - w)*(re.vtec / 9.52437*mf0 +
					 ((obsfile->all_sate[i].obs_data[j].l1*c / f1 - obsfile->all_sate[i].obs_data[j].l2*c / f2) - (obsfile->all_sate[i].obs_data[j - 1].l1*c / f1 - obsfile->all_sate[i].obs_data[j - 1].l2*c / f2)))) / (*mf);
			 }
			 /*cout.setf(ios::fixed);
			 cout << " 时间：" << obsfile->all_sate[i].obs_data[j].utime_o.hour << " " << obsfile->all_sate[i].obs_data[j].utime_o.minute << " " << obsfile->all_sate[i].obs_data[j].utime_o.second << endl;
			 cout << "p2-p1："<<obsfile->all_sate[i].obs_data[j].p2 - obsfile->all_sate[i].obs_data[j].p1 << endl;
			 cout << "系数："<<9.52437 / (*mf) << endl;*/
			 if (re.vtec<0.0){
				 qq = 0;
				 continue;
			 }
			 pt->allresult.push_back(re);
			 num++;
			 qq++;
			 mf0 = *mf;
		 }
		 //cout << "周跳个数:" << cycle_slip << endl;
	 }
	 delete sate_coor;
	 delete pep;
	 delete mf;
	 cout << "计算完毕" << endl;
	 cout << "结果个数：" << pt->allresult.size() << endl;
 }