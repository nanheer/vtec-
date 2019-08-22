#ifndef FUNCTION_H
#define FUNCTION_H
#include<string>
#include"�ṹ��.h"
#include"����.h"
//#include"matrix.h"
using namespace std;
//��sp3�ļ�
void readsp3file(string strn, psp3 sp3file);
//����vtec����ķ�ʽ��ȡo�ļ�
void readofile_vtec(string stro, pobs obsfile);
//��ȡionex�ļ�
void read_ionex(string stri, pio ionfile);
//��ȡdcb�ļ�
void read_dcb(string strd, pdcb dcbfile);
//GPST��UTC
void gpsttoutc(pgpst gt, ptc ut);
//JULIANDAY��BDT
void juliandaytobdt(pjulian ju, pbdt bt);
//UTC��BDT
void utctobdt(ptc ut, pbdt bt);
//UTC��GPST
void utctogpst(ptc ut,pgpst gt);

//UTC��JULIANDAY
void utctojulianday(ptc ut,pjulian ju);

//JULIANDAY��UTC
void juliandaytoutc(pjulian ju,ptc ut);

//GPST��JULIANDAY
void gpsttojulianday(pgpst gt, pjulian ju);

//JULIANDAY��GPST
void juliandaytogpst(pjulian ju, pgpst gt);
//�������ղ�ֵ
double deltjulianday(ptc u1, ptc u2);
//�����ձ任
void transjulian(pjulian ju1, double * dt, pjulian ju2);
//

//XYZ��BLH
void xyztoblh(pxyz px,pblh pb);
//BLH��XYZ
void blhtoxyz(pblh pb, pxyz px);
//XYZ��ENU
void xyztoenu(pxyz pxcenter, pxyz px, penu pe);
//ENU��XZY
void enutoxyz(pxyz pxcenter,penu pe, pxyz px);
//ENU��ENUPOLAR
void enutoenupolar(penu pe, penupolar pep);
//�����Ƿ�λ�Ǻ͸߶Ƚ�
void sate_azi_ele(pxyz pxcenter, pxyz px, penupolar pep);
//�������ӳٸ���
void tropo(double mjd,pxyz px1, pxyz px2,double * dx);
//����������
void get_nominal_metedata(double mjd, double lat, double lon,double dhgt, double* pres, double* temp, double* rhumi, double* undu);
//���Сֵ
int min(int a, int b);
//���ջ����߸߸���
void antena_height_correction(pxyz px1, pxyz px2, penu pe, double*p1, double*p2, double*pp1, double*pp2);
//��������������������
void cal_sate_coor(string prn, ptc ut, psp3 sp3file, bool &flag, pxyz coor);
//��ƽ����Ǽ���ƫ�����
double calE(double* M, double* e);
//��ƫ�������������
double calf(double* E, double* e);
//��ת����
//void rotamatrix(string str, double seta,Matrix r);
//������
void putresult(pvt pt);
//vtec����
void vtec(pobs obsfile, psp3 sp3file, pio ionfile, pdcb dcbfile, pvt pt);
//��ionex�ļ���ֵ����������vtec
double ionex_vtec(pio ionfile, UTC vtime, double lat, double lon);//lat��lon�Ǵ�ؾ�γ��
//��dcb�ļ����ҵ���Ӧ��dcbֵ
bool find_dcb(bool b, string prn, string station, pdcb dcbfile, double* dcb_val);
//bool find_dcb(bool b, string prn, string station, pio ionfile, double* dcb_val);
#endif