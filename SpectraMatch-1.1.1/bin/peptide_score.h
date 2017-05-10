#include<bits/stdc++.h>
using namespace std;
#define maxn 1127970 //库的大小
#define minn 1e-8
#define threshold 7.0
#define MAX 10000010

#ifndef SCORE_H
#define SCORE_H

double vectorsum(vector<double>::iterator start, vector<double>::iterator start1, vector<double>::iterator end1, vector<double>::iterator endl1);
double vectorSum(vector<double>::iterator start, vector<double>::iterator end1);
double vector_sqrt_Sum(vector<double>::iterator start, vector<double>::iterator end1);
bool compare(pair<double, double> p1, pair<double, double> p2);
bool compare1(pair<double, double> p1, pair<double, double> p2);
void Print(char* filename3, char* filename4, char* filename5);
int BinarySearchRecursive_left(int low, int high, double key);
int BinarySearchRecursive_right(int low, int high, double key);
string *FileList(char* filename);
string peptide_name(string s);

struct p
{
	char s[300];
	double mu;
};
extern p ans[maxn];
extern p ans_sort[maxn];
int cmp(p a, p b);

struct sxp
{
	double score;
	string s;
};
extern sxp list_score[MAX];
int cmp_sxp(sxp a,sxp b);
void mion_sort(char* filename1,char* filename2);
void sorted_lib_pathname2mw(char* filename);
double select(char line[]);
void* thread1(void *args);
void* thread2(void *args);
void* thread3(void *args);
void* thread4(void *args);
void* thread5(void *args);
void* thread6(void *args);
void* thread7(void *args);
void* thread8(void *args);
void* thread9(void *args);
void* thread10(void *args);
void* thread11(void *args);
void* thread12(void *args);
void* thread13(void *args);
void* thread14(void *args);
void* thread15(void *args);
void* thread16(void *args);
void* thread17(void *args);
void* thread18(void *args);
void* thread19(void *args);
void* thread20(void *args);
int score(char* Filename,int n,int m,char* file_name);

#endif   
