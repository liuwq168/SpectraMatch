#include<bits/stdc++.h>
#include "peptide_score.h"
using namespace std;
char Name[200];
char Thread1_name[200];
char Thread2_name[200];
char Thread3_name[200];
char Thread4_name[200];
char Thread5_name[200];
char Thread6_name[200];
char Thread7_name[200];
char Thread8_name[200];
char Thread9_name[200];
char Thread10_name[200];
char Thread11_name[200];
char Thread12_name[200];
char Thread13_name[200];
char Thread14_name[200];
char Thread15_name[200];
char Thread16_name[200];
char Thread17_name[200];
char Thread18_name[200];
char Thread19_name[200];
char Thread20_name[200];
int N;
int P;
const int StringNum = 200000;

namespace DenoiseCutThreshold{
struct Pos
{
	double x, y;
}arr[100010], sub[100010], ans[100010], last;

double mu = 0.0;

int CmpY(Pos a, Pos b)
{
	if (a.y == b.y)
		return a.x < b.x;
	return a.y > b.y;
}

int CmpX(Pos a, Pos b)
{
	if (a.x == b.x)
		return a.y > b.y;
	return a.x < b.x;
}

double mark[15];
int temp = 0, cnt = 0, i = 0, j = 0, t = 0, cnt1 = 0;
const int range = 100.0;

void DenoiseByCutThreshold(char* BeforeDenoiseFilename, char* AfterDenoiseFilename)
{
	
	ifstream input;
	input.open(BeforeDenoiseFilename);
	input.precision(7);
	input.setf(ios_base::fixed, ios_base::floatfield);
	
	ofstream output;
	output.open(AfterDenoiseFilename);
	output.precision(7);
	output.setf(ios_base::fixed, ios_base::floatfield);
	
	double a;
	cnt = 0;
	memset(arr, 0, sizeof(arr));
	memset(sub, 0, sizeof(sub));
	memset(mark, 0, sizeof(mark));
	temp = 0;
	while (input >> a)
	{
		if (!temp)arr[cnt].x = a;
		else arr[cnt++].y = a;
		temp = !temp;
	}
	mu = arr[cnt].x;
	int anscnt = int(0.5 * (cnt-1));
	cnt -= 1;
	last.x = arr[cnt].x, last.y = arr[cnt].y;
	
	sort(arr, arr + cnt, CmpY);
	sort(arr, arr + anscnt, CmpX);
	
	for (i = 0; i<anscnt; i++)
	{
		output << arr[i].x << " " << arr[i].y << endl;
	}
	output << last.x << " " << last.y << endl;
	
	input.close();
	output.close();
	return;
}
}
namespace ZoneDenoise{
const int MAX_POINT = 1e4 + 10;
const int RANGE = 1e2;
const double DeleteRange = 1.01;
typedef pair<double, double> pdd;
vector<pdd> point, answer;
vector<string> inString;
bool tabu[MAX_POINT];
void init_point() {
	point.clear();
	char temp[MAX_POINT];
	for (int i = 0; i<inString.size() - 1; i++) {
		double x, y;
		strcpy(temp, inString[i].c_str());
		sscanf(temp, "%lf%lf", &x, &y);
		point.push_back(make_pair(x, y));
	}
}
double Distance(int x, int y) {
	double d = point[x].first - point[y].first;
	if (d<0) d = -d;
	return d;
}
void get(int Start, int End) {
	int Count = End - Start + 1;
	int need = max(10, Count / 2);
	while (need>0) {
		int high = -1;
		for (int i = Start; i <= End; i++) {
			if (tabu[i]) continue;
			if (high == -1) {
				high = i;
				continue;
			}
			if (point[high].second<point[i].second) {
				high = i;
			}
		}
		if (high == -1) return;
		for (int i = Start; i <= End; i++) {
			if (tabu[i]) continue;
			if (Distance(high, i)<DeleteRange) {
				if (point[high].second < point[i].second) {
					high = i;
					tabu[i] = true;
				}
			}
		}
		answer.push_back(point[high]);
		tabu[high] = true;
		need--;
	}
}
void solve() {
	init_point();
	sort(point.begin(), point.end());
	for (int i = 0; i<MAX_POINT; i++) {
		tabu[i] = false;
	}
	answer.clear();
	for (int x = point[0].first; x <= point[point.size() - 1].first; x += RANGE) {
		for (int i = 0; i<point.size(); i++) {
			if (point[i].first<x) continue;
			if (point[i].first>x + RANGE) break;
			int last = i;
			for (int j = i; j<point.size(); j++) {
				if (point[j].first>x + RANGE) break;
				last = j;
			}
			get(i, last);
			break;
		}
	}
	sort(answer.begin(), answer.end());
}

void DenoiseByZone(char* filename1, char* filename2) {
	ifstream input;
	input.open(filename1);
	input.precision(10);
	input.setf(ios_base::fixed, ios_base::floatfield);
	inString.clear();
	string bufferString;
	while (getline(input, bufferString)) {
		inString.push_back(bufferString);
	}
	input.close();
	solve();
	ofstream output;
	output.open(filename2);
	output.precision(10);
	output.setf(ios_base::fixed, ios_base::floatfield);
	for (int i = 0; i<answer.size(); i++) {
		output << answer[i].first << " " << answer[i].second << endl;
	}
	//output << answer[answer.size() - 1].first << " " << answer[answer.size() - 1].second << endl;
	output << inString[inString.size()-1] << endl;
	output.close();
}
}
namespace ZoneChargeDenoise{
	const int MAX_POINT = 1e4 + 10;
const int RANGE = 1e2;
const double DeleteRange = 1.01;
double Maximum=0.0;
typedef pair<double, double> pdd;
vector<pdd> point, answer;
vector<string> inString;
bool tabu[MAX_POINT];
double chg;
void init_point() {
	point.clear();
	char temp[MAX_POINT];
	for (int i = 0; i<inString.size(); i++) {
		double x, y;
		strcpy(temp, inString[i].c_str());
		sscanf(temp, "%lf%lf", &x, &y);
		point.push_back(make_pair(x, y));
	}
	vector<double> rep;
	for (int i = 0; i < point.size(); i++)
	{
		rep.push_back(point[i].second);
	}
	chg=rep[point.size()-1];
	
	sort(rep.begin(), rep.end());
	Maximum=rep[rep.size()-1];
	
	vector<double>().swap(rep);
}
int Distance(int x, int y) {
	double d = point[x].first - point[y].first;
	if (d<0) d = -d;
	if (d >= 0.99&&d <= 1.01)return 1;
	return 0;
}
int Distance2(int x, int y) {
	double d = point[x].first - point[y].first;
	if (d<0) d = -d;
	if (d >= 0.49&&d <=0.51 )return 1;
	return 0;
}
int Distance3(int x, int y) {
	double d = point[x].first - point[y].first;
	if (d<0) d = -d;
	if (d >= 0.32&&d <= 0.34)return 1;
	return 0;
}
int Distance4(int x, int y) {
	double d = point[x].first - point[y].first;
	if (d<0) d = -d;
	if (d >= 0.24&&d <= 0.26)return 1;
	return 0;
}
void get(int Start, int End) {
	int Count = End - Start + 1;
	int need = max(10, Count / 2);
	int flag = 1;
	while (need>0||flag) {
		int high = -1;
		for (int i = Start; i <= End; i++) {
			if (tabu[i]) continue;
			if (high == -1) {
				high = i;
				continue;
			}
			if (point[high].second<point[i].second) {
				high = i;
			}
		}
		if (high == -1) return;
		for (int i = Start; i <= End; i++) {
			if (tabu[i]) continue;
			if (int(chg) == 1&&Distance(high, i)) {
				tabu[i] = true;
			}
			if (int(chg)>=2&&Distance2(high, i)) {
				tabu[i] = true;
			}
			if (int(chg)>=3&&Distance3(high, i)) {
				tabu[i] = true;
			}
			if (int(chg)>=4&&Distance4(high, i)) {
				tabu[i] = true;
			}

		}
		if (point[high].second < Maximum/5)
		{
			flag = 0;
			if (need <= 0)break;

		}
		answer.push_back(point[high]);
		tabu[high] = true;
		need--;
	}
}
void solve() {
	init_point();
	sort(point.begin(), point.end());
	for (int i = 0; i<point.size(); i++) {
		tabu[i] = false;
	}
	answer.clear();
	for (int x = point[0].first; x <= point[point.size() - 3].first; x += RANGE) {
		for (int i = 0; i<point.size()-2; i++) {
			if (point[i].first<x) continue;
			if (point[i].first>x + RANGE) break;
			int last = i;
			for (int j = i; j<point.size()-2; j++) {
				if (point[j].first>x + RANGE) break;
				last = j;
			}
			get(i, last);
			break;
		}
	}
	sort(answer.begin(), answer.end());
}

void DenoiseByCharge(char* filename1, char* filename2) {
	ifstream input;
	input.open(filename1);
	input.precision(10);
	input.setf(ios_base::fixed, ios_base::floatfield);
	inString.clear();
	string bufferString;
	while (getline(input, bufferString)) {
		inString.push_back(bufferString);
	}
	input.close();
	solve();
	ofstream output;
	output.open(filename2);
	output.precision(10);
	output.setf(ios_base::fixed, ios_base::floatfield);
	for (int i = 0; i<answer.size(); i++) {
		output << answer[i].first << " " << answer[i].second << endl;
	}
	//output << answer[answer.size() - 1].first << " " << answer[answer.size() - 1].second << endl;
	output << inString[inString.size()-1] << endl;
	output.close();
}
}

namespace FiveDotFilter{
const double uu = 0.5;
const int MAX_POINT = 1e4 + 10;
struct myclass {
  bool operator() (double i,double j) { return i+1e-6<j; }
} myobj;
typedef pair<double, double> pdd;
vector<pdd> point;
vector<string> inString;
double Leftside_one_condition(vector<double> part, int n)  
{
	if (part[n]>uu*(*max_element(part.begin()+n-2,part.begin()+n+2,myobj)))
		return part[n];
	else
		return 0;
}

double Leftside_two_condition(vector<double> part, int n) 
{
	sort(part.begin()+n-2,part.begin()+n+2,myobj);
	if (part[n]>uu*(part[n]))
		return part[n];
	else
		return 0;
}

double Rightside_one_condition(vector<double> part, int n)  
{
	if (part[n]<uu*(*max_element(part.begin()+n-2,part.begin()+n+2,myobj)))
		return part[n];
	else
		return 0;
}

double Rightside_two_condition(vector<double> part, int n)  
{
	sort(part.begin()+n-2,part.begin()+n+2,myobj);
	if (part[n]<-uu*(part[n]))
		return part[n];
	else
		return 0;
}

double Zxh_about_section(double test[][2], int begin, int end)  
{
	double zxh = 0.0;
	double t = 0.0;
	double r = 0.0;
	for (int i = begin; i<end + 1; i++)
	{
		t += test[i][0] * test[i][1];
		r += test[i][1];
	}
	zxh = t / r;
	return zxh;  
}

double Max_in_section(int start, int stop)
{
	double Max = 0.0;
	Max = point[start].second;
	for (int i = start + 1; i<stop + 1; i++)
		if (point[i].second > Max)
			Max = point[i].second;
	return Max;
}

int Max_in_index(int start, int stop)
{

	int ncount = start;
	for (int i = start + 1; i<stop + 1; i++)
		if (point[i].second > point[ncount].second)
			ncount = i;
	return ncount;
}

void init_point() {
	point.clear();
	char temp[MAX_POINT];
	for (int i = 0; i<inString.size() - 1; i++) {
		double x, y;
		strcpy(temp, inString[i].c_str());
		sscanf(temp, "%lf%lf", &x, &y);
		point.push_back(make_pair(x, y));
	}
}

void DenoiseByFiveDot(char* filename1, char* filename2)
{
	ifstream input;
	input.open(filename1);
	input.precision(10);
	input.setf(ios_base::fixed, ios_base::floatfield);
	inString.clear();
	string bufferString;
	while (getline(input, bufferString)) {
		inString.push_back(bufferString);
	}
	input.close();
	init_point();
	//cout<<"111"<<endl;
	ofstream outFile2;
	outFile2.open(filename2);
	outFile2.precision(10);
	outFile2.setf(ios_base::fixed, ios_base::floatfield);
	
	vector<double> Cacdf(point.size(),0), Cacfab(point.size(),0);
	vector<double> Cacdf_1(point.size(),0), Cacfab_2(point.size(),0);
	if (point.size() > 10){
		for (int i = 0; i < point.size() - 4; i++)
		{
			Cacdf[i+2]=(-2 * point[i].second - point[i+1].second + point[i+3].second + point[i+4].second / 10);
		}
		for (int i = 2; i < point.size() - 3; i++)
			Cacfab[i]=fabs(point[i+1].second-point[i].second);
		
		int k1, k2, ind;
		int i = 2;
		int j;
		//cout<<"222"<<endl;
		for (; i < point.size() - 6; i++)
		{
			double ls = Leftside_one_condition(Cacfab, i + 2);  
			k1 = i + 2;
			if (ls != 0)
			{
				double ls1 = Leftside_two_condition(Cacdf, k1); 
				if (ls1 != 0)
				{
					j = k1 + 1;
					double rs = Rightside_one_condition(Cacfab, j + 2);  
					k2 = j + 2;									
					if (rs != 0)
					{
						double rs1 = Rightside_two_condition(Cacdf, k2); 
																			
						if (rs1 != 0)
						{
							ind = Max_in_index(k1, k2);
							Cacdf_1[j]=point[ind].first; 
							//cout<<"333"<<endl;					 
							Cacfab_2[j]=Max_in_section(k1, k2);

							if (Cacdf_1[j] != 0 && Cacfab_2[j] != 0)
								outFile2 << Cacdf_1[j] << " " << Cacfab_2[j] << "\n";
							else
								break;
							//cout<<"444"<<endl;
						}
						else
							continue;
					}
					else
						continue;
				}
				else
					continue;
			}
			else
				continue;

			i = k2;
		}
	}
	else{
		for (int i=0;i<inString.size()-1;i++)
			outFile2<< inString[i] << endl;
	}
	outFile2 << inString[inString.size()-1] << endl;
	outFile2.close();


	vector<double>().swap(Cacdf);
	vector<double>().swap(Cacfab);
	vector<double>().swap(Cacdf_1);
	vector<double>().swap(Cacfab_2);
	vector<pair<double, double> >().swap(point);
}
}

int main(int argc,char* argv[])
{
	string usage="\nProgram: SpectraMatch (Tools for scores in the peptide searching library)\n";
	usage+="Version: 1.1.1 \n\n";
	usage+="Usage:   SpectraMatch <command> [options]\n";
	usage+="Commands:\n";
	usage+="  --Filtering\n";
	usage+="	thrs     \tusing one-cut threshold to filter the spectra noise\n";
	usage+="	sect     \tselect one section to filter the spectra noise\n";
	usage+="	chrg     \tbase on the charge state to filter the spectra noise\n";
	usage+="	five     \tbase on the five dot to filter the spectra noise\n\n";
	usage+="  --SpectraScore\n";
	usage+="	sortlib  \tconvert the library pathname to sorted mass ion table list\n";
	usage+="	score    \tusing the filtered spectra do the SpectraMatch process\n";
	
	string pep_version = "SpectraMatch ";
#ifdef PEP_VERSION  
	pep_version += PEP_VERSION;
#else
	pep_version += "v1.1.1";
#endif
	pep_version += "\n";	
	
	cin.seekg(0, ios::end);  //表示输入流的结束位置
	int cin_length = cin.tellg();  //tellg（）函数不需要带参数，它返回当前定位指针的位置，也代表着输入流的大小
	cin.seekg(0, ios::beg);  //表示输入流的开始位置

	if (argc == 2) {
		string tmp(argv[1]);  //定义string类型tmp，并赋初值argv[1],值为AGE v0.4
		if (tmp == "-version") {  //  ./age_align -version  输出位  AGE v0.4
			cout << pep_version << endl;
			return 0;
		}
		if (tmp == "thrs"){
			string show_thrs="Usage: SpectraMatch thrs [options] in.txt\n";
			show_thrs+="Input options:\n";
			show_thrs+="  --in-file-list \tthe set of spectra with noise, input the spectra path-name list\n\n";
			show_thrs+="Output options:\n";
			show_thrs+="  --out-file-list\tthe set of spectra filtered, mean no noise spectra, output the spectra path-name list\n\n";
			show_thrs+="ParameterSet:\n";
			show_thrs+="  --num-of-proc  \tthe number of input-spectra to be processed";
			cout<<show_thrs<<endl;
		}
		if (tmp == "sect"){
			string show_sect="Usage: SpectraMatch sect [options] in.txt\n";
			show_sect+="Input options:\n";
			show_sect+="  --in-file-list \tthe set of spectra with noise, input the spectra path-name list\n\n";
			show_sect+="Output options:\n";
			show_sect+="  --out-file-list\tthe set of spectra filtered, mean no noise spectra, output the spectra path-name list\n\n";
			show_sect+="ParameterSet:\n";
			show_sect+="  --num-of-proc  \tthe number of input-spectra to be processed";
			cout<<show_sect<<endl;
		}
		if (tmp == "chrg"){
			string show_chrg="Usage: SpectraMatch chrg [options] in.txt\n";
			show_chrg+="Input options:\n";
			show_chrg+="  --in-file-list \tthe set of spectra with noise, input the spectra path-name list\n\n";
			show_chrg+="Output options:\n";
			show_chrg+="  --out-file-list\tthe set of spectra filtered, mean no noise spectra, output the spectra path-name list\n\n";
			show_chrg+="ParameterSet:\n";
			show_chrg+="  --num-of-proc  \tthe number of input-spectra to be processed";
			cout<<show_chrg<<endl;
		}
		if (tmp == "five"){
			string show_five="Usage: SpectraMatch five [options] in.txt\n";
			show_five+="Input options:\n";
			show_five+="  --in-file-list \tthe set of spectra with noise, input the spectra path-name list\n\n";
			show_five+="Output options:\n";
			show_five+="  --out-file-list\tthe set of spectra filtered, mean no noise spectra, output the spectra path-name list\n\n";
			show_five+="ParameterSet:\n";
			show_five+="  --num-of-proc  \tthe number of input-spectra to be processed";
			cout<<show_five<<endl;
		}
		if (tmp == "sortlib"){
			string show_sortlib="Usage: SpectraMatch sortlib [options] in.txt\n";
			show_sortlib+="Input options:\n";
			show_sortlib+="  --lib  \tthe library path name list to search\n\n";
			show_sortlib+="Output options:\n";
			show_sortlib+="  --out  \tthe sorted path name list by mass ion\n\n";
			cout<<show_sortlib<<endl;
		}
		if (tmp == "score"){
			string show_score="Usage: SpectraMatch score [options] in.txt\n";
			show_score+="Input options:\n";
			show_score+="  --lib-for-ser\tthe set of spectra library, input the spectra library path-name list\n\n";
			show_score+="  --sam-for-sco\tthe set of sample spectra filtered, mean no noise spectra, input the spectra sample path-name list\n\n";
			show_score+="ParameterSet:\n";
			show_score+="  --num-of-proc\tthe number of input-spectra to be processed\n";
			show_score+="  --out-shape  \tthe base to shape generated thread file name ,output name ,output folder name";
			cout<<show_score<<endl;
		}
	}
	
	if (argc == 6){
		string tmp(argv[1]);
		if (tmp == "sortlib"){
			string tmp1(argv[2]);
			if (tmp1 == "--lib"){
				string tmp2(argv[4]);
				if (tmp2 == "--out"){
					char file1[200];
					char file2[200];
					string tmp3(argv[3]);
					string tmp4(argv[5]);
					strcpy(file1,tmp3.c_str());
					strcpy(file2,tmp4.c_str());
					mion_sort(file1,file2);
				}
			}
		}
	}

	if (argc == 8) {
		string tmp(argv[1]);
		if (tmp == "thrs") {
			using namespace DenoiseCutThreshold;
			string tmp1(argv[2]);
			if (tmp1 == "--in-file-list"){
				string tmp2(argv[4]);
				if (tmp2 == "--out-file-list") {
					string tmp3(argv[6]);
					if (tmp3 == "--num-of-proc") {
						int num;
						num = atoi(argv[7]);	
						string* inp;
						inp = FileList(argv[3]);
						string* oup;
						oup = FileList(argv[5]);

						for (int i = 0; i < num; i++)
						{
							string s1 = inp[i];
							string s2 = oup[i];
							
							char file1[200];
							char file2[200];
							strcpy(file1,s1.c_str());
							strcpy(file2,s2.c_str());
							DenoiseCutThreshold::DenoiseByCutThreshold(file1,file2);
						}
						
						return 0;
					}
				}
			}
		}
		if (tmp == "sect") {
			using namespace ZoneDenoise;
			string tmp1(argv[2]);
			if (tmp1 == "--in-file-list"){
				string tmp2(argv[4]);
				if (tmp2 == "--out-file-list") {
					string tmp3(argv[6]);
					if (tmp3 == "--num-of-proc") {
						int num;
						num = atoi(argv[7]);	
						string* inp;
						inp = FileList(argv[3]);
						string* oup;
						oup = FileList(argv[5]);

						for (int i = 0; i < num; i++)
						{
							string s1 = inp[i];
							string s2 = oup[i];
							
							char file1[200];
							char file2[200];
							strcpy(file1,s1.c_str());
							strcpy(file2,s2.c_str());
							ZoneDenoise::DenoiseByZone(file1,file2);
						}				
						return 0;
					}
				}
			}
		}
		if (tmp == "chrg") {
			using namespace ZoneChargeDenoise;
			string tmp1(argv[2]);
			if (tmp1 == "--in-file-list"){
				string tmp2(argv[4]);
				if (tmp2 == "--out-file-list") {
					string tmp3(argv[6]);
					if (tmp3 == "--num-of-proc") {
						int num;
						num = atoi(argv[7]);	
						string* inp;
						inp = FileList(argv[3]);
						string* oup;
						oup = FileList(argv[5]);

						for (int i = 0; i < num; i++)
						{
							string s1 = inp[i];
							string s2 = oup[i];						
							char file1[200];
							char file2[200];
							strcpy(file1,s1.c_str());
							strcpy(file2,s2.c_str());
							ZoneChargeDenoise::DenoiseByCharge(file1,file2);
						}
						
						return 0;
					}
				}
			}
		}
		if (tmp == "five") {
			using namespace ZoneChargeDenoise;
			string tmp1(argv[2]);
			if (tmp1 == "--in-file-list"){
				string tmp2(argv[4]);
				if (tmp2 == "--out-file-list") {
					string tmp3(argv[6]);
					if (tmp3 == "--num-of-proc") {
						int num;
						num = atoi(argv[7]);	
						string* inp;
						inp = FileList(argv[3]);
						string* oup;
						oup = FileList(argv[5]);

						for (int i = 0; i < num; i++)
						{
							string s1 = inp[i];
							string s2 = oup[i];						
							char file1[200];
							char file2[200];
							strcpy(file1,s1.c_str());
							strcpy(file2,s2.c_str());
							FiveDotFilter::DenoiseByFiveDot(file1,file2);
						}
						
						return 0;
					}
				}
			}
		}
	}

	if (argc==10){
		string tmp(argv[1]);
		if (tmp == "score") {
			string tmp1(argv[2]);
			if (tmp1 == "--lib-for-ser"){
				//mion_sort(argv[3]);
				sorted_lib_pathname2mw(argv[3]);
				string tmp2(argv[4]);
				string tmp3(argv[6]);
				string tmp4(argv[8]);
				string tmp5(argv[9]);
				if (tmp2 == "--sam-for-sco" && tmp3 == "--num-of-proc" && tmp4 == "--out-shape") {
					string folder=tmp5+"_Out";
	                system(("mkdir "+folder).c_str());
	                string threadname1=folder+"/"+tmp5+"1_output.txt";
	                string threadname2=folder+"/"+tmp5+"2_output.txt";
	                string threadname3=folder+"/"+tmp5+"3_output.txt";
	                string threadname4=folder+"/"+tmp5+"4_output.txt";
	                string threadname5=folder+"/"+tmp5+"5_output.txt";
	                string threadname6=folder+"/"+tmp5+"6_output.txt";
	                string threadname7=folder+"/"+tmp5+"7_output.txt";
	                string threadname8=folder+"/"+tmp5+"8_output.txt";
	                string threadname9=folder+"/"+tmp5+"9_output.txt";
	                string threadname10=folder+"/"+tmp5+"10_output.txt";
	                string threadname11=folder+"/"+tmp5+"11_output.txt";
	                string threadname12=folder+"/"+tmp5+"12_output.txt";
	                string threadname13=folder+"/"+tmp5+"13_output.txt";
	                string threadname14=folder+"/"+tmp5+"14_output.txt";
	                string threadname15=folder+"/"+tmp5+"15_output.txt";
	                string threadname16=folder+"/"+tmp5+"16_output.txt";
	                string threadname17=folder+"/"+tmp5+"17_output.txt";
	                string threadname18=folder+"/"+tmp5+"18_output.txt";
	                string threadname19=folder+"/"+tmp5+"19_output.txt";
	                string threadname20=folder+"/"+tmp5+"20_output.txt";
	                strcpy(Thread1_name,threadname1.c_str());
	                strcpy(Thread2_name,threadname2.c_str());
	                strcpy(Thread3_name,threadname3.c_str());
	                strcpy(Thread4_name,threadname4.c_str());
	                strcpy(Thread5_name,threadname5.c_str());
	                strcpy(Thread6_name,threadname6.c_str());
	                strcpy(Thread7_name,threadname7.c_str());
	                strcpy(Thread8_name,threadname8.c_str());
	                strcpy(Thread9_name,threadname9.c_str());
	                strcpy(Thread10_name,threadname10.c_str());
	                strcpy(Thread11_name,threadname11.c_str());
	                strcpy(Thread12_name,threadname12.c_str());
	                strcpy(Thread13_name,threadname13.c_str());
	                strcpy(Thread14_name,threadname14.c_str());
	                strcpy(Thread15_name,threadname15.c_str());
	                strcpy(Thread16_name,threadname16.c_str());
	                strcpy(Thread17_name,threadname17.c_str());
	                strcpy(Thread18_name,threadname18.c_str());
	                strcpy(Thread19_name,threadname19.c_str());
	                strcpy(Thread20_name,threadname20.c_str());

					N = atoi(argv[7]);
					strcpy(Name,argv[5]);
					pthread_t id1;
					int ret1;
					ret1=pthread_create(&id1,NULL,thread1,NULL);
					pthread_t id2;
					int ret2;
					ret2=pthread_create(&id2,NULL,thread2,NULL);
					pthread_t id3;
	                                int ret3;
	                                ret3=pthread_create(&id3,NULL,thread3,NULL);
	                                pthread_t id4;
	                                int ret4;
	                                ret4=pthread_create(&id4,NULL,thread4,NULL);
					pthread_t id5;
	                                int ret5;
	                                ret5=pthread_create(&id5,NULL,thread5,NULL);
	                                pthread_t id6;
	                                int ret6;
	                                ret6=pthread_create(&id6,NULL,thread6,NULL);
					pthread_t id7;
	                                int ret7;
	                                ret7=pthread_create(&id7,NULL,thread7,NULL);
	                                pthread_t id8;
	                                int ret8;
	                                ret8=pthread_create(&id8,NULL,thread8,NULL);
	                pthread_t id9;
	                int ret9;
	                ret9=pthread_create(&id9,NULL,thread9,NULL);
	                pthread_t id10;
	                int ret10;
	                ret10=pthread_create(&id10,NULL,thread10,NULL);
	                pthread_t id11;
	                int ret11;
	                ret11=pthread_create(&id11,NULL,thread11,NULL);
	                pthread_t id12;
	                int ret12;
	                ret12=pthread_create(&id12,NULL,thread12,NULL);
	                pthread_t id13;
	                int ret13;
	                ret13=pthread_create(&id13,NULL,thread13,NULL);
	                pthread_t id14;
	                int ret14;
	                ret14=pthread_create(&id14,NULL,thread14,NULL);
	                pthread_t id15;
	                int ret15;
	                ret15=pthread_create(&id15,NULL,thread15,NULL);
	                pthread_t id16;
	                int ret16;
	                ret16=pthread_create(&id16,NULL,thread16,NULL);
	                pthread_t id17;
	                int ret17;
	                ret17=pthread_create(&id17,NULL,thread17,NULL);
	                pthread_t id18;
	                int ret18;
	                ret18=pthread_create(&id18,NULL,thread18,NULL);
	                pthread_t id19;
	                int ret19;
	                ret19=pthread_create(&id19,NULL,thread19,NULL);
	                pthread_t id20;
	                int ret20;
	                ret20=pthread_create(&id20,NULL,thread20,NULL);
					pthread_join(id1,NULL);
					pthread_join(id2,NULL);
					pthread_join(id3,NULL);
	                                pthread_join(id4,NULL);
					pthread_join(id5,NULL);
	                                pthread_join(id6,NULL);
					pthread_join(id7,NULL);
	                                pthread_join(id8,NULL);
	                pthread_join(id9,NULL);
	                pthread_join(id10,NULL);
	                pthread_join(id11,NULL);
	                pthread_join(id12,NULL);
	                pthread_join(id13,NULL);
	                pthread_join(id14,NULL);
	                pthread_join(id15,NULL);
	                pthread_join(id16,NULL);
	                pthread_join(id17,NULL);
	                pthread_join(id18,NULL);
	                pthread_join(id19,NULL);
	                pthread_join(id20,NULL);
	                
	                string Out_thp=folder+"/"+tmp5+"_merge_output.txt";
	                string ping="cat "+threadname1+" "+threadname2+" "+threadname3
	                +" "+threadname4+" "+threadname5+" "+threadname6+" "+threadname7
	                +" "+threadname8+" "+threadname9+" "+threadname10+" "+threadname11
	                +" "+threadname12+" "+threadname13+" "+threadname14+" "+threadname15
	                +" "+threadname16+" "+threadname17+" "+threadname18+" "+threadname19
	                +" "+threadname20+" >> "+Out_thp;
					system(ping.c_str());
					//score(argv[4]);
					return 0;
				}
			}
		}
	}

	//if (cin_length <= 0) {  //若参数2个，cin长度小于等于0，输出usage  cerr通常用于输出错误信息与其他不属于正常逻辑的输出内容
		//cerr << usage << endl;
		//return 0;
	//}
	if (argc == 1) {  //若参数2个，cin长度小于等于0，输出usage  cerr通常用于输出错误信息与其他不属于正常逻辑的输出内容
		cerr << usage << endl;
		return 0;
	}
	return 0;
}

void *thread1(void *args)
{
        score(Name,0*N,N/20,Thread1_name);
}

void *thread2(void *args)
{
        score(Name,1*N/20,2*N/20,Thread2_name);
}

void *thread3(void *args)
{
        score(Name,2*N/20,3*N/20,Thread3_name);
}

void *thread4(void *args)
{
        score(Name,3*N/20,4*N/20,Thread4_name);
}

void *thread5(void *args)
{
        score(Name,4*N/20,5*N/20,Thread5_name);
}

void *thread6(void *args)
{
        score(Name,5*N/20,6*N/20,Thread6_name);
}

void *thread7(void *args)
{
        score(Name,6*N/20,7*N/20,Thread7_name);
}

void *thread8(void *args)
{
        score(Name,7*N/20,8*N/20,Thread8_name);
}

void *thread9(void *args)
{
        score(Name,8*N/20,9*N/20,Thread9_name);
}

void *thread10(void *args)
{
        score(Name,9*N/20,10*N/20,Thread10_name);
}

void *thread11(void *args)
{
        score(Name,10*N/20,11*N/20,Thread11_name);
}

void *thread12(void *args)
{
        score(Name,11*N/20,12*N/20,Thread12_name);
}

void *thread13(void *args)
{
        score(Name,12*N/20,13*N/20,Thread13_name);
}

void *thread14(void *args)
{
        score(Name,13*N/20,14*N/20,Thread14_name);
}

void *thread15(void *args)
{
        score(Name,14*N/20,15*N/20,Thread15_name);
}

void *thread16(void *args)
{
        score(Name,15*N/20,16*N/20,Thread16_name);
}

void *thread17(void *args)
{
        score(Name,16*N/20,17*N/20,Thread17_name);
}

void *thread18(void *args)
{
        score(Name,17*N/20,18*N/20,Thread18_name);
}

void *thread19(void *args)
{
        score(Name,18*N/20,19*N/20,Thread19_name);
}

void *thread20(void *args)
{
        score(Name,19*N/20,20*N/20,Thread20_name);
}