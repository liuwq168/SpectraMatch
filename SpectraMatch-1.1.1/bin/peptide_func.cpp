#include "peptide_score.h"
int cmp(p a, p b)
{
	return a.mu<b.mu;
}
double a,b,d;
string c;
p ans[maxn];
p ans_sort[maxn];
const int MAX_POINT = 1e4 + 10;
const int RANGE = 1e2;
const double DeleteRange = 1.01;
double Maximum=0.0;
typedef pair<double, double> pdd;
vector<pdd> point, answer;
vector<string> inString;
bool tabu[MAX_POINT];
double chg;

int cmp_sxp(sxp a,sxp b)
{
	return a.score>b.score;
}
sxp list_score[MAX];

double select(char line[])
{
	const char* split = "\t";
	char* p;
	p = strtok(line, split);
	int count = 0;
	double st;
	while (p != NULL) {
		if (count == 12) { 
			st = atof(p);
		}
		p = strtok(NULL, "\t");
		count++;
	}
	return st;
}

void mion_sort(char* filename1,char* filename2) {
	ofstream output;
	output.open(filename2);
	ifstream input;
	ifstream input1;
	input.open(filename1);
	int cnt = 0;
	while (input.getline(ans_sort[cnt].s, 300)) {		
		input1.open(ans_sort[cnt].s);
		while (input1 >> a && input1 >> b)ans_sort[cnt].mu = a;
		input1.clear();
		input1.close();	
		cnt++;
	}
	sort(ans_sort, ans_sort + cnt, cmp);
	for (int i=0;i<cnt;i++){
		output<<ans_sort[i].s<<"\t"<<ans_sort[i].mu<<"\n";
	}
	input.clear();
	input.close();
	output.close();
}

void sorted_lib_pathname2mw(char* filename)
{
	ifstream input;
	input.open(filename);
	int cnt = 0;
	while (input >> c && input >> d) {	
		strcpy(ans[cnt].s,c.c_str());
		ans[cnt].mu = d;
		cnt++;
	}
	input.clear();
	input.close();
}

double vectorsum(vector<double>::iterator start, vector<double>::iterator start1, vector<double>::iterator end1, vector<double>::iterator endl1)
{
	double sum = 0;
	while (start != end1 && start1 != endl1)
	{
		sum += (*start) * (*start1);
		start++;
		start1++;
	}
	return sum;
}

double vectorSum(vector<double>::iterator start, vector<double>::iterator end1)
{
	double sum = 0;
	while (start != end1)
	{
		sum += (*start);
		start++;
	}
	return sum;
}

double vector_sqrt_Sum(vector<double>::iterator start, vector<double>::iterator end1)
{
	double sum = 0;
	while (start != end1)
	{
		sum += (*start)*(*start);
		start++;
	}
	return sqrt(sum);
}

bool compare(pair<double, double> p1, pair<double, double> p2) {
	return p1.second > p2.second;
}

bool compare1(pair<double, double> p1, pair<double, double> p2) {
	return p1.first < p2.first;
}

string peptide_name(string s) {
	int start, end;
	start = s.find_last_of("//");
	s = s.replace(0, start + 1, "");
	end = s.find_first_of("_");
	s = s.substr(0, end);
	return s;
}

void Print(char* filename3, char* filename4, char* filename5)
{	
	ofstream out222;
	out222.open(filename5, ios::app);
	out222.precision(12);
	out222.setf(ios_base::fixed, ios_base::floatfield);
	
	ifstream fin71;
	fin71.open(filename3, ifstream::in);
	fin71.precision(10);
	fin71.setf(ios_base::fixed, ios_base::floatfield);
	vector<pair<double, double> > petide3;
	pair<double, double> pte3;
	while (fin71)
	{
		double s1; double s2;
		fin71 >> s1; fin71 >> s2;
		
		pte3.first = s1; pte3.second = s2;
		petide3.push_back(pte3);
	}
	fin71.clear();
	fin71.close();
	
	ifstream fin72;
	fin72.open(filename4, ifstream::in);
	fin72.precision(10);
	fin72.setf(ios_base::fixed, ios_base::floatfield);
	vector<pair<double, double> > petide4;
	pair<double, double> pte4;
	while (fin72)
	{
		double s11; double s12;
		fin72 >> s11; fin72 >> s12;
		pte4.first = s11; pte4.second = s12;
		petide4.push_back(pte4);
	}
	fin72.clear();
	fin72.close();
	
	int N2, P2, J2, G1, G2;
	N2 = petide3.size() - 2;
	P2 = petide4.size() - 2;
	J2 = min(N2, P2);
	G1 = N2 + 2;
	G2 = P2 + 2;
	double Max1, Max2;
	vector<pair<double, double> > petide5;
	pair<double, double> pte5;
	vector<pair<double, double> > petide6;
	pair<double, double> pte6;
	vector<double> sam;
	vector<double> lib;
	int count = 0;
	double minmin;
	
	if (fabs(petide3[petide3.size() - 1].first - petide4[petide4.size() - 1].first) / (petide4[petide4.size() - 1].first) < 1e-5)
	{
		petide3.pop_back();
		petide4.pop_back();
		Max1 = petide3[petide3.size() - 1].first;
		Max2 = petide4[petide4.size() - 1].first;
		petide3.pop_back();
		petide4.pop_back();

	

		for (int i = 0; i < petide3.size(); i++)
		{
			if (petide3[i].first - petide4[0].first < -0.010) continue;
			double mark = 0x3f3f3f3f;
			double epss = 1e-8;
			int flag=0;
			for (int j = 1; j < petide4.size(); j++)
			{
				if (petide3[i].first - petide4[j].first >= -0.010 && petide3[i].first - petide4[j].first <= 0.010)
				{
					double temp = fabs(petide3[i].first - petide4[j].first);
					if (temp+epss <= mark) {
						temp = mark;
						pte5.first = petide3[i].first;
						pte5.second = petide3[i].second;
						petide5.push_back(pte5);
						pte6.first = petide4[j].first;
						pte6.second = petide4[j].second;
						petide6.push_back(pte6);
						sam.push_back(petide3[i].second);
						lib.push_back(petide4[j].second);
						if(!flag)
							count++;
						flag = 1;
					}
				}
				else if (petide3[i].first - petide4[j].first > 0.010)continue;
				else break;
				}
			}
		}
		
		double vector_score1 = count;
		double vector_score2 = N2;

		int NN2 = petide5.size();
		int PP2 = petide6.size();
	
		if (NN2 > 1)
		{
			if (petide3.size() > 1 && petide4.size() > 1)
			{

				if (N2 > J2)
				{
					sort(petide3.begin(), petide3.end(), compare);

					sort(petide3.begin(), petide3.end(), compare1);
				}
				else
				{
					sort(petide4.begin(), petide4.end(), compare);

					sort(petide4.begin(), petide4.end(), compare1);
				}

				int vtc = petide5.size();
				
				vector<double> vcf3_1x;
				vector<double> vcf3_1y;
				vector<double> vcf4_1x;
				vector<double> vcf4_1y;
				for (int i = 0; i < vtc; i++)
				{
					vcf3_1x.push_back(petide5[i].first);
					vcf3_1y.push_back(petide5[i].second);
					vcf4_1x.push_back(petide6[i].first);
					vcf4_1y.push_back(petide6[i].second);
				}
				
				double p2osx1 = 0.0;
				double p2osx3 = 0.0;
				double p2osx4 = 0.0;
				p2osx1 = vectorsum(vcf3_1x.begin(), vcf4_1x.begin(), vcf3_1x.end(), vcf4_1x.end());
				p2osx3 = vectorsum(vcf3_1x.begin(), vcf3_1x.begin(), vcf3_1x.end(), vcf3_1x.end());
				p2osx4 = vectorsum(vcf4_1x.begin(), vcf4_1x.begin(), vcf4_1x.end(), vcf4_1x.end());
				double dot2_productx = 0.0;
				dot2_productx = pow(p2osx1, 2) / (p2osx3 * p2osx4);
				
				double p2osyy1 = 0.0;
				double p2osyy3 = 0.0;
				double p2osyy4 = 0.0;
				double ST1 = 0.0;
				double ST2 = 0.0;
				ST1 = vectorSum(vcf3_1y.begin(), vcf3_1y.end());
				ST2 = vectorSum(vcf4_1y.begin(), vcf4_1y.end());

				for (int j = 0; j < vcf3_1y.size(); j++)
				{
					p2osyy1 += vcf3_1y[j] * 1 / (1 + vcf3_1y[j] / (ST1 - 0.5)) * vcf4_1y[j] * 1 / (1 + vcf4_1y[j] / (ST2 - 0.5));
					p2osyy3 += pow(vcf3_1y[j] * 1 / (1 + vcf3_1y[j] / (ST1 - 0.5)), 2);
					p2osyy4 += pow(vcf4_1y[j] * 1 / (1 + vcf4_1y[j] / (ST2 - 0.5)), 2);
				}
				double dot2_productyy = 0.0;
				dot2_productyy = pow(p2osyy1, 2) / (p2osyy3 * p2osyy4);
    				
				
				
				
				double peak_intensity1=0.0;
				double peak_intensity2=0.0;
				double peak_intensity3=0.0;
				double peak_score = 0.0;
				peak_intensity1= vectorsum(sam.begin(),lib.begin(),sam.end(),lib.end());
				peak_intensity2= vector_sqrt_Sum(sam.begin(),sam.end());
				peak_intensity3= vector_sqrt_Sum(lib.begin(),lib.end());
				peak_score = peak_intensity1 / (peak_intensity2 * peak_intensity3);

				out222 << filename3 << '\t' << filename4 << '\t';
				out222 << dot2_productx << '\t' << dot2_productyy<<'\t';
				out222 << (vector_score1) / sqrt(vector_score1*vector_score2) << '\t';
				out222 << (vector_score1) / (vector_score2) << '\t';
				out222 << peak_score<<"\t";
				out222 << peptide_name(filename4) << "\n"; 

                                
				out222.close();
				
				vector<double>().swap(vcf3_1x);
				vector<double>().swap(vcf3_1y);
				vector<double>().swap(vcf4_1x);
				vector<double>().swap(vcf4_1y);
				vector<double>().swap(sam);
				vector<double>().swap(lib);
				
				
			}
 		
		}
	vector<pair<double, double> >().swap(petide3);
	vector<pair<double, double> >().swap(petide4);
	vector<pair<double, double> >().swap(petide5);
	vector<pair<double, double> >().swap(petide6);
	
	
}

int BinarySearchRecursive_left(int low, int high, double key) 
{
	
	if (low > high)
		return -1;
	int mid = (low + high) / 2;
	
	if (ans[mid + 1].mu >= key && ans[mid].mu < key)
		return mid;
	else if (ans[mid + 1].mu < key)
		return BinarySearchRecursive_left(mid + 1, high, key);
	else if (ans[mid].mu >= key)
		return BinarySearchRecursive_left(low, mid - 1, key);
}

int BinarySearchRecursive_right(int low, int high, double key) 
{
	if (low > high)
		return -1;
	int mid = (low + high) / 2;
	if (ans[mid - 1].mu <= key && ans[mid].mu > key)
		return mid;
	else if (ans[mid].mu <= key)
		return BinarySearchRecursive_right(mid + 1, high, key);
	else if (ans[mid - 1].mu > key)
		return BinarySearchRecursive_right(low, mid - 1, key);
}
		
string *FileList(char* filename)
{
	ifstream inFile;
	string tmpStr("");
	string* tp;
	tp = new string[200000]; //50951
	long index = 0;
	inFile.open(filename, ios::in);
	while (getline(inFile, tmpStr))
	{
		tp[index] = tmpStr;
		index += 1;
	}
	inFile.clear();
	inFile.close();
	return tp;
	delete[] tp;
}
