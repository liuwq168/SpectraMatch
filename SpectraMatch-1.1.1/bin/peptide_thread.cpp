#include<bits/stdc++.h>
#include "peptide_score.h"
using namespace std;

int score(char* Filename,int n,int m,char* file_name)
{
	char name[300];
	strcpy(name,file_name);
	
	string* sx=NULL;
        //char filename[200];
        //strcpy(filename,"./mgf_csv_rep44.txt");
        //strcpy(filename,Filename);
	sx = FileList(Filename);
	
	for (int j = n; j < m; j++)
	{
		
		string s=sx[j];
		char file[300];
  		strcpy(file,s.c_str());
		ifstream fin77;
		fin77.open(file, ifstream::in);
		fin77.precision(10);
		fin77.setf(ios_base::fixed, ios_base::floatfield);
		vector<double> vcf23_x;
		vector<double> vcf23_y;
		
		while (fin77)
		{
			double s21; double s22;
			fin77 >> s21; fin77 >> s22;
			//cout<<s21<<","<<s22<<endl;
			
			vcf23_x.push_back(s21); vcf23_y.push_back(s22);
		}
		fin77.clear();
		fin77.close();
		
		double ll = 0.0;
		ll = vcf23_x[vcf23_x.size() - 1];
		double eps = 0.5;
		double st, ed;
		st = ll - 0.5;
		ed = ll + 0.5;

		int ll1 = 0,ll2;
		
		ll1 = BinarySearchRecursive_left(0, maxn, st);
		ll2 = BinarySearchRecursive_right(0, maxn, ed);
		
  		
		
		if (ll1 != -1 && ll2 != -1)
		{
			for (int k = ll1; k < ll2; k++)
			{
				
				Print(file, ans[k].s, name);
			}
		}
		else if (ll1 != -1 && ll2 == -1)
		{
			for (int k = ll1; k < maxn; k++)
			{
				
				Print(file, ans[k].s, name);
			}
		}
		else if (ll1 == -1 && ll2 != -1)
		{
			for (int k = 0; k < ll2; k++)
			{
				
				Print(file, ans[k].s, name);
			}
		}
		else
		{
			continue;
		}

	}
	
 	return 0;
}