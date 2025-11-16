#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <algorithm>
using namespace std;

ifstream INPUT;
ofstream OUTPUT;
ifstream SIZE;

struct leap // Структура скачка
{
	double t; // время;
	char label; // label;
};

vector <leap> leaps; // скачки;

bool comp(leap A, leap B) { // Компаратор для сортировки;
	return A.t < B.t;
}

int main(){
    int MAXHIGH = 0;
    SIZE.open("size_simulate_evol.txt");
    string number;
    while(SIZE >> number){
        int n;
        SIZE >> n;

        leaps.clear();
        string file_name;
        file_name = "Result_Ripser_simulate_evol/" + number + ".persistence";
        INPUT.open(file_name.c_str());
        file_name = "Result_b0curve_simulate_evol/" + number + ".b0curve";
        OUTPUT.open(file_name.c_str());
                
        string s;
        INPUT>>s;
        while(s!="0:"){
            INPUT>>s;
        }

        char ch;
        double b,d;
        while(INPUT>>ch>>b>>ch>>d>>ch){
            leap q;
            q.t = b;
            q.label = 'b';
            leaps.push_back(q);
            q.t = d;
            q.label = 'd';
            leaps.push_back(q);
        }
        sort(leaps.begin(), leaps.end(), comp);
        int cnt=1;
        OUTPUT<<0.0<<' '<<cnt<<endl;
        for(int i=0;i<(int)leaps.size();i++){
            if(leaps[i].label == 'b'){
                cnt++;
            }
            if(leaps[i].label == 'd'){
                cnt--;
            }
            OUTPUT<<leaps[i].t<<' '<<cnt<<endl;
            MAXHIGH = max(MAXHIGH, cnt);
        }
        OUTPUT<<1.0<<' '<<cnt<<endl;

        leaps.clear();

        INPUT.close();
        OUTPUT.close();
    }
    cout<<MAXHIGH<<endl;
    return 0;
}
