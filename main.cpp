#include <iostream>
#include <string>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include <sstream>
using namespace std;


static const char base[4]={'A', 'T', 'G', 'C'};

char randomBase(){
	int randNum;
	randNum=rand() % 4;
	return base[randNum];
}

bool checkRandNum (int numBef, int numAft){
	bool test;
        if (numAft==numBef){
            test=true;
        }
        else test=false;
        return test;
}

bool checkRandLine(int * randSubsLine, int itr){
	bool test;
	for(int i=itr-1; i>=0; i--){
		if (randSubsLine[itr]==randSubsLine[i]){
			test=true;
			break;
		}
		test=false;
	}
return test; 
}

int main(){
	int i=0, j=0, k=0, m=0, subsNum=5, subsLineNum=30, simReadNum=100, randBaseLoc[subsNum], randSubsLine[subsLineNum], randRange=22, randStart=2, randEnd, lineNum=6000;
    string in="mature.hsa.dna.fa", out="output.fastq", qtyLine="";
	ifstream infile;
    ofstream outfile;
	string line, totalLine[lineNum], totalLineStr, headLine[lineNum], subsLine, tempLine, fastqStr, randStartStr, randEndStr, randBaseLocChar, randBaseLocStr, subsLineRng, totalLineRng;
    char read[30], tempChar;
	infile.open(in.c_str());
	while(getline(infile, line)){	 
		if(line[i]!='>'){
            totalLine[j]=line;
            ++j;
	}	
    else {
        headLine[k]=line;
        k++;
    }
    }
	infile.close();
    
for (i=0; i<lineNum; i++){
totalLineStr=totalLine[i];
//--------------------------------------------------------------------------------------------------------set substitution location
        subsLine="";
        if (totalLine[i]=="")break;
    for (j=0; j<subsNum; j++){
        randBaseLoc[j]=rand()% (totalLine[i].length()-1);
        if (j>0){
            while(checkRandLine(randBaseLoc, j)){
                randBaseLoc[j]=rand()% (totalLine[i].length()-1);
            }
        }
        // cout<<randBaseLoc[j]<<"\t";
    }
    int* end = randBaseLoc + subsNum;
    sort(randBaseLoc, end);

    // cout<<endl;

//-----------------------------------------------------------------------------------------------------applying substitution
    tempLine=totalLine[i];
    k=0;
    for(j=0; j<totalLine[i].length(); j++){
        if(randBaseLoc[k]==j){
            tempChar=randomBase();
            while(tempChar==tempLine[j]){
                tempChar=randomBase();
            }
            tempLine[j]=tempChar;
            ++k;
        }
        subsLine+=tempLine[j];
}
// cout<<totalLine[i]<<endl<<subsLine<<endl;

//-------------------------------------------------------------------------------------------------set number of line that contain substituted base
    for(j=0; j<subsLineNum; j++){
        randSubsLine[j]=rand() % simReadNum;
        if (j>0){
            while(checkRandLine(randSubsLine, j)){
                randSubsLine[j]=rand() % simReadNum;
            }
        }
    }

    end = randSubsLine + subsLineNum;
	sort(randSubsLine,end);
    
    // for(j=0; j<subsLineNum; j++){
    //     cout<<randSubsLine[j]<<"   ";
    // }
//-----------------------------------------------------------------------------------------------------storing fastq format sequence into variable
int l=0;
    for(j=0; j<simReadNum; j++){
        switch (totalLineStr.length()){
            case 17:
                randRange=16;
                randStart=0;
                break;
            case 18:
                randRange=rand() % 2+16;
                randStart=rand() % 2;
                break;
            case 19:
                randRange=rand() % 3+16;
                randStart=rand() % 3;
                break;
            case 20:
                randRange=rand() % 4+16;
                randStart=rand() % 4;
                break;
            case 21:
                randRange=rand() % 5+16;
                randStart=rand() % 5;
                break;
            case 22:
                randRange=rand() % 6+16;
                randStart=rand() % 6;
                break;
        }


        randEnd=randStart+randRange;
        if(randEnd>=totalLineStr.length()-1){
            randEnd=totalLineStr.length()-1;
        }
        randStartStr=to_string(randStart+1);
        randEndStr=to_string(randEnd+1);

        randBaseLocStr="";
        for(k=0; k<subsNum; k++){
            m=randBaseLoc[k]+1;
            randBaseLocChar=to_string(m);
            randBaseLocStr+=randBaseLocChar+" "+totalLineStr[m-1]+">"+subsLine[m-1]+", ";
        }

        qtyLine="";subsLineRng="";totalLineRng="";

        if (j==randSubsLine[l]){
        fastqStr+=headLine[i]+"; "+randStartStr+ ":" + randEndStr+"; "+ randBaseLocStr+"\n";
            for(k=randStart; k<=randEnd; k++){
                subsLineRng+=subsLine[k];
            }
        for(k=0; k<subsLineRng.length(); k++){
            qtyLine+="f";
        }
        fastqStr+=subsLineRng+"\n+\n"+qtyLine+"\n";
        ++l;
        }
        
        else{
        fastqStr+=headLine[i]+"; "+randStartStr+ ":" + randEndStr+"; "+"\n";
            for(k=randStart; k<=randEnd; k++){
                totalLineRng+=totalLineStr[k];
            }
        for(k=0; k<totalLineRng.length(); k++){
            qtyLine+="f";
        }
        fastqStr+=totalLineRng+"\n+\n"+qtyLine+"\n";            
        }
    }
}
outfile.open(out.c_str());
    outfile<<fastqStr;
outfile.close();
    return 0;
}