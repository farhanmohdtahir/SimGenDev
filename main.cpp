#include <iostream>
#include <string>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include <sstream>
#include <getopt.h>
using namespace std;


static const char base[4]={'A', 'T', 'G', 'C'};

void help();

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

int main(int argc, char * argv []){
	int i=0, j=0, k=0, m=0, readRange=0, subsNum=5, subsLineNum=30, simReadNum=100, randBaseLoc[subsNum], randSubsLine[subsLineNum], randRange=22, randStart=2, randEnd, lineNum=60000, opt=0;
    int totalRange=6, subsPres[totalRange], rangeNum[totalRange], noSubsRange[totalRange];
    string in="mature.hsa.dna.fa", out="output.fastq", qtyLine="", headLineTemp;
	ifstream infile;
    ofstream outfile;
	string line, totalLine[lineNum], totalLineStr, headLine[lineNum], subsLine, tempLine, fastqStr, randStartStr, randEndStr, randBaseLocChar, randBaseLocStr, subsLineRng, totalLineRng;
    char read[30], tempChar;
    double subsPerc=10, subsNumTemp;

//--------------------------------------------------------------------------------------------------------getopt
    bool option=false;
    //getopt function	
    static struct option long_options[] = {
        {"help",                       no_argument,       0,  'h' },
        {"inputFile",              required_argument,     0,  'i' },
        {"outputFile",             required_argument,     0,  'o' },                         
        {0,                               0,              0,   0  }
    };
    
        int long_index =0;
    while ((opt = getopt_long_only(argc, argv,"", 
                   long_options, &long_index )) != -1) {
        switch (opt) {
             case 'h' : help();option = true;
                 break;
             case 'i' : in = optarg;option = true;
                 break;
             case 'o' : out = optarg;option = true;
                 break;                                               
             default: help(); 
                 exit(EXIT_FAILURE);
        }
    }
    
	infile.open(in.c_str());
	while(getline(infile, line)){	 
		if(line[i]!='>'){
            totalLine[j]=line;
            ++j;
	}	
    else {
        headLineTemp=line;
        headLineTemp[0]='@';
        headLine[k]=headLineTemp;
        k++;
    }
    }
	infile.close();
    
for (i=0; i<lineNum; i++){
totalLineStr=totalLine[i];
subsNumTemp=(subsPerc/100)*totalLineStr.length();
subsNum=1;//round(subsNumTemp);
for (k=0; k<totalRange; k++){
    rangeNum[k]=0;
    subsPres[k]=0;
}
//-----------------------------------------------------------------------------------------------------set substitution location
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

        readRange=randEnd-randStart+1;
        if (readRange==17){
            ++rangeNum[0];
            for(k=0; k<subsNum; k++){
                if (randBaseLoc[k]>=randStart && randBaseLoc[k]<=randEnd){
                    ++subsPres[0];
                    break;
                }
        }
        }
        else if (readRange==18){
            ++rangeNum[1];
            for(k=0; k<subsNum; k++){
                if (randBaseLoc[k]>=randStart && randBaseLoc[k]<=randEnd){
                    ++subsPres[1];
                    break;
                }            
        }
        }
        else if (readRange==19){
            ++rangeNum[2];
            for(k=0; k<subsNum; k++){
                if (randBaseLoc[k]>=randStart && randBaseLoc[k]<=randEnd){
                    ++subsPres[2];
                    break;
                }
        }
        }
        else if (readRange==20){
            ++rangeNum[3];
            for(k=0; k<subsNum; k++){
                if (randBaseLoc[k]>=randStart && randBaseLoc[k]<=randEnd){
                    ++subsPres[3];
                    break;
                }
        }
        }
        else if (readRange==21){
            ++rangeNum[4];
            for(k=0; k<subsNum; k++){
                if (randBaseLoc[k]>=randStart && randBaseLoc[k]<=randEnd){
                    ++subsPres[4];
                    break;
                }
        }
        }
        else if (readRange==22){
            ++rangeNum[5];
            for(k=0; k<subsNum; k++){
                if (randBaseLoc[k]>=randStart && randBaseLoc[k]<=randEnd){
                    ++subsPres[5];
                    break;
                }
        }
        }
        randStartStr=to_string(randStart+1);
        randEndStr=to_string(randEnd+1);

        randBaseLocStr="";
            // for(k=0; k<subsNum; k++){
            //     m=randBaseLoc[k]+1;
            //     randBaseLocChar=to_string(m);
            //     randBaseLocStr+=randBaseLocChar+" "+totalLineStr[m-1]+">"+subsLine[m-1];
            //     if (k<subsNum-1) randBaseLocStr+=",";
            // }        

        qtyLine="";subsLineRng="";totalLineRng="";

        if (j==randSubsLine[l]){
            for(k=0; k<subsNum; k++){
                if (randBaseLoc[k]>=randStart&&randBaseLoc[k]<=randEnd){
                m=randBaseLoc[k]-randStart+1;
                randBaseLocChar=to_string(m);
                randBaseLocStr+=randBaseLocChar+" "+totalLineStr[randBaseLoc[k]]+">"+subsLine[randBaseLoc[k]];
                if (k<subsNum-1) randBaseLocStr+=",";
                }
            }    
        fastqStr+=headLine[i]+" "+randStartStr+ ":" + randEndStr+"  "+ randBaseLocStr+"\n";
            for(k=randStart; k<=randEnd; k++){
               subsLineRng+=subsLine[k];
            }
        for(k=0; k<subsLineRng.length(); k++){
            qtyLine+="h";
        }
        fastqStr+=subsLineRng+"\n+\n"+qtyLine+"\n";
        ++l;
        }
        
        else{
        fastqStr+=headLine[i]+" "+randStartStr+ ":" + randEndStr+" "+"\n";
            for(k=randStart; k<=randEnd; k++){
                totalLineRng+=totalLineStr[k];
            }
        for(k=0; k<totalLineRng.length(); k++){
            qtyLine+="h";
        }
        fastqStr+=totalLineRng+"\n+\n"+qtyLine+"\n";            
        }
    }
    for (j=0; j<totalRange; j++){
        noSubsRange[j]=rangeNum[j]-subsPres[j];
    }
    cout<<headLine[i]<<endl<<totalLine[i].length()<<"+"<<simReadNum<<"+"<<subsNum<<"+"<<rangeNum[0]<<"+"<<rangeNum[1]<<"+"<<rangeNum[2]<<"+"<<rangeNum[3]<<
    "+"<<rangeNum[4]<<"+"<<rangeNum[5]<<"+"<<noSubsRange[0]<<"+"<<noSubsRange[1]<<"+"<<noSubsRange[2]<<"+"<<noSubsRange[3]<<
    "+"<<noSubsRange[4]<<"+"<<noSubsRange[5]<<"+"<<subsPres[0]<<"+"<<subsPres[1]<<"+"<<subsPres[2]<<"+"<<subsPres[3]<<
    "+"<<subsPres[4]<<"+"<<subsPres[5]<<endl<<endl;
}
outfile.open(out.c_str());
    outfile<<fastqStr;
outfile.close();
    return 0;
}

void help(){
	cout<<"FASTQSimGen"<<endl;
	cout<<"Usage: ./ar -i inputFile -o outputFile -m mutationPercentage"<<endl;
	cout<<"-i [required argument] - name of input file (.fa/.fasta)\n";
	cout<<"-o [required argument] - name of output file (.fq/fastq)\n";
	cout<<"Thank You\n";
}