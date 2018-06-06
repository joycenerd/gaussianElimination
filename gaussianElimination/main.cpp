//
//  main.cpp
//  gaussianElimination
//
//  Created by Joyce Chin on 2018/6/5.
//  Copyright © 2018年 Joyce Chin. All rights reserved.
//

#include <cstdio>
#include <cmath>
#include <cassert>
#define N 120
#define DEBUG 0
using namespace std;

FILE *outputFile;
float augmentMatrix[N][N],determinant;
int n,cases;
int flag,swapTimes;

//錊後的驗算副程式-->把項諒解帶回去Ax看會不會等於b
int checkValid(float x[]){
    float ans=0;
    int check=1;
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            ans+=augmentMatrix[i][j]*x[j];
        }
        if(abs(ans-augmentMatrix[i][n])>0.0001){
            check=0;
            break;
        }
    }
    if(check) return 1;
    else return -1;
}

//尋找和行行互換
void swapRow(int row){
    int i,j;
    float temp[N];
    flag=0;
    for(i=row;i<n;i++){
        if(augmentMatrix[i][row]!=0){
            for(j=0;j<n+1;j++){
                temp[j]=augmentMatrix[row][j];
                augmentMatrix[row][j]=augmentMatrix[i][j];
                augmentMatrix[i][j]=temp[j];
            }
            flag=1;
            break;
        }
    }
}

//將augmentMatrix[row][row]變成1
void toOne(int row){
    float divisor;  //要除divisor才會變成1
    int j;
    swapTimes=0;  //行行互換的次數，如果最後為奇數次行列式值要多加個負號
    if(abs(augmentMatrix[row][row])>0.0001){
        flag=1;
        divisor=augmentMatrix[row][row];
    }
    //如果那個算過小趨近於0，那就從下面剩下的行數相對位置中去尋找大於0的元素
    else if(abs(augmentMatrix[row][row])<=0.0001){
        swapRow(row);
        swapTimes++;
        if(flag) toOne(row);
    }
    //printf("divisor=%f\n",divisor);
    if(flag){
        determinant*=divisor;  //邊做高斯消去邊算行列式值（原本的對角線元素相乘）
        for(j=0;j<n+1;j++) augmentMatrix[row][j]/=divisor;
    }
}

//高斯消去法的主體
void gaussianElimination(){
    fprintf(outputFile,"Case %d:\n",cases);
    int row=0;
    int i,j;
    float divisor;
    determinant=1;
    //一行一行的做消去
    while(row<n){
        //判斷是不是那一行的係數全為0
        for(j=0;j<n;j++){
            if(abs(augmentMatrix[row][j])>0.0001){
                flag=1;
                break;
            }
        }
        if(j==n) flag=0;
        toOne(row);  //將augmentMatrix[row][row]這個元素變為1
        //如果行列式值不為0就接著做高斯消去（由上面的toOne和swapRow去判斷）
        //最後回消成只剩對角線為1，向量解會是擴增矩陣的最後一個column
        if(flag){
            for(i=0;i<n;i++){
                if(i!=row){
                    divisor=-augmentMatrix[i][row];
                    for(j=0;j<n+1;j++)
                        augmentMatrix[i][j]+=augmentMatrix[row][j]*divisor;
                }
            }
        }
        else{
            fprintf(outputFile,"The determinant is 0\nThis equation is linear dependence\n");
            break;
        }
#if DEBUG //== 1
        for(i=0;i<n;i++){
            for(j=0;j<n+1;j++){
                if(j==n) printf("%f\n",augmentMatrix[i][j]);
                else printf("%f ",augmentMatrix[i][j]);
            }
        }
#endif
        row++;
    }
    float x[N];
    //如果有解救把行列式值、向量解和是否為正確解輸出到outputFile
    if(row==n){
        for(i=0;i<n;i++) x[i]=augmentMatrix[i][n];
        if(swapTimes%2==1) determinant=-determinant;
        fprintf(outputFile,"The determinant is %f\n",determinant);
        fprintf(outputFile,"The solution vector is: ");
        for(i=0;i<n;i++) fprintf(outputFile,"%f ",x[i]);
        fprintf(outputFile,"\n");
        int valid=checkValid(x);
        if(valid) fprintf(outputFile,"This is a valid solution\n");
        else fprintf(outputFile,"This is an invalid solution\n");
    }
    fprintf(outputFile,"======================================================\n");
}

int main()
{
    FILE *inputFile;  //輸入檔
    //開檔讀檔
    inputFile=fopen("/Users/joycechin/Desktop/EX7/gaussianElimination/ex7_sample.in","r");
    assert(inputFile!=NULL);
    outputFile=fopen("/Users/joycechin/Desktop/EX7/gaussianElimination/gaussianElimination/output.out","w");
    assert(outputFile!=NULL);
    float coefficient[N][N];
    float b[N];
    int i,j;
    cases=0;  //記錄現在是第幾比測資
    while(1){
        fscanf(inputFile,"%d",&n);  //n是我的dimension（未知數）
        if(n==0) break;
        //Ａx=b的A矩陣，就是係數構成的矩陣
        for(i=0;i<n;i++) for(j=0;j<n;j++) fscanf(inputFile,"%f",&coefficient[i][j]);
        //Ａx=b的b
        for(i=0;i<n;i++) fscanf(inputFile,"%f",&b[i]);
        //把A和b做成擴增矩陣
        for(i=0;i<n;i++) for(j=0;j<n;j++) augmentMatrix[i][j]=coefficient[i][j];
        for(i=0;i<n;i++) augmentMatrix[i][n]=b[i];
        gaussianElimination();
        cases++;
    }
    fclose(inputFile);
    fclose(outputFile);
    return 0;
}
