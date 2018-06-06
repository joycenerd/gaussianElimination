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

void toOne(int row){
    float divisor;
    int j;
    swapTimes=0;
    if(abs(augmentMatrix[row][row])>0.0001) flag=1;
    if(augmentMatrix[row][row]>0.0001) divisor=augmentMatrix[row][row];
    else if(augmentMatrix[row][row]<-0.0001) divisor=augmentMatrix[row][row];
    else if(abs(augmentMatrix[row][row])<=0.0001){
        swapRow(row);
        swapTimes++;
        if(flag) toOne(row);
    }
    //printf("divisor=%f\n",divisor);
    if(flag){
        determinant*=divisor;
        for(j=0;j<n+1;j++) augmentMatrix[row][j]/=divisor;
    }
}

void gaussianElimination(){
    /*outputFile=fopen("/Users/joycechin/Desktop/EX7/gaussianElimination/gaussianElimination/output.out","w");
    assert(outputFile!=NULL);*/
    fprintf(outputFile,"Case %d:\n",cases);
    int row=0;
    int i,j;
    float divisor;
    determinant=1;
    while(row<n){
        for(j=0;j<n;j++){
            if(abs(augmentMatrix[row][j])>0.0001){
                flag=1;
                break;
            }
        }
        if(j==n) flag=0;
        toOne(row);
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
    FILE *inputFile;
    inputFile=fopen("/Users/joycechin/Desktop/EX7/gaussianElimination/ex7_sample.in","r");
    assert(inputFile!=NULL);
    outputFile=fopen("/Users/joycechin/Desktop/EX7/gaussianElimination/gaussianElimination/output.out","w");
    assert(outputFile!=NULL);
    float coefficient[N][N];
    float b[N];
    int i,j;
    cases=0;
    while(1){
        fscanf(inputFile,"%d",&n);
        if(n==0) break;
        for(i=0;i<n;i++) for(j=0;j<n;j++) fscanf(inputFile,"%f",&coefficient[i][j]);
        for(i=0;i<n;i++) fscanf(inputFile,"%f",&b[i]);
        for(i=0;i<n;i++) for(j=0;j<n;j++) augmentMatrix[i][j]=coefficient[i][j];
        for(i=0;i<n;i++) augmentMatrix[i][n]=b[i];
        gaussianElimination();
        cases++;
    }
    fclose(inputFile);
    fclose(outputFile);
    /*while(1){
        scanf("%d",&n);
        if(n==0) break;
        for(i=0;i<n;i++) for(j=0;j<n;j++) scanf("%f",&coefficient[i][j]);
        for(i=0;i<n;i++) scanf("%f",&b[i]);
        for(i=0;i<n;i++) for(j=0;j<n;j++) augmentMatrix[i][j]=coefficient[i][j];
        for(i=0;i<n;i++) augmentMatrix[i][n]=b[i];
#if DEBUG==1
        for(i=0;i<n;i++){
            for(j=0;j<n+1;j++){
                if(j==n) printf("%f\n",augmentMatrix[i][j]);
                else printf("%f ",augmentMatrix[i][j]);
            }
        }
#endif
        gaussianElimination();
    }*/
    return 0;
}
