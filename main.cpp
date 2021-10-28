//
//  main.cpp
//  inverse
//
//  Created by 张凯 on 2021/10/21.
//  Copyright © 2021 张凯. All rights reserved.
//

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <ctime>
#define pi 3.141592653
using namespace std;
struct complex{
    double re=0.0;
    double im=0.0;
    friend struct complex operator *(const struct complex&x,const struct complex&y){
        struct complex ans={0.0,0.0};
        ans.re = x.re*y.re-x.im*y.im;
        ans.im = x.im*y.re+x.re*y.im;
        return ans;
    }
    friend struct complex operator +(const struct complex&x,const struct complex&y){
        struct complex ans={0.0,0.0};
        ans.re = x.re+y.re;
        ans.im = x.im+y.im;
        return ans;
    }
    friend struct complex operator -(const struct complex&x,const struct complex&y){
        struct complex ans={0.0,0.0};
        ans.re = x.re-y.re;
        ans.im = x.im-y.im;
        return ans;
    }
    void print(){
        cout<<this->re<<"+"<<this->im<<"i"<<endl;
    }
};
vector<struct complex>x,y,W;
// caculating the W(N,K) using Euler's formula
struct complex WNK(int n,int k){
    struct complex ans={0.0,0.0};
    ans.re = cos(2*pi*k/n);
    ans.im = -sin(2*pi*k/n);
    ans.re = abs(ans.re)>1e-6?ans.re:0;
    ans.im = abs(ans.im)>1e-6?ans.im:0;
    return ans;
}
// caculating the n points butterfly graph nets of x[n]
vector<struct complex> butterfly(vector<struct complex> &x,int stage=1){
    int N = x.size();
    int group = int(pow(2, stage-1));
    int step = N/group;
    for (int i=1;i<=group;i++){
        int start = (i-1)*step;
        int stop = i*step-1;
        int half = (start + stop + 1)/2-1;
        for (int n=start; n<=half; n++)
            y[n] = x[n] + x[n+step/2];
        for (int n=half+1; n<=stop; n++)
            y[n] = (x[n-step/2]-x[n])*W[n%(step/2)*group];
    }
    return y;
}
// inverse x[n]
vector<struct complex> inverse(vector<struct complex> &x){
    int j=0;
    int N = x.size();
    for(int i = 0; i < N-1; i ++){
        if(i < j){
            auto temp = x[i];
            x[i] = x[j];
            x[j] = temp;
        }
        int k = N >> 1;
        while(k<=j)
            j = j-k,k>>=1;
        j = j+k;
    }
    return x;
}
// caculating the n points fft of x[n]
vector<struct complex> FFT(vector<struct complex> &x){
    y = x;
    int N = x.size();
    for (int k = 0; k<N/2; k++)
        W.push_back(WNK(N, k));
    int m = int(log2(N));
    for (int stage = 1; stage<=m; stage++)
        x = butterfly(x,stage);
    x = inverse(x);
    return x;
}
// caculating the n points DFT of x[n]
vector<struct complex> DFT(vector<struct complex> &x){
    y = x;
    int N = x.size();
    for (int k = 0; k<N; k++) {
        struct complex tmp = {0,0};
        for (int n=0; n<N; n++)
            tmp = tmp + x[n]*WNK(N, n*k);
        y[k] = tmp;
    }
    return y;
}
int main(int argc, const char * argv[]) {
    cout<<"Please inter the points of DFT:"<<endl;
    int N;
    cin>>N;
    fstream file1("re.txt");
    fstream file2("im.txt");
    double re,im;
    for (int i = 0 ; i<N; i++){
        file1>>re;
        file2>>im;
        struct complex temp = {re,im};
        x.push_back(temp);
    }
    file1.close();
    file2.close();
    fstream ans1("ans_re.txt");
    fstream ans2("ans_im.txt");
    clock_t starttime,endtime;
    starttime = clock();
    x = FFT(x);
    endtime = clock();
    cout << "The run time is: " <<(double)(endtime - starttime)/ CLOCKS_PER_SEC << "s" << endl;
    for (int i = 0; i<x.size(); i++){
        ans1<<x[i].re<<endl;
        ans2<<x[i].im<<endl;
    }
    cout<<x.size()<<endl;
    ans1.close();
    ans2.close();
}
