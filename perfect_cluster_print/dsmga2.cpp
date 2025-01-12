/***************************************************************************
 *   Copyright (C) 2015 Tian-Li Yu and Shih-Huan Hsu                       *
 *   tianliyu@ntu.edu.tw                                                   *
 ***************************************************************************/

#include <list>
#include <vector>
#include <algorithm>
#include <iterator>
#include <random>
#include <iostream>
#include "chromosome.h"
#include "dsmga2.h"
#include "fastcounting.h"
#include "statistics.h"

#include <iomanip>
using namespace std;


DSMGA2::DSMGA2 (int n_ell, int n_nInitial, int n_maxGen, int n_maxFe, int fffff) {

    
    previousFitnessMean = -INF;
    ell = n_ell;
    nCurrent = (n_nInitial/2)*2;  // has to be even

    Chromosome::function = (Chromosome::Function)fffff;
    Chromosome::nfe = 0;
    Chromosome::lsnfe = 0;
    Chromosome::hitnfe = 0;
    Chromosome::hit = false;

    selectionPressure = 2;
    maxGen = n_maxGen;
    maxFe = n_maxFe;

    graph.init(ell);
    graph_size.init(ell);
     
    bestIndex = -1;
    masks = new list<int>[ell];
    selectionIndex = new int[nCurrent];
    orderN = new int[nCurrent];
    orderELL = new int[ell];
    population = new Chromosome[nCurrent];
    copy_population = new Chromosome[nCurrent+2]; 
    temp_population = new Chromosome[nCurrent+2];
    fastCounting = new FastCounting[ell];

    for (int i = 0; i < ell; i++)
        fastCounting[i].init(nCurrent);


    pHash.clear();
    for (int i=0; i<nCurrent; ++i) {
        population[i].initR(ell);
        double f = population[i].getFitness();
        pHash[population[i].getKey()] = f;
    }

    if (GHC) {
        for (int i=0; i < nCurrent; i++)
            population[i].GHC();
    }
    // cout << "DSMGA2 constructor" << endl;
}


DSMGA2::~DSMGA2 () {
    delete []masks;
    delete []orderN;
    delete []orderELL;
    delete []selectionIndex;
    delete []population;
    delete []fastCounting;
}



bool DSMGA2::isSteadyState () {

    if (stFitness.getNumber () <= 0)
        return false;

    if (previousFitnessMean < stFitness.getMean ()) {
        previousFitnessMean = stFitness.getMean () + EPSILON;
        return false;
    }

    return true;
}

bool DSMGA2::converged() {
    if (stFitness.getMax() == lastMax &&
        stFitness.getMean() == lastMean &&
        stFitness.getMin() == lastMin)
        convergeCount++;
    else
        convergeCount = 0;

    lastMax = stFitness.getMax();
    lastMean = stFitness.getMean();
    lastMin = stFitness.getMin();
    return (convergeCount > 300) ? true : false;
}

int DSMGA2::doIt (bool output) {
    generation = 0;
    while (!shouldTerminate ()) {
        oneRun (output);
       #ifdef DEBUG
       cin.get();
        #endif
    }
    return generation;
}


void DSMGA2::oneRun (bool output) {
    cout << "########################## NEW ONERUN ###############################" << endl;
    bool trace_process = false;
    // cout << "Generation: " << generation << endl;
    if (trace_process){
    cout << "Generation: " << generation << endl;}


    int m = ell / 3; // 分成幾組
    int original_nCurrent = nCurrent;

    bool conv0 = false;
    bool conv1 = false;
    bool conv2 = false;

    vector<int> i_count100(m, 0);
    vector<int> i_count011(m, 0);
    vector<int> i_count010(m, 0);
    vector<int> i_count101(m, 0);
    vector<int> i_count001(m, 0);
    vector<int> i_count110(m, 0);
    
    if (original_nCurrent >= 10)
    {
     


        for (int group = 0; group < m; ++group)
        {

            for (int i = 0; i < nCurrent; i++) {
                if (population[i].getVal(group * 3) == 1 && population[i].getVal(group * 3 + 1) == 0 && population[i].getVal(group * 3 + 2) == 0)
                {
                    i_count100[group]+=1;
                }
                else if (population[i].getVal(group * 3) == 0 && population[i].getVal(group * 3 + 1) == 1 && population[i].getVal(group * 3 + 2) == 1)
                {
                    i_count011[group]+=1;
                }
                else if (population[i].getVal(group * 3) == 0 && population[i].getVal(group * 3 + 1) == 1 && population[i].getVal(group * 3 + 2) == 0)
                {
                    i_count010[group]+=1;
                }
                else if (population[i].getVal(group * 3) == 1 && population[i].getVal(group * 3 + 1) == 0 && population[i].getVal(group * 3 + 2) == 1)
                {
                    i_count101[group]+=1;
                }
                else if (population[i].getVal(group * 3) == 0 && population[i].getVal(group * 3 + 1) == 0 && population[i].getVal(group * 3 + 2) == 1)
                {
                    i_count001[group]+=1;
                }
                else if (population[i].getVal(group * 3) == 1 && population[i].getVal(group * 3 + 1) == 1 && population[i].getVal(group * 3 + 2) == 0)
                {
                    i_count110[group]+=1;
                }
            }

            cout << endl;
            if (i_count011[group] == nCurrent)
            {
                // conv0 = true;
                cout << "第 " << group << " BB收斂到次佳解 011" << endl;
                // cout << "no!!!" << endl;
            }else if (i_count101[group] == nCurrent)
            {
                // conv1 = true;
                cout << "第 " << group << " BB收斂到次佳解 101" << endl;
                // cout << "no!!!" << endl;
            }else if (i_count110[group] == nCurrent)
            {
                // conv2 = true;
                cout << "第 " << group << " BB收斂到次佳解 110" << endl;
                // cout << "no!!!" << endl;
            }else{
                cout << "沒有分群前" << endl;
                cout << "BB " << group << " 沒有收斂" << endl;
            }

        }  

        cout << "===========" << endl;

    }

    // i_count100.assign(m, 0);
    // i_count011.assign(m, 0);
    // i_count010.assign(m, 0);
    // i_count101.assign(m, 0);
    // i_count001.assign(m, 0);
    // i_count110.assign(m, 0);
    int nCurrent1 = nCurrent;
    int nCurrent2;
    // cout << "nCurrent1: " << nCurrent1 << endl;


    
    vector<vector<int>> groups(3);
   
    int count100_011 = 0;
    int count100 = 0;
    int count011 = 0;
    int count010_101 = 0;
    int count010 = 0;
    int count101 = 0;
    int count001_110 = 0;
    int count001 = 0;
    int count110 = 0;
    


    // 定義三個群組的容器
    for (int i = 0; i < nCurrent; ++i) 
    {
        count100_011 = 0;
        count010_101 = 0;
        count001_110 = 0;

        count100 = 0;
        count011 = 0;
        count010 = 0;
        count101 = 0;
        count001 = 0;
        count110 = 0;

        // 計算每個 pattern 的出現次數
        for (int group = 0; group < m; ++group) {
            if (population[i].getVal(group * 3) == 1 && population[i].getVal(group * 3 + 1) == 0 && population[i].getVal(group * 3 + 2) == 0)
            {
                count100+=1;
            }
            else if (population[i].getVal(group * 3) == 0 && population[i].getVal(group * 3 + 1) == 1 && population[i].getVal(group * 3 + 2) == 1)
            {
                count011+=1;
            }
            else if (population[i].getVal(group * 3) == 0 && population[i].getVal(group * 3 + 1) == 1 && population[i].getVal(group * 3 + 2) == 0)
            {
                count010+=1;
            }
            else if (population[i].getVal(group * 3) == 1 && population[i].getVal(group * 3 + 1) == 0 && population[i].getVal(group * 3 + 2) == 1)
            {
                count101+=1;
            }
            else if (population[i].getVal(group * 3) == 0 && population[i].getVal(group * 3 + 1) == 0 && population[i].getVal(group * 3 + 2) == 1)
            {
                count001+=1;
            }
            else if (population[i].getVal(group * 3) == 1 && population[i].getVal(group * 3 + 1) == 1 && population[i].getVal(group * 3 + 2) == 0)
            {
                count110+=1;
            }
        }
        count100_011 = count100 + count011;
        count010_101 = count010 + count101;
        count001_110 = count001 + count110;

        if (trace_process) {
        cout << "count100 = "  << count100 << endl;
        cout << "count011 = "  << count011 << endl;
        cout << "count010 = "  << count010 << endl;
        cout << "count101 = "  << count101 << endl;
        cout << "count001 = "  << count001 << endl;
        cout << "count110 = "  << count110 << endl;

        cout << "count100_011 = "  << count100_011 << endl;
        cout << "count010_101 = "  << count010_101 << endl;
        cout << "count001_110 = "  << count001_110 << endl;
        }

        // 根據最多 pattern 分群
        int maxCount = std::max({count100_011, count010_101, count001_110});
        if (trace_process){
        cout << "maxCount: " << maxCount << endl;
        }


        if (maxCount == count100_011) {
            if (trace_process){
            cout << "groups[0].push_back index = "<< i << endl;
            for (int k = 0; k < ell; k++) {
                cout << population[i].getVal(k);
            }
            cout << endl;
            }

            groups[0].push_back(i);
        } else if (maxCount == count010_101) {
            if (trace_process){
            cout << "groups[1].push_back index = "<< i << endl;
            for (int k = 0; k < ell; k++) {
                cout << population[i].getVal(k);
            }
            cout << endl;
            }
            // cout << "count010_101: " << count010_101 << endl;
            groups[1].push_back(i);
        } else if (maxCount == count001_110){
            if (trace_process){
            cout << "groups[2].push_back index = "<< i << endl;
            for (int k = 0; k < ell; k++) {
                cout << population[i].getVal(k);
            }
            cout << endl;
            }
            // cout << "count001_110: " << count001_110 << endl;
            groups[2].push_back(i);
        }


        if (trace_process){

       

        // Output the indices of the three groups
        cout << "Group 0: ";
        for (int idx : groups[0]) {
            cout << idx << " ";
        }
        cout << endl;

        cout << "Group 1: ";
        for (int idx : groups[1]) {
            cout << idx << " ";
        }
        cout << endl;

        cout << "Group 2: ";
        for (int idx : groups[2]) {
            cout << idx << " ";
        }
        cout << endl;  

        }


    }

    if (generation == 0)
    {
        cout << endl;
        cout << "Gen = " << generation <<" & GHC+分群後"<< endl; 
        
    }else{
        cout << "Gen = " << generation << " 分群後 & 還沒執行 三個G執行 RM+BM 前"<< endl; 
    }

    conv0 = false;
    conv1 = false;
    conv2 = false;

    i_count100.assign(m, 0);
    i_count011.assign(m, 0);
    i_count010.assign(m, 0);
    i_count101.assign(m, 0);
    i_count001.assign(m, 0);
    i_count110.assign(m, 0);
    
    if (original_nCurrent >= 10)
    {
        // cout << "========================"<< endl; 
        cout<<"n0 = " << groups[0].size() << " n1 = " << groups[1].size() << " n2 = " << groups[2].size() << endl;


        for (int group = 0; group < m; ++group) {

            for (int i = 0; i < groups[0].size(); i++) {
                if (population[groups[0][i]].getVal(group * 3) == 1 && population[groups[0][i]].getVal(group * 3 + 1) == 0 && population[groups[0][i]].getVal(group * 3 + 2) == 0)
                {
                    i_count100[group]+=1;
                }
                else if (population[groups[0][i]].getVal(group * 3) == 0 && population[groups[0][i]].getVal(group * 3 + 1) == 1 && population[groups[0][i]].getVal(group * 3 + 2) == 1)
                {
                    i_count011[group]+=1;
                }
            }

            if (i_count011[group] == groups[0].size())
            {
                conv0 = true;
                cout << "G0 在 Gen" << generation << " 的第 " << group << " BB收斂到次佳解" << endl;
            }
        }  

        cout << "===========" << endl;

    

        for (int group = 0; group < m; ++group) {

            for (int i = 0; i < groups[1].size(); i++) {
                if (population[groups[1][i]].getVal(group * 3) == 0 && population[groups[1][i]].getVal(group * 3 + 1) == 1 && population[groups[1][i]].getVal(group * 3 + 2) == 0)
                {
                    i_count010[group]+=1;
                }
                else if (population[groups[1][i]].getVal(group * 3) == 1 && population[groups[1][i]].getVal(group * 3 + 1) == 0 && population[groups[1][i]].getVal(group * 3 + 2) == 1)
                {
                    i_count101[group]+=1;
                }
            }

            if (i_count101[group] == groups[1].size())
            {
                conv1 = true;
                cout << "G1 在 Gen" << generation << " 的第 " << group << " BB收斂到次佳解" << endl;
                


                // for (int i = 0; i < groups[1].size(); i++) {
                //     for (int j = 0; j < ell; j++) {
                //         cout << population[groups[1][i]].getVal(j);
                //     }
                //     cout << endl;        
                // }  


            }
        }
        cout << "===========" << endl;
  

        for (int group = 0; group < m; ++group) {

            for (int i = 0; i < groups[2].size(); i++) {
                if (population[groups[2][i]].getVal(group * 3) == 0 && population[groups[2][i]].getVal(group * 3 + 1) == 0 && population[groups[2][i]].getVal(group * 3 + 2) == 1)
                {
                    i_count001[group]+=1;
                }
                else if (population[groups[2][i]].getVal(group * 3) == 1 && population[groups[2][i]].getVal(group * 3 + 1) == 1 && population[groups[2][i]].getVal(group * 3 + 2) == 0)
                {
                    i_count110[group]+=1;
                }
            }

            if (i_count110[group] == groups[2].size())
            {
                conv2 = true;
                cout << "G2 在 Gen" << generation << " 的第 " << group << " BB收斂到次佳解" << endl;
               


                // for (int i = 0; i < groups[2].size(); i++) {
                //     for (int j = 0; j < ell; j++) {
                //         cout << population[groups[2][i]].getVal(j);
                //     }
                //     cout << endl;        
                // }  
            }
            
        }

        cout << "===========" << endl;

        if (conv0 && conv1 && conv2)
        {
            cout << "Gen = " << generation << " 三群都有BB收斂到次佳解" << "n=" << original_nCurrent<< endl;
            cout << "這是在分群好，還沒執行 RM+BM" << endl;
            // cout << "=======Pop G0========"<< endl;

            // for (int i = 0; i < groups[0].size(); i++) {
            //     for (int j = 0; j < ell; j++) {
            //         cout << population[groups[0][i]].getVal(j);
            //     }
            //     cout << endl;        
            // }  


            // cout << "=======Pop G1========"<< endl;

            // for (int i = 0; i < groups[1].size(); i++) {
            //     for (int j = 0; j < ell; j++) {
            //         cout << population[groups[1][i]].getVal(j);
            //     }
            //     cout << endl;        
            // } 

            // cout << "=======Pop G2========"<< endl;
            // for (int i = 0; i < groups[2].size(); i++) {
            //     for (int j = 0; j < ell; j++) {
            //         cout << population[groups[2][i]].getVal(j);
            //     }
            //     cout << endl;        
            // } 
            // generation = maxGen;
        }
       
    

    }

    for (int i = 0; i < nCurrent; ++i) {
        copy_population[i] = population[i];
    }   

    // 分好群 >>> 執行獨立 RM + BM
    int i_count = 0;
    for (int G = 0; G < 3; ++G) {
        nCurrent = groups[G].size();

        // if (nCurrent % 2 == 1){
        //      // Add a chromosome to this group
        //      population[nCurrent] = population[groups[G][0]];
        //      groups[G].push_back(nCurrent);
        //      nCurrent++;
        //  }
           
        
        for (int tt = 0; tt < groups[G].size(); ++tt) {
            // cout << "groups[G][tt] = " << groups[G][tt] << endl;
            population[tt] = copy_population[groups[G][tt]]; 
        }

        // for (int tt = 0; tt < groups[G].size(); ++tt) {
        //     population[tt] = copy_population[tt];
        // }


        if (CACHE)
            Chromosome::cache.clear();

        if (trace_process){
        cout << "++++++++++++++++++++++Mixing++++++++++++++++++++++" << endl;
        


        cout << "Group " << G << ":" << endl;   
        for (int tt = 0; tt < nCurrent; ++tt) {
            for (int ttt = 0; ttt < ell; ++ttt) {
                cout << population[tt].getVal(ttt);
            }
            cout << endl;
        }       
       
        }

        
        
        // cout <<"還沒執行RM+BM的 G" << G << "pop" << endl;
        // cout << "G" << G << "size = " <<groups[G].size()<< endl;
        // for (int i = 0; i < groups[G].size(); i++) {
        //     for (int j = 0; j < ell; j++) {
        //         cout << population[i].getVal(j);
        //     }
        //     cout <<endl;        
        // }  
        // cout << "-----" << endl; 
        
        
        // cout <<  endl;
        cout << "G = " << G << " nCurrent = " << nCurrent << endl; 
        cout << "mixing start G = " << G<< endl;
        mixing();
        cout << "mixing end G = " << G<< endl;

       
        // cout <<"執行 RM+BM 後的 G" << G << "pop" << endl;
        // cout << "G" << G << "size = " <<groups[G].size()<< endl;
        // for (int i = 0; i < groups[G].size(); i++) {
        //     for (int j = 0; j < ell; j++) {
        //         cout << population[i].getVal(j);
        //     }
        //     cout <<endl;        
        // }  
        // cout << "-----" << endl; 
    



        i_count100.assign(m, 0);
        i_count011.assign(m, 0);
        i_count010.assign(m, 0);
        i_count101.assign(m, 0);
        i_count001.assign(m, 0);
        i_count110.assign(m, 0);

        if (G == 0){

            conv0 = false;
   


            for (int group = 0; group < m; ++group) {

                for (int current_i = 0; current_i < groups[0].size(); current_i++) {
                    if (population[current_i].getVal(group * 3) == 1 && population[current_i].getVal(group * 3 + 1) == 0 && population[current_i].getVal(group * 3 + 2) == 0)
                    {
                        i_count100[group]+=1;
                    }
                    else if (population[current_i].getVal(group * 3) == 0 && population[current_i].getVal(group * 3 + 1) == 1 && population[current_i].getVal(group * 3 + 2) == 1)
                    {
                        i_count011[group]+=1;
                    }
                }

                if (i_count011[group] == groups[0].size())
                {
                    conv0 = true;
                    cout << "G0 在 Gen" << generation << " 的第 " << group << " BB收斂到次佳解" << endl;
                    
                }
            }

        }else if (G == 1){
            
            conv1 = false;
        

            for (int group = 0; group < m; ++group) {

                for (int current_i = 0; current_i < groups[1].size(); current_i++) {
                    if (population[current_i].getVal(group * 3) == 0 && population[current_i].getVal(group * 3 + 1) == 1 && population[current_i].getVal(group * 3 + 2) == 0)
                    {
                        i_count010[group]+=1;
                    }
                    else if (population[current_i].getVal(group * 3) == 1 && population[current_i].getVal(group * 3 + 1) == 0 && population[current_i].getVal(group * 3 + 2) == 1)
                    {
                        i_count101[group]+=1;
                    }
                }

                if (i_count101[group] == groups[1].size())
                {
                    conv1 = true;
                    cout << "G1 在 Gen" << generation << " 的第 " << group << " BB收斂到次佳解" << endl;

                    // for (int i = 0; i < groups[1].size(); i++) {
                    //     for (int j = 0; j < ell; j++) {
                    //         cout << population[i].getVal(j);
                    //     }
                    //     cout << endl;        
                    // }  
                }
            }


        }else if (G == 2){
            
            conv2 = false;
            for (int group = 0; group < m; ++group) {

                for (int current_i = 0; current_i < groups[2].size(); current_i++) {
                    if (population[current_i].getVal(group * 3) == 0 && population[current_i].getVal(group * 3 + 1) == 0 && population[current_i].getVal(group * 3 + 2) == 1)
                    {
                        i_count001[group]+=1;
                    }
                    else if (population[current_i].getVal(group * 3) == 1 && population[current_i].getVal(group * 3 + 1) == 1 && population[current_i].getVal(group * 3 + 2) == 0)
                    {
                        i_count110[group]+=1;
                    }
                }

                if (i_count110[group] == groups[2].size())
                {
                    conv2 = true;
                    cout << "G2 在 Gen" << generation << " 的第 " << group << " BB收斂到次佳解" << endl;

                    // for (int i = 0; i < groups[2].size(); i++) {
                    //     for (int j = 0; j < ell; j++) {
                    //         cout << population[i].getVal(j);
                    //     }
                    //     cout << endl;        
                    // }  
                }
                
            }
        }

        if (G == 2 ){
            if (conv0 && conv1 && conv2)
            {
                cout << "Gen = " << generation << " After BM 三群都有BB收斂到次佳解!!" << "n=" << original_nCurrent<< endl;
                // cout << "Here" << endl;
                // cout << "=======Pop G0========"<< endl;

                // for (int i = 0; i < groups[0].size(); i++) {
                //     for (int j = 0; j < ell; j++) {
                //         cout << population[groups[0][i]].getVal(j);
                //     }
                //     cout << endl;        
                // }
            }else{
                cout << "Gen = " << generation << " After BM 三群沒有都有BB收斂到次佳解" << "n=" << original_nCurrent <<endl;
                if(conv0){
                    cout << "存在 G0 BB收斂到次佳解" << endl;
                }
                if(conv1){
                    cout << "存在 G1 BB收斂到次佳解" << endl;
                }
                if(conv2){
                    cout << "存在 G2 BB收斂到次佳解" << endl;
                }
            }
        }

        cout << "xxxxxxxxxxxxxxxxxxxxx" <<"所有 mixing 結束" << "xxxxxxxxxxxxxxxxxxxxx" << endl;
        
    
        for (int ii = 0; ii < nCurrent; ++ii) {
            if (ii+i_count >= original_nCurrent) break;
            temp_population[ii+i_count] = population[ii];
        }
        i_count += nCurrent;

        if (G == 0) {
            nCurrentA = nCurrent;
        } else if (G == 1) {
            nCurrentB = nCurrent;
        } else {
            nCurrentC = nCurrent;
        }

    }

    nCurrent = original_nCurrent;
    // 將群組中的染色體放回 population
    for (int i = 0; i < nCurrent; ++i) {
        population[i] = temp_population[i];
    }



    // // Output the indices of the three groups
    // cout << "Group 0: ";
    // for (int idx : groups[0]) {
    //     cout << idx << " ";
    // }
    // cout << endl;

    // cout << "Group 2: ";
    // for (int idx : groups[1]) {
    //     cout << idx << " ";
    // }
    // cout << endl;

    // cout << "Group 2: ";
    // for (int idx : groups[2]) {
    //     cout << idx << " ";
    // }
    // cout << endl;
    // cout <<"!!!!!!!!!" <<endl;    


    // nCurrent =  groups[0].size() + groups[1].size() + groups[2].size();
    // nCurrent2 = nCurrent;
    // if (nCurrent2 - nCurrent1 > 2){
    //     cout << "nCurrent2 - nCurrent1: " << nCurrent2 - nCurrent1 << endl;
    //     // cout << "nCurrent1: " << nCurrent1 << endl;
    //     // cout << "nCurrent2: " << nCurrent2 << endl;
    // }
    

    // if (nCurrent - original_nCurrent == 2)
    // {
    //     // cout << "nCurrent - original_nCurrent = " << nCurrent - original_nCurrent << endl;
    //     if (temp_population[nCurrent-1].getFitness()>=temp_population[nCurrent-2].getFitness()&& temp_population[nCurrent-1].getFitness()>=temp_population[nCurrent-3].getFitness())
    //     {
    //         population[original_nCurrent-1] = temp_population[nCurrent-1];
    //         nCurrent = original_nCurrent;
    //         for (int i = 0; i < original_nCurrent-1; ++i) {
    //             population[i] = temp_population[i];
    //         }
    //     }else if (temp_population[nCurrent-2].getFitness()>=temp_population[nCurrent-1].getFitness()&& temp_population[nCurrent-2].getFitness()>=temp_population[nCurrent-3].getFitness())
    //     {
    //         population[original_nCurrent-1] = temp_population[nCurrent-2];
    //         nCurrent = original_nCurrent;
    //         for (int i = 0; i < original_nCurrent-1; ++i) {
    //             population[i] = temp_population[i];
    //         }
    //     }else if (temp_population[nCurrent-3].getFitness()>=temp_population[nCurrent-1].getFitness()&& temp_population[nCurrent-3].getFitness()>=temp_population[nCurrent-2].getFitness())
    //     {
    //         population[original_nCurrent-1] = temp_population[nCurrent-3];
    //         nCurrent = original_nCurrent;
    //         for (int i = 0; i < original_nCurrent-1; ++i) {
    //             population[i] = temp_population[i];
    //         }
    //     }else{
    //         nCurrent = original_nCurrent;
    //         // 將群組中的染色體放回 population
    //         for (int i = 0; i < nCurrent; ++i) {
    //             population[i] = temp_population[i];
    //         }            
    //     }
    // }
    // else{
    //         nCurrent = original_nCurrent;
    //         // 將群組中的染色體放回 population
    //         for (int i = 0; i < nCurrent; ++i) {
    //             population[i] = temp_population[i];
    //          }
    // }
    


    // nCurrent = original_nCurrent;
    // // 將群組中的染色體放回 population
    // for (int i = 0; i < nCurrent; ++i) {
    //     population[i] = temp_population[i];
    // }


    double max = -INF;
    stFitness.reset ();

    for (int i = 0; i < nCurrent; ++i) {
        double fitness = population[i].getFitness();
        if (fitness > max) {
            max = fitness;
            bestIndex = i;
        }
        stFitness.record (fitness);

    }

    if (output)
        showStatistics ();

    ++generation;
}


bool DSMGA2::shouldTerminate () {
    bool  termination = false;

    if (maxFe != -1) {
        if (Chromosome::nfe > maxFe)
            termination = true;
    }

    if (maxGen != -1) {
        if (generation > maxGen)
            termination = true;
    }

    if (population[0].getMaxFitness() <= stFitness.getMax() )
        termination = true;


    if (stFitness.getMax() - EPSILON <= stFitness.getMean() )
        termination = true;

    return termination;

}


bool DSMGA2::foundOptima () {
    return (stFitness.getMax() > population[0].getMaxFitness());
}


void DSMGA2::showStatistics () {

    printf ("Gen:%d  Fitness:(Max/Mean/Min):%f/%f/%f \n ",
            generation, stFitness.getMax (), stFitness.getMean (),
            stFitness.getMin ());
    fflush(NULL);
}



void DSMGA2::buildFastCounting() {

    if (SELECTION) {
        for (int i = 0; i < nCurrent; i++)
            for (int j = 0; j < ell; j++) {
                fastCounting[j].setVal(i, population[selectionIndex[i]].getVal(j));
            }

    } else {
        for (int i = 0; i < nCurrent; i++)
            for (int j = 0; j < ell; j++) {
                fastCounting[j].setVal(i, population[i].getVal(j));
            }
    }

}

int DSMGA2::countOne(int x) const {

    int n = 0;

    for (int i=0; i<fastCounting[0].lengthLong; ++i) {
        unsigned long val = 0;

        val = fastCounting[x].gene[i];

        n += myBD.countOne(val);
    }

    return n;
}


int DSMGA2::countXOR(int x, int y) const {

    int n = 0;

    for (int i=0; i<fastCounting[0].lengthLong; ++i) {
        unsigned long val = 0;


        val = fastCounting[x].gene[i];

        val ^= fastCounting[y].gene[i];

        n += myBD.countOne(val);
    }

    return n;
}

//2016-03-09
// Almost identical to DSMGA2::findClique
// except check 00 or 01 before adding connection
void DSMGA2::findMask(Chromosome& ch, list<int>& result,int startNode){
    result.clear();

    
	DLLA rest(ell);
	genOrderELL();
	for( int i = 0; i < ell; i++){
		if(orderELL[i] == startNode)
			result.push_back(orderELL[i]);
		else
			rest.insert(orderELL[i]);
	}

	double *connection = new double[ell];

	for(DLLA::iterator iter = rest.begin(); iter != rest.end(); iter++){
	    pair<double, double> p = graph(startNode, *iter);
		int i = ch.getVal(startNode);
		int j = ch.getVal(*iter);
		if(i == j)//p00 or p11
			connection[*iter] = p.first;
		else      //p01 or p10
			connection[*iter] = p.second;
	}
   
    while(!rest.isEmpty()){

        // cout << "A" << endl;

	    double max = -INF;
		int index = -1;
		for(DLLA::iterator iter = rest.begin(); iter != rest.end(); iter++){
		    if(max < connection[*iter]){
                // cout << "B" << endl;
			    max = connection[*iter];
				index = *iter;
			}
		}

		rest.erase(index);
		result.push_back(index);

		for(DLLA::iterator iter = rest.begin(); iter != rest.end(); iter++){

            if(index == -1)
            {
                break;
            }
            
            // cout << "index= " <<index<< endl;
			pair<double, double> p = graph(index, *iter);
			int i = ch.getVal(index);
			int j = ch.getVal(*iter);

			if(i == j)//p00 or p11
            {
                // cout << "C" << endl;
				connection[*iter] += p.first;
                // cout << "D" << endl;
            }
			else      //p01 or p10
            {
                // cout << "E" << endl;
				connection[*iter] += p.second;
                // cout << "F" << endl;
            }

            // if (*iter >= 0 && *iter < ell) {
            //     if (i == j) // p00 or p11
            //         connection[*iter] += p.first;
            //     else        // p01 or p10
            //         connection[*iter] += p.second;
            // } else {
            //     std::cerr << "Error: *iter out of bounds! Value: " << *iter << std::endl;
            // }
		}

        if(index == -1)
        {
            break;
        }
	}



	delete []connection;
  
}

void DSMGA2::restrictedMixing(Chromosome& ch) {
    
    int startNode = myRand.uniformInt(0, ell - 1);    
    


    list<int> mask;
	findMask(ch, mask,startNode);
    size_t size = findSize(ch, mask);
   
    list<int> mask_size; 
    findMask_size(ch,mask_size,startNode,size);
    size_t size_original = findSize(ch,mask_size);

    if (size > size_original)
        size = size_original;
    while (mask.size() > size)
        mask.pop_back();
   
    bool taken = restrictedMixing(ch, mask);


    EQ = true;
    if (taken) {
    
        genOrderN();

        // cout << "The Donor: " << endl;  
        // for (int i=0; i<ell; ++i) {
        //     cout << ch.getVal(i);
        // }
        // cout<< endl;       


        // cout << "After RM pop: " << endl;
        // for (int _index=0; _index<nCurrent; ++_index) {
        //     for (int i=0; i<ell; ++i) {
        //         cout << population[orderN[_index]].getVal(i);
        //     }
        //     cout<< endl;
        // }


        // cout << "Donor" << endl;
        // for (auto it = mask.begin(); it != mask.end(); ++it) {
        //     cout << "("<< *it << " ";
        //     cout<< ch.getVal(*it) << ")"; 
        // }
        // cout << endl;

        for (int i=0; i<nCurrent; ++i) {

            if (EQ)
                backMixingE(ch, mask, population[orderN[i]]);
            else
                backMixing(ch, mask, population[orderN[i]]);
        }

        //  cout << "After BM pop:" << endl;
        // for (int i=0; i<nCurrent; ++i) {
        //     for (int j=0; j<ell; ++j) {
        //         cout << population[i].getVal(j);
        //     }
        // cout << endl;
        // }

    }

}
void DSMGA2::findMask_size(Chromosome& ch, list<int>& result,int startNode,int bound){
    result.clear();

    
	DLLA rest(ell);

	for( int i = 0; i < ell; i++){
		if(orderELL[i] == startNode)
			result.push_back(orderELL[i]);
		else
			rest.insert(orderELL[i]);
	}

	double *connection = new double[ell];

	for(DLLA::iterator iter = rest.begin(); iter != rest.end(); iter++){
	    pair<double, double> p = graph_size(startNode, *iter);
		int i = ch.getVal(startNode);
		int j = ch.getVal(*iter);
		if(i == j)//p00 or p11
			connection[*iter] = p.first;
		else      //p01 or p10
			connection[*iter] = p.second;
	}
    bound--;
    while(!rest.isEmpty()&&bound>0){
        bound--;
	    double max = -INF;
		int index = -1;
		for(DLLA::iterator iter = rest.begin(); iter != rest.end(); iter++){
		    if(max < connection[*iter]){
			    max = connection[*iter];
				index = *iter;
			}
		}

		rest.erase(index);
		result.push_back(index);

		for(DLLA::iterator iter = rest.begin(); iter != rest.end(); iter++){
			pair<double, double> p = graph_size(index, *iter);
			int i = ch.getVal(index);
			int j = ch.getVal(*iter);
			if(i == j)//p00 or p11
				connection[*iter] += p.first;
			else      //p01 or p10
				connection[*iter] += p.second;
		}
	}

    delete []connection;
  
}

void DSMGA2::backMixing(Chromosome& source, list<int>& mask, Chromosome& des) {
    // cout << "Hey!!" << endl;
    Chromosome trial(ell);
    trial = des;
    for (list<int>::iterator it = mask.begin(); it != mask.end(); ++it)
        trial.setVal(*it, source.getVal(*it));

    if (trial.getFitness() > des.getFitness()) {
        cout << "BM 成功了" << endl;
        pHash.erase(des.getKey());
        pHash[trial.getKey()] = trial.getFitness();
        des = trial;
          
        return;
    }

}

void DSMGA2::backMixingE(Chromosome& source, list<int>& mask, Chromosome& des) {

    Chromosome trial(ell);
    trial = des;
    for (list<int>::iterator it = mask.begin(); it != mask.end(); ++it)
        trial.setVal(*it, source.getVal(*it));

    if (trial.getFitness() > des.getFitness()) {
        cout << "BM 成功了" << endl;
        pHash.erase(des.getKey());
        pHash[trial.getKey()] = trial.getFitness();

        EQ = false;
        des = trial;

    return;
    }

    //2016-10-21
    if (trial.getFitness() >= des.getFitness() - EPSILON) {
        pHash.erase(des.getKey());
        pHash[trial.getKey()] = trial.getFitness();

        des = trial;
        return;
    }

}

bool DSMGA2::restrictedMixing(Chromosome& ch, list<int>& mask) {

    bool taken = false;
    size_t lastUB = 0;

    for (size_t ub = 1; ub <= mask.size(); ++ub) {

        size_t size = 1;
        Chromosome trial(ell);
        trial = ch;
	    
		//2016-03-03
	    vector<int> takenMask;

        for (list<int>::iterator it = mask.begin(); it != mask.end(); ++it) {
            
            //2016-03-03
			takenMask.push_back(*it);

            trial.flip(*it);

            ++size;
            if (size > ub) break;
        }

        //if (isInP(trial)) continue;
        //2016-10-21
        if (isInP(trial)) break;

        if (trial.getFitness() >= ch.getFitness() - EPSILON) {
            pHash.erase(ch.getKey());
            pHash[trial.getKey()] = trial.getFitness();

            taken = true;
            ch = trial;
        }

        if (taken) {
            lastUB = ub;
            break;
        }
    }

    if (lastUB != 0) {
        while (mask.size() > lastUB)
            mask.pop_back();
    }

    return taken;

}

size_t DSMGA2::findSize(Chromosome& ch, list<int>& mask) const {

    DLLA candidate(nCurrent);
    for (int i=0; i<nCurrent; ++i)
        candidate.insert(i);

    size_t size = 0;
    for (list<int>::iterator it = mask.begin(); it != mask.end(); ++it) {

        int allele = ch.getVal(*it);

        for (DLLA::iterator it2 = candidate.begin(); it2 != candidate.end(); ++it2) {
            if (population[*it2].getVal(*it) == allele)
                candidate.erase(*it2);

            if (candidate.isEmpty())
                break;
        }

        if (candidate.isEmpty())
            break;

        ++size;
    }

    return size;


}

size_t DSMGA2::findSize(Chromosome& ch, list<int>& mask, Chromosome& ch2) const {

    size_t size = 0;
    for (list<int>::iterator it = mask.begin(); it != mask.end(); ++it) {
        if (ch.getVal(*it) == ch2.getVal(*it)) break;
        ++size;
    }
    return size;
}

void DSMGA2::mixing() {

    

    if (SELECTION)
        selection();

    //* really learn model
    buildFastCounting();
    buildGraph();
    buildGraph_sizecheck();
    for (int i=0; i<ell; ++i)
       findClique(i, masks[i]); // replaced by findMask in restrictedMixing

    int repeat = (ell>50)? ell/50: 1;

    // cout << "Population: before restrictedMixing(population[orderN[i]])" << endl;
    // for (int i=0; i<nCurrent; ++i) {
    //     for (int j=0; j<ell; ++j) {
    //         cout << population[i].getVal(j);
    //     }
    //     cout << endl;
    // }





    for (int k=0; k<repeat; ++k) {

        genOrderN();
        for (int i=0; i<nCurrent; ++i) {
            restrictedMixing(population[orderN[i]]);
            if (Chromosome::hit) break;
        }
        if (Chromosome::hit) break;
    }



    // cout << "mixing end" << endl << endl;


}

inline bool DSMGA2::isInP(const Chromosome& ch) const {

    unordered_map<unsigned long, double>::const_iterator it = pHash.find(ch.getKey());
    return (it != pHash.end());
}

inline void DSMGA2::genOrderN() {
    myRand.uniformArray(orderN, nCurrent, 0, nCurrent-1);
}

inline void DSMGA2::genOrderELL() {
    myRand.uniformArray(orderELL, ell, 0, ell-1);
}

void DSMGA2::buildGraph() {

    int *one = new int [ell];
    for (int i=0; i<ell; ++i) {
        one[i] = countOne(i);
    }

    for (int i=0; i<ell; ++i) {

        for (int j=i+1; j<ell; ++j) {

            int n00, n01, n10, n11;
            int nX =  countXOR(i, j);

            n11 = (one[i]+one[j]-nX)/2;
            n10 = one[i] - n11;
            n01 = one[j] - n11;
            n00 = nCurrent - n01 - n10 - n11;

            double p00 = (double)n00/(double)nCurrent;
            double p01 = (double)n01/(double)nCurrent;
            double p10 = (double)n10/(double)nCurrent;
            double p11 = (double)n11/(double)nCurrent;
            double p1_ = p10 + p11;
            double p0_ = p00 + p01;
            double p_0 = p00 + p10;
            double p_1 = p01 + p11;

            double linkage = computeMI(p00,p01,p10,p11);
            
            //2016-04-08_computeMI_entropy
            double linkage00 = 0.0, linkage01 = 0.0;
            if (p00 > EPSILON)
                linkage00 += p00*log(p00/p_0/p0_);
            if (p11 > EPSILON)
                linkage00 += p11*log(p11/p_1/p1_);
            if (p01 > EPSILON)
                linkage01 += p01*log(p01/p0_/p_1);
            if (p10 > EPSILON)
                linkage01 += p10*log(p10/p1_/p_0);
           
            if(Chromosome::nfe < 0){
                pair<double, double> p(linkage, linkage);
                graph.write(i, j, p);
            }
            else{
                pair<double, double> p(linkage00, linkage01);
                graph.write(i, j, p);
            }
				
        }
    }


    delete []one;

}
void DSMGA2::buildGraph_sizecheck() {

    int *one = new int [ell];
    for (int i=0; i<ell; ++i) {
        one[i] = countOne(i);
    }

    for (int i=0; i<ell; ++i) {

        for (int j=i+1; j<ell; ++j) {

            int n00, n01, n10, n11;
            int nX =  countXOR(i, j);

            n11 = (one[i]+one[j]-nX)/2;
            n10 = one[i] - n11;
            n01 = one[j] - n11;
            n00 = nCurrent - n01 - n10 - n11;

            double p00 = (double)n00/(double)nCurrent;
            double p01 = (double)n01/(double)nCurrent;
            double p10 = (double)n10/(double)nCurrent;
            double p11 = (double)n11/(double)nCurrent;
            double p1_ = p10 + p11;
            double p0_ = p00 + p01;
            double p_0 = p00 + p10;
            double p_1 = p01 + p11;

            double linkage = computeMI(p00,p01,p10,p11);
            
            //2016-04-08_computeMI_entropy
            double linkage00 = 0.0, linkage01 = 0.0;
            if (p00 > EPSILON)
                linkage00 += p00*log(p00/p_0/p0_);
            if (p11 > EPSILON)
                linkage00 += p11*log(p11/p_1/p1_);
            if (p01 > EPSILON)
                linkage01 += p01*log(p01/p0_/p_1);
            if (p10 > EPSILON)
                linkage01 += p10*log(p10/p1_/p_0);
        
	
            pair<double, double> p(linkage, linkage);
            graph_size.write(i, j, p);
			
        }
    }


    delete []one;

}


// from 1 to ell, pick by max edge
void DSMGA2::findClique(int startNode, list<int>& result) {
    
   }

    
    double DSMGA2::computeMI(double a00, double a01, double a10, double a11) const {

    double p0 = a00 + a01;
    double q0 = a00 + a10;
    double p1 = 1-p0;
    double q1 = 1-q0;

    double join = 0.0;
    if (a00 > EPSILON)
        join += a00*log(a00);
    if (a01 > EPSILON)
        join += a01*log(a01);
    if (a10 > EPSILON)
        join += a10*log(a10);
    if (a11 > EPSILON)
        join += a11*log(a11);

    double p = 0.0;
    if (p0 > EPSILON)
        p += p0*log(p0);
    if (p1 > EPSILON)
        p += p1*log(p1);


    double q = 0.0;
    if (q0 > EPSILON)
        q += q0*log(q0);
    if (q1 > EPSILON)
        q += q1*log(q1);

    return -p-q+join;

}


void DSMGA2::selection () {
    tournamentSelection ();
}


// tournamentSelection without replacement
void DSMGA2::tournamentSelection () {
    int i, j;

    int randArray[selectionPressure * nCurrent];

    for (i = 0; i < selectionPressure; i++)
        myRand.uniformArray (randArray + (i * nCurrent), nCurrent, 0, nCurrent - 1);

    for (i = 0; i < nCurrent; i++) {

        int winner = 0;
        double winnerFitness = -INF;

        for (j = 0; j < selectionPressure; j++) {
            int challenger = randArray[selectionPressure * i + j];
            double challengerFitness = population[challenger].getFitness();

            if (challengerFitness > winnerFitness) {
                winner = challenger;
                winnerFitness = challengerFitness;
            }

        }
        selectionIndex[i] = winner;
    }
}


int DSMGA2::get_ch_dist(Chromosome x, Chromosome y) {
    
    int result(0);
    
    // for (int i=0; i!=ell; ++i) {
    //     if (x.getVal(i) != y.getVal(i))
    //         ++result;
    // }
    
    // assert(result >=0 && result <= ell);

    for (int i=0; i<ell; ++i) {
        result += myBD.countOne(x.getVal(i) ^ y.getVal(i));
    }
    
    return result;
}
