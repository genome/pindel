//
//  smith_waterman_alignment.cpp
//  
//  Taken from: https://wiki.uni-koeln.de/biologicalphysics/index.php/Implementation_of_the_Smith-Waterman_local_alignment_algorithm
//

#include <iostream>
#include <limits.h>
#include "smith_waterman_alignment.h"


const int MAX_MISMATCHES = 2;



const float mu = 0.33;
const float delta = 1.33;
int ind;


static double similarity_score(char a,char b){
    
    double result;
    if(a==b){
        result=1.;
    }
    else{
        result=-mu;
    }
    return result;
}


static double find_array_max(double array[],int length){
    
    double max = array[0];            // start with max = first element
    ind = 0;
    
    for(int i = 1; i<length; i++){
        if(array[i] > max){
            max = array[i];
            ind = i; 
        }
    }
    return max;                    // return highest value in array
}




// Do Smith-Waterman alignment for given sequences.
int get_alignment_pos(std::string seq_a, std::string seq_b, bool& viable_alignment) {
    // string s_a=seq_a,s_b=seq_b;
    int N_a = seq_a.length();                     // get the actual lengths of the sequences
    int N_b = seq_b.length();

    // initialize H
    double H[N_a+1][N_b+1];     
    for(int i=0;i<=N_a;i++){
        for(int j=0;j<=N_b;j++){
            H[i][j]=0.;
        }
    } 

    double temp[4];
    int I_i[N_a+1][N_b+1],I_j[N_a+1][N_b+1];     // Index matrices to remember the 'path' for backtracking


    for(int i=1;i<=N_a;i++){
        for(int j=1;j<=N_b;j++){
            temp[0] = H[i-1][j-1]+similarity_score(seq_a[i-1],seq_b[j-1]); 
            temp[1] = H[i-1][j]-delta;                  
            temp[2] = H[i][j-1]-delta;                 
            temp[3] = 0.;
            H[i][j] = find_array_max(temp,4);
            switch(ind){
                case 0:                                  // score in (i,j) stems from a match/mismatch
                    I_i[i][j] = i-1;
                    I_j[i][j] = j-1;
                    break;
                case 1:                                  // score in (i,j) stems from a deletion in sequence A
                    I_i[i][j] = i-1;
                    I_j[i][j] = j;
                    break;
                case 2:                                  // score in (i,j) stems from a deletion in sequence B
                    I_i[i][j] = i;
                    I_j[i][j] = j-1;
                    break;
                case 3:                                  // (i,j) is the beginning of a subsequence
                    I_i[i][j] = i;
                    I_j[i][j] = j;	
                    break;
            }
        }
    }

    // Print the matrix H to the console
    /*
    std::cout<<"**********************************************"<<std::endl;
    std::cout<<"The scoring matrix is given by  "<<std::endl<<std::endl;
    for(int i=1;i<=N_a;i++){
        for(int j=1;j<=N_b;j++){
            std::cout<<H[i][j]<<" ";
        }
        std::cout<<std::endl;
    }
     */

    // search H for the maximal score
    double H_max = 0.;
    int i_max=0,j_max=0;
    for(int i=1;i<=N_a;i++){
        for(int j=1;j<=N_b;j++){
            if(H[i][j]>H_max){
                H_max = H[i][j];
                i_max = i;
                j_max = j;
            }
        }
    }

    //cout<<H_max<<endl;
//    std::cout << N_a << ", " << N_b << ", " << i_max << ", " << j_max << std::endl;
//    std::cout << "max score: " << H_max << std::endl;
//    std::cout << "alignment index: " << (i_max - j_max) << std::endl;

    
    // Backtracking from H_max
    int current_i=i_max,current_j=j_max;
    int next_i=I_i[current_i][current_j];
    int next_j=I_j[current_i][current_j];
    int tick=0;
    char consensus_a[N_a+N_b+2],consensus_b[N_a+N_b+2];

    while(((current_i!=next_i) || (current_j!=next_j)) && (next_j>=0) && (next_i>=0)){
    
        if(next_i==current_i)  consensus_a[tick] = '-';                  // deletion in A
        else                   consensus_a[tick] = seq_a[current_i-1];   // match/mismatch in A
    
        if(next_j==current_j)  consensus_b[tick] = '-';                  // deletion in B
        else                   consensus_b[tick] = seq_b[current_j-1];   // match/mismatch in B
    
        // mkroon: fix for adding first character of the alignment.
        if (next_i == 0) {
            next_i = -1;
        } else if (next_j == 0) {
            next_j = -1;
        } else {
            current_i = next_i;
            current_j = next_j;
            next_i = I_i[current_i][current_j];
            next_j = I_j[current_i][current_j];
        }
        tick++;
    }

    /*
    // Output of the consensus motif to the console
    std::cout<<std::endl<<"***********************************************"<<std::endl;
    std::cout<<"The alignment of the sequences"<<std::endl<<std::endl;
    for(int i=0;i<N_a;i++){std::cout<<seq_a[i];}; std::cout<<"  and"<<std::endl;
    for(int i=0;i<N_b;i++){std::cout<<seq_b[i];}; std::cout<<std::endl<<std::endl;
    std::cout<<"is for the parameters  mu = "<<mu<<" and delta = "<<delta<<" given by"<<std::endl<<std::endl;  
    for(int i=tick-1;i>=0;i--) std::cout<<consensus_a[i]; 
    std::cout<<std::endl;
    for(int j=tick-1;j>=0;j--) std::cout<<consensus_b[j];
    std::cout<<std::endl;
    */
    
    int num_mismatches = 0;
    for (int i=tick-1; i>=0; i--) {
        if (consensus_a[i] != consensus_b[i]) {
            num_mismatches++;
        }
    }

    viable_alignment = num_mismatches <= MAX_MISMATCHES;
    return i_max - j_max;
}

