//2. Estimate mobility, then measure the correlation function in both x and y axis.
//1. Introduce mobility (cells can exchange with neighbors in x axis) to reduce rho_y to 0.17


#include <iostream>
#include <utility> 

#include <string>
#include <list>
#include <queue>
#include <set>

#include <unordered_map>
#include <unordered_set>
#include <iterator> 
#include <limits>
#include <ctime>

#include "mersenne_twister.hpp"
using namespace std;

const long MAXRAND=2147483647;
const double INF= numeric_limits<double>::infinity(); 
MTRand RanGen( (unsigned)time(NULL)); //Generate an instance from object MTRand, initiate it with the current time tt.


//Define the Biofilm class
template<class T>
class Biofilm{    
    
protected:
    int hei_bio;
    int wid_bio;
    int hei_cell;
    int wid_cell;
    T val_on;
    T val_off;
    vector< vector<T> > biofilm;
    
public:
    Biofilm(int hei_bio,int wid_bio, T val_on, T val_off, int hei_cell=2,int wid_cell=1): hei_bio(hei_bio), wid_bio(wid_bio), val_on(val_on), val_off(val_off), hei_cell(hei_cell), wid_cell(wid_cell){}
    //Initiate the biofilm 2D vector with given phi, rho and mobility.
    bool generateBiofilm(double phi,double rho,double mobility);
    //Visualize the spatial correlation strength
    double getSpaCorY();    
    //measure the covariance between cells along y axis
    vector<double> measureCoVarY( );  //return by value or return by reference.
    //measure the covariance between cells along both x axis
    vector<double> measureCoVarX( );
    //print the biofilm to a 2D matrix in a file
    bool print( FILE* fpt);
    //introduce mobility: every cell has the same possiblity to exchange with its neighbor.
    vector<T>& exchangeRow(vector<T>& row, int wid_bio, double mobility);
    //mother row gives birth to daughter row with correlation integrated, but not mobility.
    bool giveBirth(vector<T>& mother_row, vector<T>& daughter_row, double ponon, double ponoff);
    };

template<class T>
bool Biofilm<T>::generateBiofilm(double phi,double rho,double mobility){
    double pon=phi;
    double poff=1-pon;
	double ponon=rho+ (1-rho)* pon;   //p(d:on|m:on): rho=0, no spatial correlation. rho=1, perfect lineage inheritance
	double poffon=1-ponon;	
	double ponoff=(1- rho)* pon;        //p(d:on|m:off)
    double poffoff=1-ponoff;
    
    vector<T> appetizer( this->wid_bio );
    for(int j=0; j< this->wid_bio; ++j){
        if( RanGen.randExc()< phi )
            appetizer[j]= val_on;
        else
            appetizer[j]= val_off;   
    }
    
    for(int i=1; i<  30; ++i){ //Initiate appetizer for 30 times so that it reach a steady state.
        vector<T> tem( this->wid_bio );
        for(int j=0; j<this->wid_bio; ++j){
            giveBirth( appetizer, tem, ponon, ponoff);
            exchangeRow(tem, this->wid_bio, mobility);           
            for(int j=0; j<this->wid_bio; ++j)
                appetizer[j]= tem[j];
        }
    }
    
    biofilm.push_back( vector<T>(this->wid_bio) );
    for(int j=0; j<this->wid_bio; ++j)
        biofilm[0][j]= appetizer[j];
            
    for(int i=1; i< this-> hei_bio; ++i){
        biofilm.push_back( vector<T>(this->wid_bio) );
        giveBirth( biofilm[i-1], biofilm[i], ponon, ponoff);
        exchangeRow( biofilm[i], this->wid_bio, mobility);             
        }
    return true;
    }
    

    

template<class T>
double Biofilm<T>::getSpaCorY(){
    int n_onon=0, n_offon=0, n_onoff=0, n_offoff=0;  //n_daughter,mother
    for(int i=0; i< this->hei_bio -1; ++i){
        for(int j=0; j< this->wid_bio; ++j){
            if(biofilm[i][j]== val_on and biofilm[i+1][j]== val_on )
                n_onon+=1;
            if(biofilm[i][j]== val_on and biofilm[i+1][j]== val_off )
                n_offon+=1;
            if(biofilm[i][j]== val_off and biofilm[i+1][j]== val_on )
                n_onoff+=1;
            if(biofilm[i][j]== val_off and biofilm[i+1][j]== val_off )
                n_offoff+=1;           
        }
        
    double p_onon= n_onon/(n_onon+ n_offon);
    double p_onoff=n_onoff/(n_onoff+ n_offoff);
    double rho=p_onon- p_onoff;
    return rho;
    }
}

template<class T>    
vector<double> Biofilm<T>::measureCoVarY(  ){
    vector<double> covariance_y( this->hei_bio, 0 );
    vector<double> mean_y(this->hei_bio, 0);
    for (int j=0; j< this->wid_bio; ++j){
        for( int d=0; d< this->hei_bio; ++d){
            mean_y[d]+= this->biofilm[d][j];             
        }        
    }
    
    for( int d=0; d< this->hei_bio; ++d)
            mean_y[d]/= this->wid_bio;
            
    for (int j=0; j< this->wid_bio; ++j){
        for( int d=0; d< this->hei_bio; ++d){
            covariance_y[d]+= (this->biofilm[0][j] -mean_y[0] )*( this->biofilm[d][j]- mean_y[d] );             
        }        
    }
    
    for( int d=0; d< this->hei_bio; ++d)
            covariance_y[d]/= this->wid_bio;
    
    return covariance_y;          
}

template<class T>    
vector<double> Biofilm<T>::measureCoVarX(  ){
    vector<double> covariance_x( this->wid_bio, 0 );
    vector<double> mean_x(this->wid_bio, 0);
    for (int i=0; i< this->hei_bio; ++i){
        for( int d=0; d< this->wid_bio; ++d){
            mean_x[d]+= this->biofilm[i][d];             
        }        
    }
    
    for( int d=0; d< this->wid_bio; ++d)
            mean_x[d]/= this->hei_bio;
            
    for (int i=0; i< this->hei_bio; ++i){
        for( int d=0; d< this->wid_bio; ++d){
            covariance_x[d]+= (this->biofilm[i][0] -mean_x[0] )*( this->biofilm[i][d]- mean_x[d] );             
        }        
    }
    
    for( int d=0; d< this->hei_bio; ++d)
            covariance_x[d]/= this->hei_bio;
    
    return covariance_x;          
}

template<class T>
bool Biofilm<T>::print(FILE* fpt){
    if ( fpt.is_open()== false )
        return false

    for(int i=0;i< this->hei_bio; ++i){        
        for(int j=0; j< this->wid_bio; ++j )
            fpt << this->biofilm[i][j] <<"\t";
        fpt <<"\n";
    }
    
    return true;
}
    

template<class  T>
vector<T>& Biofilm<T>::exchangeRow(vector<T>& row, int wid_bio, double mobility){
    for(int j=0; j< wid_bio; ++j){
        if (RanGen.randExc()< 0.5 ){
            if(RanGen.randExc()< mobility ){
                int left= (j==0)? (wid_bio-1): j-1;
                T tem_swap= row[j];
                row[j]= row[left];
                row[left]= tem_swap;    
            }    
        }
        else{
            if(RanGen.randExc()< mobility ){
                int right= (j== (wid_bio-1) )? (0): j+1;
                T tem_swap= row[j];
                row[j]= row[right];
                row[right]= tem_swap;    
            }             
        }
    }
    return row;
}

template<class  T>
bool Biofilm<T>::giveBirth(vector<T>& mother_row, vector<T>& daughter_row, double ponon, double ponoff){
    for(int j=0; j<this->wid_bio; ++j){
        if( mother_row[j]==val_on){
            if( RanGen.randExc()< ponon )
                daughter_row[j]= val_on;
            else
                daughter_row[j]= val_off;
            }
        else{
            if( RanGen.randExc()< ponoff )
                daughter_row[j]= val_on;
            else
                daughter_row[j]= val_off;
            }
        }
    return true;    
    }
    
int main(){
    int hei_bio= 35, wid_bio= 230, repeats=100;
    double phi=0.43, rho_dynamic=0.38, rho_spatial=0.17, mobilities[]= {0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49}; 
    
    double rho_sparials[ sizeof(mobilities)/sizeof(mobilities[0] ) ];
    
    
    Biofilm<int> biofilm(hei_bio, wid_bio, 1, 0, 2, 1);
    
    
    biofilm.generateBiofilm(1, 0, phi, rho_dynamic, 0);
    
    FILE* fpt;
    fpt=fopen("./test/biofilm.dat","w");
    biogilm.print();
    double tem= biofilm.getSpaCorY();
    
    cout<<tem<<"The calculated rho" <<endl;
    return 0;   
    }






















