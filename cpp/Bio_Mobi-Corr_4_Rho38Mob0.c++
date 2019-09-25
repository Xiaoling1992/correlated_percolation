//_4. New approach to gorw the biofilm to introduce in x axis: giveBirthOdd, exchangeOdd.
//_3. Different definition of correlation function in percolation theory: at a distance d, the probability that cell(d) is the same with cell(0). same: 1, different, -1;
//_2. when exchange rate is 0.4, rho_spatial=0.17. Now set exchange rate to be 0.4, visualize the correlation length.
//2. Estimate mobility, then measure the correlation function in both x and y axis.
//1. Introduce mobility (cells can exchange with neighbors in x axis) to reduce rho_y to 0.17


#include <iostream>
#include <fstream>
#include <cstdio>
#include <utility> 

#include <vector>
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

template<class T>
bool printVector(vector<T> output, string output_path){
    
    ofstream fpt;
    fpt.open(output_path);
    for(int i=0; i< output.size(); ++i )
        fpt<< output[i] <<"\t";
    
    fpt<<endl;
    fpt.close();    
    return true;
}

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
    
    vector<double> measureCorFunY( double phi );
    vector<double> measureCorFunX(  double phi);
    //print the biofilm to a 2D matrix in a file
    bool print( string output_path);
    //introduce mobility: every cell has the same possiblity to exchange with its neighbor.
    bool exchangeRow(vector<T>& row, int wid_bio, double mobility);
    bool exchangeRowOdd(vector<T>& above, vector<T>& below, int wid_bio, double mobility);
    bool exchangeRowEven(vector<T>& row, int wid_bio, double mobility);
    //mother row gives birth to daughter row with correlation integrated, but not mobility.
    bool giveBirth(vector<T>& mother_row, vector<T>& daughter_row, double ponon, double ponoff);
    bool giveBirthOdd(vector<T>& mother_row, vector<T>& daughter_row, double ponon, double ponoff);
    bool giveBirthEven(vector<T>& mother_row, vector<T>& daughter_row, double ponon, double ponoff);
    };

template<class T>
bool Biofilm<T>::generateBiofilm(double phi,double rho,double mobility){
    double pon=phi;
    double poff=1-pon;
	double ponon=rho+ (1-rho)* pon;   //p(d:on|m:on): rho=0, no spatial correlation. rho=1, perfect lineage inheritance
	double poffon=1-ponon;	
	double ponoff=(1- rho)* pon;        //p(d:on|m:off)
    double poffoff=1-ponoff;
    
    vector<T> appetizer( this->wid_bio ); //Initiate appetizer
    for(int j=0; j< this->wid_bio; ++j){
        if( RanGen.randExc()< phi )
            appetizer[j]= val_on;
        else
            appetizer[j]= val_off;   
    }
    
    vector<T> tem( this->wid_bio );
    for(int i=1; i<  30; ++i){ //Initiate appetizer for 30 times so that it reach a steady state.
        giveBirthOdd( appetizer, tem, ponon, ponoff);
        exchangeRowOdd(appetizer, tem, this->wid_bio, mobility);    
        
        giveBirthEven( appetizer, tem, ponon, ponoff);
        exchangeRowEven( tem, this->wid_bio, mobility);     
           
        for(int j=0; j<this->wid_bio; ++j)
            appetizer[j]= tem[j];
    }

    biofilm.push_back( vector<T>(this->wid_bio) );
    for(int j=0; j<this->wid_bio; ++j)  //use deep copy
        biofilm[0][j]= appetizer[j];
            
    for(int i=1; i< this-> hei_bio; ++i){
        biofilm.push_back( vector<T>(this->wid_bio) );
        
        giveBirthOdd( biofilm[i-1], biofilm[i], ponon, ponoff);
        exchangeRowOdd(biofilm[i-1], biofilm[i], this->wid_bio, mobility);    
        
        giveBirthEven( biofilm[i-1], biofilm[i], ponon, ponoff);
        exchangeRowEven( biofilm[i], this->wid_bio, mobility);     
        }
            
    return true;
    }
    

    

template<class T>
double Biofilm<T>::getSpaCorY(){
    double n_onon=0, n_offon=0, n_onoff=0, n_offoff=0;  //n_daughter,mother
    
    for(int i=0; i< this->hei_bio -1; ++i){
        for(int j=0; j< this->wid_bio; ++j){
            if(biofilm[i][j]== val_on && biofilm[i+1][j]== val_on )
                n_onon+=1;
            else if(biofilm[i][j]== val_on && biofilm[i+1][j]== val_off )
                n_offon+=1;
            else if(biofilm[i][j]== val_off && biofilm[i+1][j]== val_on )
                n_onoff+=1;
            else if(biofilm[i][j]== val_off && biofilm[i+1][j]== val_off )
                n_offoff+=1;           
        }
    }
        
    //cout<<n_onon<<" offon "<<n_offon <<" onoff "<<n_onoff  <<" offoff "<<n_offoff <<endl; 
    double p_onon= n_onon/(n_onon+ n_offon);
    double p_onoff=n_onoff/(n_onoff+ n_offoff);
    double rho=p_onon- p_onoff;
    return rho;
  
}

template<class T>    
vector<double> Biofilm<T>::measureCoVarY(  ){ //calculate the covariance between the row 0 and row i. Used the sample mean of each row instead of phi.
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
vector<double> Biofilm<T>::measureCoVarX(  ){  //calculate the covariance between the column 0 and column d. Used the sample mean of each column instead of phi.
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
vector<double> Biofilm<T>::measureCorFunY( double phi ){ //calculate the auto-covariance function between row 0 and row d.
    vector<double> covariance_y( this->hei_bio, 0 );
    
    for( int d=0; d< this->hei_bio; ++d){        
        for (int j=0; j< this->wid_bio; ++j)        
            covariance_y[d]+= (this->biofilm[0][j] )*( this->biofilm[d][j] );             
        covariance_y[d]= covariance_y[d]/this->wid_bio- phi*phi;      
    }
    return covariance_y;          
}

template<class T>    
vector<double> Biofilm<T>::measureCorFunX( double phi ){ //calculate the auto-covariance function between column 0 and column d.
    vector<double> covariance_x( this->wid_bio, 0 );
    
    for( int d=0; d< this->wid_bio; ++d){
        
        for (int i=0; i< this->hei_bio; ++i)    //y= 1, 3...2h-1       
            covariance_x[d]+= (this->biofilm[i][0] )*( this->biofilm[i][d]);   
        for (int i=1; i< this->hei_bio; ++i)    //y=2, 4, 6, 2h-2   
            covariance_x[d]+= (this->biofilm[i-1][0] )*( (d%2== 0)? this->biofilm[i-1][d]: this->biofilm[i][d] );
            
        covariance_x[d]= covariance_x[d]/(2* this->hei_bio-1)- phi*phi;
    }
    
    return covariance_x;          
}

template<class T>
bool Biofilm<T>::print( string output_path){
    ofstream fpt;
    fpt.open(output_path);
    for(int i=0;i< biofilm.size(); ++i){        
        for(int j=0; j< biofilm[i].size() ; ++j )
            fpt<< this->biofilm[i][j] <<"\t";
        fpt<<endl;
    }
    
    fpt.close();
    
    return true;
}
    

template<class  T>
bool Biofilm<T>::exchangeRow(vector<T>& row, int wid_bio, double mobility){
    for(int j=0; j< this->wid_bio; ++j){
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
    return true;
}

template<class  T>
bool Biofilm<T>::exchangeRowOdd(vector<T>& above, vector<T>& below, int wid_bio, double mobility){
    for(int j=0; j< this->wid_bio; ++j){
        if(j%2!= 1) //not a odd column
            continue;
        if (RanGen.randExc()< 0.5 ){
            if(RanGen.randExc()< mobility ){
                int left= (j==0)? (wid_bio-1): j-1;
                T tem_swap= below[j];
                below[j]= above[left];
                above[left]= tem_swap;    
            }    
        }
        else{
            if(RanGen.randExc()< mobility ){
                int right= (j== (wid_bio-1) )? (0): j+1;
                T tem_swap= below[j];
                below[j]= above[right];
                above[right]= tem_swap;    
            }             
        }
    }
    return true;
}

template<class  T>
bool Biofilm<T>::exchangeRowEven(vector<T>& row, int wid_bio, double mobility){
    for(int j=0; j< this->wid_bio; ++j){
        if(j%2!= 0) //not a even column
            continue;
            
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
    return true;
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
    
template<class T>    
bool Biofilm<T>::giveBirthOdd(vector<T>& mother_row, vector<T>& daughter_row, double ponon, double ponoff){
    for(int j=0; j<this->wid_bio; ++j){
        if(j%2 != 1) //If it is not an odd column
            continue;
            
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

template<class T>    
bool Biofilm<T>::giveBirthEven(vector<T>& mother_row, vector<T>& daughter_row, double ponon, double ponoff){
    for(int j=0; j<this->wid_bio; ++j){
        if(j%2 != 0) //If it is not an even column
            continue;
            
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
    int hei_bio= 35, wid_bio= 230, repeats=10000;
    double phi=0.43, rho_dynamic=0.38, rho_spatial=0.17, mobilities[]= { 0 }; //mobility=0.36, the rho_spatial is 0.1712
    vector<double> rho_spatials(sizeof(mobilities)/sizeof(mobilities[0] ), 0 );
    vector<double> covariance_y(hei_bio, 0);
    vector<double> covariance_x(wid_bio, 0);
    //for( int i=0; i< sizeof(mobilities)/sizeof(mobilities[0] ); ++i){
        
        for(int irepeat=0; irepeat<repeats; ++irepeat){
            vector<double> covariance_y_tem;
            vector<double> covariance_x_tem;
            Biofilm<int> biofilm(hei_bio, wid_bio, 1, 0, 2, 1); //
            
            biofilm.generateBiofilm(phi, rho_dynamic, mobilities[0] );
            covariance_y_tem= biofilm.measureCorFunY( phi);
            covariance_x_tem= biofilm.measureCorFunX( phi);
            
            for(int i=0; i<hei_bio; ++i)
                covariance_y[i]+= covariance_y_tem[i];
            for(int j=0; j<wid_bio; ++j)
                covariance_x[j]+= covariance_x_tem[j];
            
            if(irepeat==0)
                biofilm.print( "/home/xiaoling/correlated_percolation/data/Bio4_biofilm_Rho38Mob0.dat" );
        }
        for(int i=0; i<hei_bio; ++i)
            covariance_y[i]/= repeats;
        for(int j=0; j<wid_bio; ++j)
            covariance_x[j]/= repeats;
    
    
    printVector<double>(covariance_y, "/home/xiaoling/correlated_percolation/data/Bio4_CF_y_r10000_Rho38Mob0.dat");
    printVector<double>(covariance_x, "/home/xiaoling/correlated_percolation/data/Bio4_CF_x_r10000_Rho38Mob0.dat");
    //biofilm.print( "/home/xiaoling/correlated_percolation/test/biofilm1.dat" );
    
    //cout<<tem<<"The calculated rho" <<endl;
    return 0;   
    }






















