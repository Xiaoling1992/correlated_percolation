// Bio6: If the neighbor is occupied but the neighbor's neighbor is free, the daughter can push the neighbor cell to the neighbor's neighboring location.
// 5.2 _bot-per. Answer the question: sigma->0, whey is there no correlation in x-axis. Assumption: sigma->0: most daughter cells are in the bottom locations.
// 5.2: The daughter cell could be any free neighboring cells of the mother cell.
//_5. Draw the growth time of a cell from a gaussian distribution
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
double assignPhi(double phi, T on, T off){
     if( RanGen.randExc()< phi )
            return on;
}

double gaussrand_K()
{
	static double V1, V2, S;
	static int phase = 0;
	double X;

	if(phase == 0) {
		do {
			double U1 = RanGen.rand();   //rand() [0, 1]
			double U2 = RanGen.rand();
            //cout<< "U1 U2 are" <<U1 <<U2<<endl;
			V1 = 2 * U1 - 1;
			V2 = 2 * U2 - 1;
			S = V1 * V1 + V2 * V2;
			} while(S >= 1 || S == 0);

		X = V1 * sqrt(-2 * log(S) / S);
	} else
		X = V2 * sqrt(-2 * log(S) / S);

	phase = 1 - phase;

	return X;
}

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

//Define a (location, division time) priority queue
class KDivision{		  
	private:
	  int k;  //k= i*wid_bio+ j;
	  double division_time;
	
	public:
	  KDivision(int k, double division_time): k(k), division_time(division_time){};
	  
	  int getK(){
		  return this->k;
	  }
	  
	  double getDivision_time(){
	     return this->division_time;
	  }
	  
	  bool setDivision_time(double division_time){
		  this->division_time= division_time;
		  return true;
	  }
	  
	  friend bool operator<(const KDivision& point1, const KDivision& point2);  //1. Define <, so set can make it sorted.
};
//overload <: so that the minimum KDivision will be on the top of the priority_queue .	
bool operator<(const KDivision& point1, const KDivision& point2){//use const to protect the input variable.
	return point1.division_time> point2.division_time;  //We can return point1.division_time> point2.distance to reverse the sorting. (descending sorting or minimum top priorityqueue)
}
//Defien the priority queue: pop (top), insert, search, update.
//map has the complete information. pq is used to pop the minimum.
//If the top element in pq is illegal (not in map or not consistent), pop and top the new minimum.
class MyPQ{
  private:
	unordered_map<int, double> loc_time;   //key could be int, string.
    priority_queue<KDivision> pq;  //k, division time, minimum top.
	
  public:    
    MyPQ(){}
    
    bool empty(){
        return loc_time.size()==0;
        }
        
    int getSize(){
        return loc_time.size();
        }
    
    int top(){  //If the top is not in the map or is not consistent with the lement in the map, pop and get a new top.
        KDivision minimum= pq.top();
        int k= minimum.getK();
        double div_time= minimum.getDivision_time();
        
        while (loc_time.find(k) == loc_time.end() || loc_time[k] != div_time){ //not in map or not consistent
                
            pq.pop();
            minimum= pq.top();
            k= minimum.getK();
            div_time= minimum.getDivision_time();            
            }
            
        
        return k;
    }
    
    bool pop(){  //If the top is not in the map or is not consistent with the lement in the map, the top element in pq isn't a legal record. Pop the next one.
        KDivision minimum= pq.top();
        int k= minimum.getK();
        double div_time= minimum.getDivision_time();
        
        while (loc_time.find(k) == loc_time.end() || loc_time[k] != div_time){
            pq.pop();
            minimum= pq.top();
            k= minimum.getK();
            div_time= minimum.getDivision_time();            
            }
        pq.pop();
        loc_time.erase( k);      //erase by key or by iterator.     
        
        for(auto& x : loc_time)
            cout<<"    k is "<<x.first <<" time is "<<x.second;
        cout<<" size of pq is "<< pq.size()<<endl;
        
        return true;
        }
        
    bool erase (int k){
        if (loc_time.find(k) != loc_time.end() )
            loc_time.erase(k);
        return true;        
        }
        
    bool push(int k, double time){  //If k is already in the map, we can't push the new (k, time) in.
        pq.push( KDivision( k, time) );
        loc_time.insert( pair<int, double>(k, time) );  //insert( { {key1, value1}, {key2, value2}  } )
        return true;
        }
        
    bool search(int k){
        return loc_time.find(k)!= loc_time.end();   
        }
    
    double getTime(int k){
        return loc_time[k];
        }
        
    bool updateTime(int k, double time){ //update the time in map, push the new (k, time) to the priority_queue
        loc_time[k]= time;
        pq.push( KDivision(k, time) );
        return true;        
        }
	
	};

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
    //Initiate the biofilm with growth time introduced.
    double generateBiofilm_GR(double phi,double rho, double mean, double sigma);
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
 


double getGrowthTime(double mean, double sigma){
    
    double growth_time;
    do{
        double gaussian= gaussrand_K();
        //cout<< "gaussian is " << gaussian <<endl;
        growth_time= gaussian*sigma+mean;
        //cout<< "growth time is, inside do while " << growth_time <<endl;
    } while( growth_time<= 0 );
    //cout<< "growth time is " << growth_time <<endl;
    return growth_time;
    
    } 
 
//Introduce the growth rate: marked by the time need to grow before division.
template<class T>
double Biofilm<T>::generateBiofilm_GR(double phi,double rho, double mean, double sigma){
    double pon=phi;
    double poff=1-pon;
	double ponon=rho+ (1-rho)* pon;   //p(d:on|m:on): rho=0, no spatial correlation. rho=1, perfect lineage inheritance
	double poffon=1-ponon;	
	double ponoff=(1- rho)* pon;        //p(d:on|m:off)
    double poffoff=1-ponoff;
    //int hei_bio;
    //int wid_bio;
    //int hei_cell;
    //int wid_cell;
    //T val_on;
    //T val_off;
    //vector< vector<T> > biofilm;
    
    int hei_more=100, hei_cut=5; //cut down the last 5 rows;
    int length_bio_tem= (hei_bio+ hei_more)* wid_bio;
    double num_daugh_bottom=0;
    vector< T > biofilm_tem( length_bio_tem, MAXRAND ); //MAXRAND means this space is free.
    priority_queue<KDivision> pq;
    
    //Initiate the first raw.
    //for(int j=0; j<this->wid_bio; ++j){
        //biofilm_tem[j]= assignPhi<T>(phi, val_on, val_off); //has phi to be on
        //pq.push( KDivision( j, getGrowthTime(mean, sigma) ) );
    //}
    
    for(int j=0; j< this->wid_bio; j=j+1){
        biofilm_tem[j]= assignPhi<T>(phi, val_on, val_off); //has phi to be on
        pq.push( KDivision( j, getGrowthTime(mean, sigma) ) );
    }
    
    //Pop the mother whose division time is the earliest. Decide whether it has free neighboring spaces. 
    int left, right, bottom, free_space, daughter;
    int neighbors[6];
    vector<int> avai_neighbors;
    while(pq.empty() == false){  //There are still cells which didn't divide.
        KDivision mother= pq.top();
        pq.pop();
        int k= mother.getK();  //location and division time of mother cell.
        double mother_div_time= mother.getDivision_time();
        //cout<<"location "<<k<< " mother's cell, disivion time is "<<mother_div_time <<endl;
        int i= k /wid_bio;  //row i
        int j= k %wid_bio;  //column j
                
        neighbors[0]= k- wid_bio; //top
        neighbors[1]= k+ wid_bio; //bottom
        
        neighbors[2]= (j!=0) ? (k- 1) : (k-1+ wid_bio);  
        neighbors[3]= (j!= wid_bio-1) ? (k+ 1) : (k+1- wid_bio);  //right, k+1
        
        if (j%2 ==1 ){
            neighbors[4]= neighbors[2]- wid_bio;
            neighbors[5]= neighbors[3]- wid_bio;
        }
        else{
            neighbors[4]= neighbors[2]+ wid_bio;
            neighbors[5]= neighbors[3]+ wid_bio;            
            }
            
        avai_neighbors.resize(0);
        for(int i=0; i< 6; ++i){
            if ( neighbors[i]>=0 && neighbors[i]< length_bio_tem && biofilm_tem[ neighbors[i] ]==MAXRAND ) //the location is legal and free.
                avai_neighbors.push_back(  neighbors[i] );        
            }
            
        free_space= avai_neighbors.size();
        if (free_space==0)
            continue;
            
        daughter= avai_neighbors[ RanGen.randInt( free_space-1 ) ]; //pick up a daughter: [0, #daughters -1]
        
        if (daughter== neighbors[1])// If the daughter cell is the bottom cell
            num_daugh_bottom+= 1;
        
        //if the mother has free neighboring space, determine a location for its daughter.           
        if(biofilm_tem[k]== val_on)
            biofilm_tem[ daughter ]= assignPhi<T>(ponon, val_on, val_off); //has phi to be on
        else
            biofilm_tem[ daughter ]= assignPhi<T>(ponoff, val_on, val_off); //has phi to be on    
        
        if(free_space>= 2)  //If free_space==1, then the mother has no more avaiable neighboring locations
            pq.push( KDivision( k, mother_div_time+ getGrowthTime(mean, sigma) ) );  
            
        pq.push( KDivision( daughter, mother_div_time+ getGrowthTime(mean, sigma) ) );              
        }
        
    
    //int num_cells=0;
    //for(int i=0; i< length_bio_tem; ++i)
    //    num_cells+= (biofilm_tem[i] != MAXRAND); 
   
       
    biofilm= vector< vector<T> >( hei_bio, vector<T>(wid_bio ) );
    for (int i=0; i< hei_bio; ++i){
        for (int j=0; j<wid_bio; ++j)
            biofilm[i][j]= biofilm_tem[ (i+ hei_more- hei_cut)* wid_bio+ j];
    }
            
    return num_daugh_bottom/static_cast<double>(length_bio_tem- wid_bio);
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
    //int hei_bio= 35, wid_bio= 230, repeats=100;  //1. size;
    //double phi=0.43, rho_dynamic=0.38, rho_spatial=0.17, sigmas[]= { 0.0, 0.01, 0.03, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0 }; 
    //vector<double> rho_spatials( sizeof(sigmas)/sizeof(sigmas[0] ), 0 );
    //vector<double> perc_daugh_bottom( sizeof(sigmas)/sizeof(sigmas[0] ), 0 );
    ////vector<double> covariance_y(hei_bio, 0);
    ////vector<double> covariance_x(wid_bio, 0);
        
     //for( int i=0; i< sizeof(sigmas)/sizeof(sigmas[0] ); ++i){
         //cout<<"sigma is "<< sigmas[i]<<endl;
         
         ////vector<double> covariance_y(hei_bio, 0);
         ////vector<double> covariance_x(wid_bio, 0);
        
        //for(int irepeat=0; irepeat<repeats; ++irepeat){
            ////vector<double> covariance_y_tem;
            ////vector<double> covariance_x_tem;
            
            //Biofilm<int> biofilm(hei_bio, wid_bio, 1, 0, 2, 1);  //hei_bio, wid_bio, on, off, hei_cell, wid_cell
            ////generateBiofilm_GR(double phi,double rho,double mobility, double mean, double sigma, double mobi_left, double mobi_right);
            //perc_daugh_bottom[i]+= biofilm.generateBiofilm_GR(phi, rho_dynamic, 1.0, sigmas[i] );
            
            ////if (num_cells < (100+ hei_bio)* wid_bio )
            ////    cout<< num_cells <<" cells out of "<< (hei_bio+100)*wid_bio <<endl;
            
            ////biofilm.print("/home/xiaoling/test/Bio5_test.dat");
            ////rho_spatials[i]+= biofilm.getSpaCorY();
            
            ////covariance_y_tem= biofilm.measureCorFunY( phi);
            ////covariance_x_tem= biofilm.measureCorFunX( phi);  //return by value
            
            ////for(int i=0; i<hei_bio; ++i)
                ////covariance_y[i]+= covariance_y_tem[i];
            ////for(int j=0; j<wid_bio; ++j)
                ////covariance_x[j]+= covariance_x_tem[j];
        //}
        //perc_daugh_bottom[i]/= repeats;
        ////rho_spatials[i]/= repeats;
        ////for(int i=0; i<hei_bio; ++i)
            ////covariance_y[i]/= repeats;
        ////for(int j=0; j<wid_bio; ++j)
            ////covariance_x[j]/= repeats;   
    
        ////printVector<double>(covariance_y, "/home/xiaoling/correlated_percolation/data/Bio5_CF_y_r3000_Rho38Sigma"+ to_string( sigmas[i] )  +".dat");
        ////printVector<double>(covariance_x, "/home/xiaoling/correlated_percolation/data/Bio5_CF_x_r3000_Rho38Sigma" +  to_string( sigmas[i] ) +  ".dat");
    //}    
    
    //printVector<double>(perc_daugh_bottom, "/home/xiaoling/correlated_percolation/data/Bio5.2_perdaubot_sigma.dat");
    
    ////cout<<tem<<"The calculated rho" <<endl;
    //return 0;   
    
    MyPQ pq;
    pq.push(1, 1.5);
    pq.push(3, 1.1);
    pq.push(2, 1.6);
    pq.push(2, 1.3);
    
    while(pq.empty()==false ){
        cout<< pq.top() <<" time is "<<pq.getTime( pq.top() ) <<endl;
        pq.pop();       
        }
        
    cout<<"size is "<<pq.getSize()<<endl;
    
    pq.push(1, 1.5);
    pq.push(3, 1.1);
    pq.push(2, 1.6);
    
    if (pq.search(1) ){
  
        pq.updateTime(1, 0.5);
    }
    if (pq.search(3) ){

        pq.updateTime(3, 0.9);
        }
    if (pq.search(2) ){

        pq.updateTime(2, 1.9);
    }
    while(pq.empty()==false ){
        cout<< "current size is" << pq.getSize() <<endl;
        cout<< pq.top() <<" time is "<<pq.getTime( pq.top() ) <<endl;
        pq.pop() ;      
        }
        
    cout<<"size is "<<pq.getSize()<<endl;
    
    cout<<"-8%10 "<< -8 %10<< " -18%10 "<< -18%10<<endl;
    cout<<" -8/10 "<< -8 /10<< " -18/10 "<< -18/10<<endl;
    return 0;
    }
