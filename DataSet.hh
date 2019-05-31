#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
using namespace std;

//this class stores a bag of data 
//in this application it would be metagenomic reads information
class DataSet{
    
  public:
 
  //an iterator over all of the data occurence
  //tested by aaa-test-iteration.cc
  class DataOccurrenceIterator 
  {
    public:
    //initializes the DataOccurrenceIterator for a DataSet
    explicit DataOccurrenceIterator( DataSet* parent_ );
    ~DataOccurrenceIterator();
    
    //Returns true if iteration is done
    bool Done();
    
    //Advances to the next data occurrence
    void Next();
    
    //Returns the assignment of the current data occurrence
    int Assignment();
    
    //Changes the assignment of the current data occurrence
    void Reassignment( int new_assignment );
    
    //Returns the current data occurrence
    int Data();
    
    private:
    
    DataSet* parent_;
    long data_index_;
    long data_assignment_index_;   
    long current_data_end_;
  };
  
  friend class DataOccurrenceIterator;
  
  DataSet() { }
    
  //We use SAM/BAM file with all hits for initiation, which will contribute to constructing
  //a probability matrix in StatSet.cc. And then the initial dataset will be condensed by 
  //CondenseInitialData() to represent the original real reads set.
  DataSet( const char* filename );
  
  ~DataSet() { }  
  
  //Sums the counts of each assignment after burn-in period
  void AssembleAssignmentDistribution();
  
  //Calculates the abundance according to Charlie's formula and output to a file
  void CalculateAbundance( const char* filename );
  
  //In our algorithm, we simply assume real metagenomic data as non-redundant reads set,
  //because illumina reads are rarely the same in real metagenomic data (150bp reads).
  //During Initialization, we use SAM file with multi-hits, so we need condense the read dataset 
  void CondenseInitialData();
  
  const vector<int>& assignment_distributions() const 
  {
    return assignment_distribution_;  
  }
  
  const vector<int>& data_assignments() const 
  {
    return data_assignments_;
  }
  
  const vector<int>& assignment_distribution_assembly() const
  {
    return assignment_distribution_assembly_;
  }
  
  const long get_data_number() const
  {
    return data_number_;
  }
  
  const int get_assignments_number() const
  {
    return assignments_number_;
  }
  
  const vector<int> get_assignment_alpha() const
  {
    cout << "the priors are: " << endl;
    for ( int i = 0; i < assignment_alpha_.size(); ++i )
    {
      cout<<assignment_alpha_[i] << "\t";
    }
    cout << "\n\n";
    
    return assignment_alpha_;
  }
  
  protected:
  //Records each data occurrence's assignment
  vector<int> data_assignments_;
  //Records how many times each data occurs 
  vector<int> data_counts_;
  //Records each assignment's occurrence number, 
  //which may change during each Gibbs run
  vector<int> assignment_distribution_;
  //Records each assignment's pseudo counts
  vector<int> assignment_alpha_;
  //Records each assignment's character length
  vector<long> assignment_length_;
  vector<string> assignment_name_;
  //Adds distributions after burn-in period
  vector<int> assignment_distribution_assembly_;
  //Records distribution in initial SAM file
  vector<int> assignment_distribution_initial_;
  
  long data_number_;
  int assignments_number_;
  
  //Counts the assignments each during Gibbs Sampling
  void CountAssignmentDistribution();
  
};
