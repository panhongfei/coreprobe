#include <vector>
#include <cstring>
#include "DataSet.hh"

//This class stores useful statistics of data 
//In this application the data would be metagenomic reads
class StatSet {
  
  public:
  //this class pointers to elements in probability matrix 
  class AssignmentDistributionPointer {
    
    public:
    
    AssignmentDistributionPointer(): distribution_(NULL), size_(0) {}
    
    AssignmentDistributionPointer(int* distribution, int size)
                        : distribution_(distribution), size_(size) {}
                        
    void Reset(int* distribution, int size){
      distribution_=distribution;
      size_=size;
    }
    
    int size() const { return size_; }
    inline int& operator[](int index) const { return distribution_[index]; }
    void clear() { memset(distribution_, 0, sizeof(*distribution_) * size_); }
    
    private:
    //Pointer to assignment distribution in memory_alloc_
    int* distribution_;
    //
    int size_; //may be static???
  };
  
  friend class AssignmentDistributionPointer;
  
  StatSet( int num_assignments, long num_data );
  
  ~StatSet() {}
  
  void InitiateStat( DataSet &dataset );
  
  const AssignmentDistributionPointer& GetAssignmentDistribution( int data ) const;
  
  const AssignmentDistributionPointer& GetGlobalDistribution() const;
  
  void IncrementAssignment( int data, int assignment, int count );
  
  void ReStat( int data, int old_assignment, int new_assignment );
  
  int num_assignments() const { return global_distribution_.size(); }
  
  int num_data() const { return assignment_distributions_.size(); }
  
    
  protected:
  //Memory allocated for assignment counts for each data
  vector<int> memory_alloc_;
  
  private:
  //Vector storing pointers to memory location of assignment counts, and its size
  vector<AssignmentDistributionPointer> assignment_distributions_;
  //Pointer to memory location of global assignment counts
  AssignmentDistributionPointer global_distribution_; 
  
};
