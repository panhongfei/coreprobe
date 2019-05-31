#include <vector>
#include "StatSet.hh"

  StatSet::StatSet( int num_assignments, long num_data )
  {
    memory_alloc_.resize( num_assignments * (num_data + 1), 0 );
    
    assignment_distributions_.resize( num_data );    
    global_distribution_.Reset( &memory_alloc_[0] + num_assignments * num_data, num_assignments );
    for (int i = 0; i < num_data; ++i) 
    {
      assignment_distributions_[i] = AssignmentDistributionPointer(&memory_alloc_[0] + num_assignments * i, num_assignments); 
    }
  }
  
  void StatSet::InitiateStat( DataSet& dataset )
  {
    for (DataSet::DataOccurrenceIterator iter(&dataset); !iter.Done(); iter.Next()) 
    {
      IncrementAssignment(iter.Data(), iter.Assignment(), 1);
    }
  }
  
  const StatSet::AssignmentDistributionPointer& StatSet::GetAssignmentDistribution( int data ) const 
  {
    return assignment_distributions_[ data ];
  }
  
  const StatSet::AssignmentDistributionPointer& StatSet::GetGlobalDistribution() const 
  {
    return global_distribution_;
  }
  
  void StatSet::IncrementAssignment( int data, int assignment, int count ) 
  {
    if ( num_assignments() >= assignment ) {
      if ( num_data() >= data ) {
        assignment_distributions_[data][assignment] += count;
        global_distribution_[assignment] += count;
      }
    }
  }
  
  void StatSet::ReStat( int data, int old_assignment, int new_assignment ) 
  {
    IncrementAssignment(data, old_assignment, -1);
    IncrementAssignment(data, new_assignment, 1);
  }
  
  
  
