#include "DataSet.hh"

  DataSet::DataOccurrenceIterator::DataOccurrenceIterator(DataSet* parent)
  {
    parent_ = parent;
    data_index_ = 0;
    data_assignment_index_ = 0;
    current_data_end_ = parent_->data_counts_[ 0 ];
  }
  
  DataSet::DataOccurrenceIterator::~DataOccurrenceIterator() { }
  
  bool DataSet::DataOccurrenceIterator::Done()
  {
    return data_index_ >= parent_->data_counts_.size();
  } 
  
  void DataSet::DataOccurrenceIterator::Next()
  {
    ++( data_assignment_index_ );
    if ( data_assignment_index_ >= current_data_end_ )
    {
      ++data_index_;
      current_data_end_ += parent_->data_counts_[ data_index_ ];
    }      
  }
  
  int DataSet::DataOccurrenceIterator::Assignment(){
    return parent_->data_assignments_[ data_assignment_index_ ];
  }
  
  void DataSet::DataOccurrenceIterator::Reassignment( int new_assignment )
  {
    if ( 0 <= new_assignment && new_assignment <= parent_->assignment_distribution_.size() )
    {
      parent_->assignment_distribution_[ Assignment() ] -= 1;
      parent_->assignment_distribution_[ new_assignment ] += 1;
      parent_->data_assignments_[ data_assignment_index_ ] = new_assignment;
    }      
  }
  
  int DataSet::DataOccurrenceIterator::Data()
  {
    return data_index_;
  }
  
  //Condenses the vector of data_assignments_, and makes all data_counts_ into 1  
  //tested in aaa-test-iteration.cc
  void DataSet::CondenseInitialData()
  {
    long current_data_end = 0;
    for ( long i = 0; i < data_counts_.size(); ++i )
    {
      //cout << "\n";
      //cout << current_data_end << "\t";
      if ( data_counts_[ i ] > 1 )
      {
        data_assignments_.erase( data_assignments_.begin() + current_data_end, 
                                 data_assignments_.begin() + current_data_end - 1 + data_counts_[i] );
      } 
      current_data_end += 1;
      data_counts_[ i ] = 1;
    } 
    data_number_ = data_counts_.size();
  }
  
/*  DataSet::DataSet( const vector<int> &assignments, const vector<int> &start_indecies, int num_assignments )
  {
    data_number_ = assignments.size();
    assignments_number_ = num_assignments;
    data_assignments_ = assignments;
    data_counts_ = start_indecies;
    assignment_distribution_.resize( num_assignments );
    CountAssignmentDistribution();
    assignment_distribution_assembly_.resize( num_assignments );
  } 
*/  
  DataSet::DataSet( const char* filename )
  {
    using std::cerr;
    using std::ifstream;
    using std::ofstream;
    using std::istringstream;
    using std::string;
    ifstream fin( filename );
  
    if (!fin)
    {
      cout << "Error opening the file!...\n";
      exit(1);
    }
    else
    {    
      string line;    
      int assignment;
      int cocurrence;
      long assignment_length;
      string assignment_name;
      
      //the 1st line should be total hits number
      fin >> assignments_number_;
      getline( fin, line );
      
      //the 2nd line should be reference-names
      getline( fin, line );
      istringstream assignments_names( line );
      while ( assignments_names >> assignment_name )
      {
        assignment_name_.push_back( assignment_name );
      }
      
      //the 3rd line should be reference-lengths
      getline( fin, line );
      istringstream assignments_lengths( line );
      while ( assignments_lengths >> assignment_length )
      {
        assignment_length_.push_back( assignment_length );
      }
      
      //the 4th line should be the reference-id of those hits
      //the reference-id will be set in bam_process.py
      getline( fin, line );
      istringstream assignments_istream( line ); 
      while ( assignments_istream >> assignment )
      {
        data_assignments_.push_back( assignment );
      }
      
      //the 5th line should be the multi-hit counts of each read
      getline( fin, line );
      istringstream cocurrences_istream( line ); 
      while ( cocurrences_istream >> cocurrence )
      {
        data_counts_.push_back( cocurrence );
      }
                  
      data_number_ = data_counts_.size();
      assignment_distribution_.resize( assignments_number_ );
      assignment_distribution_initial_.resize( assignments_number_ );
      assignment_distribution_assembly_.resize( assignments_number_ );
      CountAssignmentDistribution();
      //assignment_distribution_initial_ = assignment_distribution_;
      for ( int i = 0; i < assignment_distribution_initial_.size(); ++i )
      {
        assignment_distribution_initial_[i] = assignment_distribution_[i];
      }
    //*dataset = DataSet( assignments, cocurrences, read_number );
    }  
    
    //calculate prior
    double assignment_init_sum = 0.0;
    long assignment_counts = 0;
    for ( int i = 0; i < assignment_distribution_.size(); ++i )
    {
      assignment_counts += assignment_distribution_initial_[ i ];
      assignment_init_sum += (double) assignment_distribution_initial_[i] / assignment_length_[i];
      //
      //cout << assignment_distribution_[ i ] << "\t" << assignment_distribution_initial_[ i ] << "\t" << assignment_counts << "\n";
    }
    //
    //cout << "\nthe assignments count(contains multi-alignments) is:" << assignment_counts << "\n\n";
    assignment_alpha_.resize( assignment_distribution_.size() );
    for ( int i = 0; i < assignment_distribution_.size(); ++i )
    {
      double temp = (double)assignment_distribution_initial_[i] / assignment_length_[i] / assignment_init_sum * assignment_counts;
      //cout << "\t" << temp;
      assignment_alpha_.push_back( (int)temp );
    }
    //cout << "\n\n";
  }
  
  void DataSet::CountAssignmentDistribution()
  {
    for (int i = 0; i < assignment_distribution_.size(); ++i) 
    {
      assignment_distribution_[i] = 0;
    }
    for (DataOccurrenceIterator iter(this); !iter.Done(); iter.Next()) 
    {
      assignment_distribution_[iter.Assignment()] += 1;
    }
  }

  
  //sum the assignments distribution after burn-in period, in order to calculate means
  void DataSet::AssembleAssignmentDistribution()
  {
    for ( int i = 0; i < assignment_distribution_.size(); ++i )
    {
      assignment_distribution_assembly_[i] += assignment_distribution_[i];
    }
  }
  
  void DataSet::CalculateAbundance( const char* filename )
  { 
    //first calculate mixture coefficient thetas   
    vector<double> assignment_theta;
    assignment_theta.resize( assignments_number_ );
    double assignment_count_sum;
    for ( int i = 0; i < assignment_distribution_.size(); ++i )
    {
      assignment_count_sum += assignment_distribution_[i];
      //assignment_count_sum += assignment_alpha_[i];
    }
    for ( int i = 0; i < assignment_distribution_.size(); ++i )
    {
      assignment_theta[i] = ( assignment_distribution_[i] /*+ assignment_alpha_[i]*/ ) / assignment_count_sum;
    }
   
    //then calculate relative abundance according to thetas and ref_lengths
    double assignment_relative_abundance1;
    double assignment_relative_abundance2;
    double assignment_relative_abundance3;
    double assignment_theta_sum;
    double assignment_align_sum;
    double assignment_init_sum;
    ofstream OutputFile( filename );
    for ( int i = 0; i < assignment_distribution_.size(); ++i )
    {
      assignment_theta_sum += assignment_theta[i] / assignment_length_[i];
    }
    for ( int i = 0; i < assignment_distribution_.size(); ++i )
    {
      assignment_align_sum += (double) assignment_distribution_[i] / assignment_length_[i];
      assignment_init_sum += (double) assignment_distribution_initial_[i] / assignment_length_[i];
      //std::cout<< (double)assignment_distribution_[i] / assignment_length_[i] << "\n";
    }
/*    
    ofstream OutputFile( filename );
    for ( int i = 0; i < assignment_name_.size(); ++i )
    {
      OutputFile << assignment_name_[ i ];
    }
    
    for ( int i = 0; i < assignment_distribution_.size(); ++i )
    {
      assignment_relative_abundance = assignment_theta[i] / assignment_length_[i] / assignment_theta_sum;
      std::cout<< assignment_relative_abundance << "\n";
      OutputFile << assignment_relative_abundance << "\t";
    }
    
    //for comparation
    //just assignment distribution from alignment
    OutputFile << "\n";
    for ( int i = 0; i < assignment_distribution_.size(); ++i )
    {
      OutputFile << assignment_distribution_[i] << "\t";
    }
    
    //distribution ajusted by ref_lengths
    OutputFile << "\n";
    double assignment_align_sum;
    for ( int i = 0; i < assignment_distribution_.size(); ++i )
    {
      assignment_align_sum += (double) assignment_distribution_[i] / assignment_length_[i];
      std::cout<< (double)assignment_distribution_[i] / assignment_length_[i] << "\n";
    }
    std::cout << "\n";
    for ( int i = 0; i < assignment_distribution_.size(); ++i )
    {
      assignment_relative_abundance = (double)assignment_distribution_[i] / assignment_length_[i] / assignment_align_sum;
      std::cout<< assignment_relative_abundance << "\n";
      OutputFile << assignment_relative_abundance << "\t";
    }
    //for comparation
    
*/
    //
    //OutputFile << "\n\n";
    OutputFile << "Speicies\testimated_abundance\tassigned_counts\tnormalized_assigned_counts\tmap_counts\tnormalized_map_counts\n";
    
    for ( int i = 0; i < assignment_distribution_.size(); ++i )
    {
      assignment_relative_abundance1 = assignment_theta[i] / assignment_length_[i] / assignment_theta_sum;
      std::cout<< assignment_relative_abundance1 << "\n";
      OutputFile << assignment_name_[i] << "\t";
      OutputFile << assignment_relative_abundance1 << "\t";
      OutputFile << assignment_distribution_[i] << "\t";
      assignment_relative_abundance2 = (double)assignment_distribution_[i] / assignment_length_[i] / assignment_align_sum;
      OutputFile << assignment_relative_abundance2 << "\t";
      OutputFile << assignment_distribution_initial_[i] << "\t";
      assignment_relative_abundance3 = (double)assignment_distribution_initial_[i] / assignment_length_[i] / assignment_align_sum;
      OutputFile << assignment_relative_abundance3 << "\n";
      
    }
    
    OutputFile.close();
  }

