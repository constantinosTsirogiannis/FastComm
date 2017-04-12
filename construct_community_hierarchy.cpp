
#include<Rcpp.h>
#include<vector>
#include<set>
#include<string>
#include<iostream>
#include<sstream>
#include<fstream>
#include<algorithm>
#include<list>
#include<map>
#include<queue>
#include<cstdlib>
#include<cctype>

///////////////////////////////////////////////////////////////////////////
//////////////// Function objects on Multipartite Networks ////////////////
///////////////////////////////////////////////////////////////////////////

namespace MultipartiteNetworkFunctors{

class Edge_list_inner_product
{
 public:

  template<class NODE_TYPE>
  double operator()( NODE_TYPE& n1, NODE_TYPE& n2 ) const
  {
    typedef NODE_TYPE  Node;
    typedef typename   Node::Neighbour_element     Neighbour_element;
    typedef typename   Node::Neighbour_iterator    Neighbour_iterator;

    int s1 = n1.number_of_neighbours(),
        s2 = n2.number_of_neighbours();

    double product(0.0);

    Neighbour_iterator nit_1 = n1.neighbours_begin(),
                       nit_2 = n2.neighbours_begin();

    while( nit_1 !=n1.neighbours_end() && nit_2 != n2.neighbours_end() )
    {
      Neighbour_element nlm1 = *nit_1, nlm2 = *nit_2;

      if(nlm1 == nlm2)
      {
        product += (nlm1.weight())*(nlm2.weight());
        nit_1++;
        nit_2++;
      }
      else if( nlm1 < nlm2 )
        nit_1++;
      else
        nit_2++;
    } 

    // Check also if the two nodes are neighbours

    Neighbour_element slm(n1.group(),n1.index(),double(0.0),0);

    Neighbour_iterator cit = n2.find_neighbour(slm);

    if( cit != n2.neighbours_end() )
    {
      product += n1.compute_predicate_value_1()*(cit->weight());
      product += n2.compute_predicate_value_1()*(cit->weight());
    }
   
    return product;

  } // double operator()(  Node& n1,  Node& n2 ) const


}; // class Edge_list_inner_product

 
class Tanimoto_coefficient
{
 public:

  template<class NODE_TYPE>
  double operator()( NODE_TYPE& n1, NODE_TYPE& n2 ) const 
  {
    double inner_product = Edge_list_inner_product()(n1,n2);

    return inner_product/(   n1.compute_predicate_value_2() 
                           + n2.compute_predicate_value_2() -inner_product ); 

  } // double operator()( const Node& n1, const Node& n2 ) const

}; // Tanimoto_coefficient

class Minimum_weight_sum
{
 public:

  template<class NODE_TYPE>
  double operator()( NODE_TYPE& n1, NODE_TYPE& n2 ) const 
  { return std::min(n1.compute_predicate_value_1(), n2.compute_predicate_value_1());} 

}; // Minimum_weight_sum

} // MultipartiteNetworkFunctors


///////////////////////////////////////////////////////////////////////////
//////////////////// Multipartite Network Edge Traits /////////////////////
///////////////////////////////////////////////////////////////////////////


namespace MultipartiteNetwork{

template < class NELM > 
struct Is_greater_element
{
  typedef NELM Neigh_element;

  bool operator()( const Neigh_element& ne1, const Neigh_element& ne2 ) const
  {
    if(ne1.group() < ne2.group())
      return true;

    if(ne1.group() > ne2.group())
      return false;

    if(ne1.index() < ne2.index())
      return true;

    return false;
  }
}; // Is_greater_element

class Edge_traits_base
{
 public:

  class Neighbour_element_base
  {
   public:
 
    Neighbour_element_base(int group, int ind, int edge_ind ):
    _group(group),_index(ind),_edge_index(edge_ind){}

    int index(void) const
    { return _index; }

    int group(void) const
    { return _group; }

    int edge_index(void) const
    { return _edge_index; }

    bool operator==(const Neighbour_element_base &nlm) const
    {
      if(nlm.group() != this->group() || nlm.index() != this->index())
        return false;

      return true;
    }

    bool operator<(const Neighbour_element_base &nlm) const
    { return Is_greater_element<Neighbour_element_base>()(*this,nlm); }

    bool operator>(const Neighbour_element_base &nlm) const
    { return Is_greater_element<Neighbour_element_base>()(nlm,*this); }

   private:

    int _group, _index, _edge_index;

  }; // class Neighbour_element_base

  class Primary_edge_base
  {
    public: 

      Primary_edge_base(){}

      int index(void) const
      { return _index; }

      void set_index(int index)
      { _index = index; }

    private:

     int _index;

  }; // Primary_edge_base

}; // class Edge_traits_base

class Tanimoto_node_traits
{
 public:

  typedef MultipartiteNetworkFunctors::Tanimoto_coefficient  Similarity_function;

  class Unweighted_node_base
  {
   protected:
    
    Unweighted_node_base(){}
    
    template<class Neighbour_iterator> 
    double _compute_predicate_value_1( Neighbour_iterator rbegin, Neighbour_iterator rend, int size ) 
    { return double(1.0);  }

    template<class Neighbour_iterator>
    double _compute_predicate_value_2( Neighbour_iterator rbegin, Neighbour_iterator rend, int size ) 
    { return double(size+1);  }

  }; // Node_base

  class Weighted_node_base
  {
   protected:
    
    Weighted_node_base():_avg(-1.0),_squared_norm(-1.0){}
    
    template<class Neighbour_iterator> 
    double _compute_predicate_value_1( Neighbour_iterator rbegin, Neighbour_iterator rend, int size )
    { 
      if( _avg == double(-1.0) )
      {
        _avg = double(0.0);

        for( Neighbour_iterator nit = rbegin; nit != rend; nit++ )
          _avg += nit->weight(); 

        _avg = _avg/double(size);
      }

      return _avg;  
    }

    template<class Neighbour_iterator>
    double _compute_predicate_value_2( Neighbour_iterator rbegin, Neighbour_iterator rend, int size )
    { 
      if( _squared_norm == double(-1.0) )
      {
        _squared_norm = double(0.0);

        for( Neighbour_iterator nit = rbegin; nit != rend; nit++ )
          _squared_norm += (nit->weight())*(nit->weight()); 

        double avgg = _compute_predicate_value_1(rbegin,rend,size);
        _squared_norm += avgg*avgg;
      }

      return _squared_norm;  
    }

   private:

    double _avg, _squared_norm;

  }; // Weighted_node_base

}; // class Tanimoto_node_traits


class Minimum_weight_sum_node_traits
{
 public:

  typedef MultipartiteNetworkFunctors::Minimum_weight_sum  Similarity_function;

  class Unweighted_node_base
  {
   protected:
    
    Unweighted_node_base(){}
     
    template<class Neighbour_iterator>
    double _compute_predicate_value_1( Neighbour_iterator rbegin, Neighbour_iterator rend, int size ) 
    { return double(size);  }

    template<class Neighbour_iterator>
    double _compute_predicate_value_2( Neighbour_iterator rbegin, Neighbour_iterator rend, int size ) 
    { return double(-1.0);  }

  }; // Node_base

  class Weighted_node_base
  {
   protected:
    
    Weighted_node_base():_weight_sum(-1.0){}
     
    template<class Neighbour_iterator>
    double _compute_predicate_value_1( Neighbour_iterator rbegin, Neighbour_iterator rend, int size )
    { 
      if( _weight_sum == double(-1.0) )
      {
        _weight_sum = double(0.0);

        for( Neighbour_iterator nit = rbegin; nit != rend; nit++ )
          _weight_sum += nit->weight(); 
      }

      return _weight_sum;  
    }

    template<class Neighbour_iterator>
    double _compute_predicate_value_2( Neighbour_iterator rbegin, Neighbour_iterator rend, int size )
    {return double(-1.0);  }

   private:

    double _weight_sum;

  }; // Weighted_node_base

}; // class Minimum_weight_sum_node_traits

template< class PNT>
class Unweighted_edge_traits: public Edge_traits_base, PNT
{
 private:

  typedef PNT                                           Predicate_node_traits;
  typedef Edge_traits_base                              Traits_base;
  typedef typename Traits_base::Neighbour_element_base  Neighbour_base;
  typedef typename Traits_base::Primary_edge_base       Primary_edge_base; 

 public:

  typedef typename Predicate_node_traits::Unweighted_node_base  Node_base;
  typedef typename Predicate_node_traits::Similarity_function   Similarity_function;

  class Neighbour_element: public Neighbour_base
  {
    typedef Neighbour_base Base;

   public:
 
    Neighbour_element(int group, int ind, double w, int edge_ind ):
    Base(group,ind,edge_ind){}

    double weight(void) const
    { return double(1.0); }
  };

  typedef Is_greater_element<Neighbour_element>           Compare_pred;
  typedef std::set<Neighbour_element,Compare_pred>        Neighbour_set;
  typedef typename Neighbour_set::iterator                Neighbour_iterator;

  class Edge_base : public Primary_edge_base
  {
    typedef Primary_edge_base  Base;

    public: 

      Edge_base():Base(){}

      void set_weight( const double& w )
      {}

      double weight(void) const
      { return 1.0; }

  }; // Edge_base

}; // Unweighted_edge_traits


template< class PNT >
class Weighted_edge_traits : public Edge_traits_base, PNT
{
  typedef PNT                                           Predicate_node_traits;
  typedef Edge_traits_base                              Traits_base;
  typedef typename Traits_base::Neighbour_element_base  Neighbour_base;
  typedef typename Traits_base::Primary_edge_base       Primary_edge_base;

 public:

  typedef typename Predicate_node_traits::Weighted_node_base   Node_base;
  typedef typename Predicate_node_traits::Similarity_function  Similarity_function;

  class Neighbour_element: public Neighbour_base
  {
    typedef Neighbour_base Base;

   public:
 
    Neighbour_element(int group, int ind, double w, int edge_ind ):
    _weight(w), Base(group,ind,edge_ind){}

    double weight(void) const
    { return _weight; }

   private:

    double  _weight;

  }; // class Neighbour_element

  typedef Is_greater_element<Neighbour_element>            Compare_pred;
  typedef std::set<Neighbour_element, Compare_pred>        Neighbour_set;
  typedef typename Neighbour_set::iterator                 Neighbour_iterator;

  class Edge_base : public Primary_edge_base
  {
    typedef Primary_edge_base  Base;

    public: 

     Edge_base():_w(-1.0),Base(){}

     void set_weight( const double& w )
     { _w = w; }

     double weight(void) const
     { return _w; }

    private:

     double _w;
    
  }; // Edge_base

}; // Weighted_edge_traits


typedef Unweighted_edge_traits<Tanimoto_node_traits>            Standard_unweighted_edge_traits;
typedef Weighted_edge_traits<Tanimoto_node_traits>              Standard_weighted_edge_traits;
typedef Unweighted_edge_traits<Minimum_weight_sum_node_traits>  MSW_unweighted_edge_traits;
typedef Weighted_edge_traits<Minimum_weight_sum_node_traits>    MSW_weighted_edge_traits;

} // namespace MultipartiteNetwork


///////////////////////////////////////////////////////////////////////////
/////////////////////// Multipartite Network Types ////////////////////////
///////////////////////////////////////////////////////////////////////////

namespace MultipartiteNetwork{

template < class TRAITS >
class Multipartite_network_node : 
      public TRAITS::Edge_weight_traits::Node_base
{
 public:
  
  typedef TRAITS                                           Network_traits;
  typedef typename Network_traits::Edge_weight_traits      Edge_weight_traits; 
  typedef typename Edge_weight_traits::Neighbour_element   Neighbour_element; 
  typedef typename Edge_weight_traits::Neighbour_set       Neighbour_set;
  typedef typename Edge_weight_traits::Neighbour_iterator  Neighbour_iterator;
  
 public:

  Multipartite_network_node(){};

  void set_name(const std::string& name)
  { _name = name;}

  void set_group(int group)
  { _group = group;}

  void set_index(int index)
  { _index = index;}

  std::string name(void) const
  { return _name;}

  int group(void) const
  { return _group;}

  int index(void) const
  { return _index;}

  double compute_predicate_value_1()
  { return _compute_predicate_value_1( neighbours_begin(), neighbours_end(), number_of_neighbours() );  }

  double compute_predicate_value_2()
  { return _compute_predicate_value_2(neighbours_begin(), neighbours_end(), number_of_neighbours() );  }


  void insert_neighbour(Neighbour_element ne)
  { _neighbours.insert(ne); }

  int number_of_neighbours(void) const
  { return _neighbours.size(); }

  Neighbour_iterator find_neighbour( const Neighbour_element nelm ) const
  { return _neighbours.find(nelm); }

  Neighbour_iterator neighbours_begin()
  { return _neighbours.begin(); }

  Neighbour_iterator neighbours_end()
  { return _neighbours.end(); }

  bool is_marked()
  { return _mark;}

  void set_mark(bool m)
  { _mark=m; }
  
 private:

  int            _index;
  int            _group;
  Neighbour_set  _neighbours;
  std::string    _name;  
  bool           _mark;

}; // Multipartite_network_node

 
template< class TRAITS >
class Multipartite_network_edge: 
      public TRAITS::Edge_weight_traits::Edge_base
{
 public:

  typedef TRAITS                              Network_traits;
  typedef typename Network_traits::Node_type  Node;

 private:

  typedef typename Network_traits::Edge_weight_traits  Edge_weight_traits;

 public:

  Multipartite_network_edge()
  {

   _mark = false;
  }

  void set_adjacent_nodes( int group1, int ind1, int group2, int ind2 )
  { 
    _group1 = group1; 
    _group2 = group2;
    _ind1 = ind1; 
    _ind2 = ind2;
  }

  std::pair<int,int> source_node(void) const
  { return std::make_pair(_group1,_ind1); }

  std::pair<int,int> target_node(void) const
  { return std::make_pair(_group2,_ind2); }

  bool mark()
  { return _mark;}

  void set_mark(bool m)
  { _mark = m;}

 private:
  
  int _group1, _group2, _ind1,_ind2;


  bool _mark;

}; // Multipartite_network_edge

template< class TRAITS >
class Multipartite_network
{
 public:

  typedef TRAITS                                       Network_traits;
  typedef typename Network_traits::Node_type           Node;
  typedef typename Network_traits::Edge_type           Edge;
  typedef typename Network_traits::Neighbour_element   Neighbour_element;
  typedef typename Node::Neighbour_iterator            Neighbour_iterator;
 
 public:

  Multipartite_network():_number_of_nodes(-1), _is_connected(-1){}

  void extract_bipartite_network_from_csv( char *filename );

  void construct_bipartite_network_from_edge_list
  ( std::vector<std::string> nodes_1, std::vector<std::string> nodes_2,
    std::vector<int> endpoints_1, std::vector<int> endpoints_2, std::vector<double> weights); 

  void construct_bipartite_network_from_edge_list( int *n_1, char **nodes_1, int *n_2, char **nodes_2,
                                                   int *e, int *endpoints_1, int *endpoints_2, double *weights);

  void export_as_txt();

  // void Extract_multipartite_network(int i, std::string filename);

  int number_of_node_groups(void) const
  { return _nodes.size(); }

  int number_of_nodes(void)
  { 
    if(_number_of_nodes == -1)
    {
      _number_of_nodes = 0;

      for( int i=0; i<_nodes.size(); i++ )
        _number_of_nodes += _nodes[i].size();
    }

    return _number_of_nodes; 
  }

  int number_of_nodes_in_group(int group) const
  { return _nodes[group].size(); }

  int number_of_edges(void) const
  { return _edges.size(); }

  Node& node(int group, int index)
  { return _nodes[group][index]; }

//  const Node& node(int group, int index) const
//  { return _nodes[group][index]; }

  Edge& edge(int ind)
  { return _edges[ind]; }  

  void unmark_all_nodes()
  {
    for(int i=0; i<_nodes.size(); i++)
      for(int j=0; j<_nodes[i].size(); j++)
        _nodes[i][j].set_mark(false);
  }

  bool is_connected()
  {
    if(_is_connected != -1)
      return _is_connected;

    if(number_of_nodes() == 0)
      return false;

    std::queue< std::pair<int, int> > node_queue;
    unmark_all_nodes();

    node_queue.push(std::make_pair(0,0));
    
    while(node_queue.size() != 0)
    {
      std::pair<int, int> node_ip = node_queue.front();
      node_queue.pop();

      if( !node(node_ip.first,node_ip.second).is_marked() )
      {
        node(node_ip.first,node_ip.second).set_mark(true);

        for( Neighbour_iterator nit  = node(node_ip.first,node_ip.second).neighbours_begin(); 
                                nit != node(node_ip.first,node_ip.second).neighbours_end(); nit++ )
        {
          int g = nit->group(), i = nit->index();
          
          if( !node(g,i).is_marked() )
            node_queue.push(std::make_pair(g,i));
        }

      } // if( !node(node_ip.first,node_ip.second).is_marked() )

    } // while(node_queue.size() != 0)

    _is_connected = 1;

    for(int i=0; i<_nodes.size(); i++)
      for(int j=0; j<_nodes[i].size(); j++)
        if( !_nodes[i][j].is_marked() )
        {
          _is_connected = 0;
          return false;
        }

    return true;

  } // bool is_connected() 

  void print_edge(int ind, std::ostream &os = std::cout)
  {
    std::pair<int,int> snd = _edges[ind].source_node(),
                       tnd = _edges[ind].target_node();

    os << _nodes[snd.first][snd.second].name() << " - " 
       << _nodes[tnd.first][tnd.second].name() << std::endl;

  }
    

 private:

  int _number_of_nodes;
  std::vector< std::vector<Node> >  _nodes;
  std::vector<Edge>                 _edges;
  int _is_connected;

}; // class Multipartite_network



// Auxiliary type for the main clustering class.
// Does not appear in the traits class.

template <class TRAITS>
class Edge_cluster_node
{
 private:

  typedef TRAITS                                              Network_traits;
  typedef typename Network_traits::Edge_type                  Edge;
  typedef typename Network_traits::Node_type                  Node;
  typedef typename std::map< int, std::set<int>* >::iterator  Container_iterator;

 public:

  Edge_cluster_node( int cluster_index, const Edge& ed, int number_of_node_groups)
  {
    _nodes = new std::map< int, std::set<int>* >();
 
    std::pair<int,int> source_nd = ed.source_node(),
                       target_nd = ed.target_node();
   
    insert_network_node(source_nd.first, source_nd.second);
    insert_network_node(target_nd.first, target_nd.second);

    if( source_nd.first == target_nd.first )
      _sum_of_group_size_squares = 4;
    else
      _sum_of_group_size_squares = 2;

    _edge = ed;
    _cluster_index = cluster_index;
    _father_index  = cluster_index;
    _number_of_edges = 1;
    _number_of_nodes = 2;
    _depth = 0;
  }

  double cluster_partition_density(void) const
  {
    if( _number_of_edges+1-_number_of_nodes == 0 )
      return double(0.0);

    return double(_number_of_edges+1-_number_of_nodes)/
           ( (double(_number_of_nodes)*double(_number_of_nodes)) 
           - _sum_of_group_size_squares - double(2*_number_of_nodes)+double(2) );
  }

  bool is_representative(void) const
  { return _cluster_index == _father_index;}

  void set_father_index(int new_index)
  { _father_index = new_index; }

  int cluster_index(void) const
  { return _cluster_index; }

  int father_index(void) const
  { return _father_index; }

  bool insert_network_node( int group, int index )
  {
    std::set<int>* set_p;
    Container_iterator s_it;

    if( (s_it = _nodes->find(group)) == _nodes->end() )
      set_p = (*_nodes)[group] = new std::set<int>();
    else
      set_p = s_it->second;
 
    bool inserted = set_p->insert(index).second;

    if(inserted)
    {
      _number_of_nodes++;
      int size = set_p->size();
      _sum_of_group_size_squares += (double(size)*double(size) - double(size-1)*double(size-1)) ;       
    }

    return inserted;
  }

  std::map< int, std::set<int>* >*  node_container_pointer() const
  { return _nodes; }

  void  set_node_container_pointer( std::map< int, std::set<int>* >* cont )
  { _nodes = cont; }

  int depth(void) const
  { return _depth;}

  Edge edge(void) const
  { return _edge; }

  void set_number_of_edges( int num )
  { _number_of_edges = num; }

  void set_depth( int new_depth )
  { _depth = new_depth; }

  int number_of_nodes(void) const
  { return _number_of_nodes; }

  int number_of_edges(void) const
  { return _number_of_edges; }
 
  double sum_of_group_size_squares(void) const
  { return _sum_of_group_size_squares;}

  void set_number_of_nodes( int num )
  { _number_of_nodes = num; }

  void set_sum_of_group_size_squares( double sgsq )
  { _sum_of_group_size_squares = sgsq;}
  

 private:
 
  int _cluster_index, _father_index, _depth, _number_of_edges, _number_of_nodes;
  std::map< int, std::set<int>* >* _nodes;
  double _sum_of_group_size_squares; // In reality, always an integer
                                     // but double here to avoid overflow for large values.
                                     // This value is maintained for the reason
                                     // of a more efficient computation of
                                     // the partition density.
  Edge _edge; 

}; // class Edge_cluster_node


template <class TRAITS>
class Single_linkage_clustering
{
 private:
 
  typedef TRAITS                                         Network_traits;
  typedef Edge_cluster_node<Network_traits>              Union_find_node;
  typedef typename Network_traits::Multipartite_network  Multipartite_network;
  typedef typename Network_traits::Node_type             Node;
  typedef typename Network_traits::Edge_type             Edge;
  typedef typename Node::Neighbour_element               Neighbour_element;
  typedef typename Node::Neighbour_iterator              Neighbour_iterator;

  typedef typename Network_traits::Similarity_function   Similarity_function;

  struct Edge_link
  {
    int edge_index_1, edge_index_2;
    double similarity;

  }; // Edge_link

  struct Merge_event
  {
    std::map<int,std::set<int> >  groupings;
    double similarity, partition_density;
    int added_links;
  };


  struct Is_greater_edge_link
  {
    bool operator()( const Edge_link& el1, const Edge_link& el2 )
    {
      if(el1.similarity > el2.similarity)
        return true;

      if(el1.similarity < el2.similarity)
        return false;

      if(el1.edge_index_1 < el2.edge_index_1)
        return true;

      if(el1.edge_index_1 > el2.edge_index_1)
        return false;

      if(el1.edge_index_2 < el2.edge_index_2)
        return true;

      return false;
    }
  }; // Is_greater_edge_link


  struct Is_greater_cluster
  {
    bool operator()( const std::vector<Edge>& ev1, const std::vector<Edge>& ev2 )
    {
      if(ev2.size() == 0 && ev1.size() != 0)
        return true;

      if(ev1.size() == 0)
        return false;

      if( ev1[0].source_node().first < ev2[0].source_node().first )
        return true;

      if( ev1[0].source_node().first > ev2[0].source_node().first )
        return false;

      if( ev1[0].source_node().second < ev2[0].source_node().second )
        return true;

      if( ev1[0].source_node().second > ev2[0].source_node().second )
        return false;

      if( ev1[0].target_node().first < ev2[0].target_node().first )
        return true;

      if( ev1[0].target_node().first > ev2[0].target_node().first )
        return false;

      if( ev1[0].target_node().second < ev2[0].target_node().second )
        return true;

      return false;
    }

  }; // Is_greater_cluster

  struct Merge_tree_node
  {
    std::vector<int> children;
    std::string name;
    int edge_weight, total_depth, parent, edges_in_subtree, edge_range_begin, edge_range_end;

    int number_of_children()
    { return children.size(); }

  };

 private:

  Union_find_node& find(int cluster_index)
  {
    Union_find_node &ufn = _union_find_ds[cluster_index];

    if( ufn.is_representative() )
      return ufn;

    Union_find_node &new_father_node = find(ufn.father_index());
    ufn.set_father_index( new_father_node.cluster_index() );

    return new_father_node;

  } // Union_find_node& find(int cluster_index)


  double unionn( int cluster_index1, int cluster_index2 )
  {
    Union_find_node &ufn1 = find(cluster_index1),
                    &ufn2 = find(cluster_index2);

    if( ufn1.depth() >= ufn2.depth() )
      return unionn( ufn1, ufn2 );
    else
      return unionn( ufn2, ufn1 );           

  } // double unionn( int ..., int ... )


  double unionn( Union_find_node& ufn1, Union_find_node& ufn2 )
  {
    ufn2.set_father_index(ufn1.cluster_index());
 
    if(ufn1.depth() == ufn2.depth())
      ufn1.set_depth( ufn1.depth()+1 );
 
    ufn1.set_number_of_edges( ufn1.number_of_edges()+ufn2.number_of_edges() );

    Union_find_node  *big, *small;

    if( ufn1.number_of_nodes() >= ufn2.number_of_nodes() )
    {
      big = &ufn1;
      small = &ufn2;
    }
    else
    {
      big = &ufn2;
      small = &ufn1;
    }

    std::map<int, std::set<int>* >* nd_ptr = small->node_container_pointer();
    typename std::map<int, std::set<int>* >::iterator it;

    for( it = nd_ptr->begin(); it != nd_ptr->end(); it++ )
    {
      std::set<int>* ns = it->second;
      typename std::set<int>::iterator sit;

      for( sit = ns->begin(); sit != ns->end(); sit++ )
        big->insert_network_node(it->first, *sit);

      ns->clear();
      delete ns;
    }

    nd_ptr->clear();
    delete nd_ptr;

    ufn1.set_node_container_pointer(big->node_container_pointer());
    ufn1.set_number_of_nodes(big->number_of_nodes());
    ufn1.set_sum_of_group_size_squares(big->sum_of_group_size_squares());

    return ufn1.cluster_partition_density();

  } // double unionn( Union_find_node& ..., Union_find_node& ... )


  // Only for debugging

  int find_smallest_edge_in_cluster( const Union_find_node &ufn )
  {
    int smallest = -1;

    for( int i=0; i<_union_find_ds.size(); i++ )
      if( _union_find_ds[i].father_index() == ufn.cluster_index() )
        if(_union_find_ds[i].edge().index() < smallest || smallest == -1 )
          smallest = _union_find_ds[i].edge().index();

    return smallest;

  } // find_smallest_edge_in_cluster( const Union_find_node &ufn )

  int number_of_clusters(void)
  {
    int num=0;

    for(int i=0; i<_union_find_ds.size(); i++)
    {
      find(i);
    
      if( _union_find_ds[i].is_representative() )
        num++;
    }

    return num;
  } 

  void print_cluster( const Union_find_node &ufn )
  {
    for( int i=0; i<_union_find_ds.size(); i++ )
    {
      find(i);

      if( _union_find_ds[i].father_index() == ufn.cluster_index() )
      {
        Edge e = _union_find_ds[i].edge();

        std::pair<int,int> snd = e.source_node(),
                           tnd = e.target_node();

        Node nd1 = _mnw.node(snd.first, snd.second),
             nd2 = _mnw.node(tnd.first, tnd.second);
         
        std::cerr << nd1.name() 
                  << " - " << nd2.name() << std::endl;
      }

    }

  } // void print_cluster( ... )

  void _arrange_edges_in_subtree(int index, int e1, int e2);

  void _arrange_edges_according_to_cluster_hierarchy();


  void _print_merge_tree(Merge_tree_node root, std::ofstream &outf)
  {
    if(root.number_of_children() == 0)
      outf << root.name << ":" << root.edge_weight;
    else
    {
      outf << "(";

      for( int i=0; i<root.number_of_children(); i++)
      { 
        _print_merge_tree(_merge_nodes[root.children[i]], outf);

        if(i < root.number_of_children()-1 )
          outf << ",";
        else
          outf << ")";
      }

      if(root.name.size() != 0 && root.edge_weight > 0 )
        outf << root.name << ":" << root.edge_weight;
      else if( root.edge_weight > 0)
        outf << root.edge_weight;
      else
        outf << root.name;
    }
  }

 public:

  Single_linkage_clustering():_max_partition_density(0.0)
  {}

  void operator()( Multipartite_network &mnw )
  {

    Similarity_function     similarity;
    std::vector<Edge_link>  links;
    double partition_density(0.0);

    _max_partition_density = 0.0;

    _mnw = mnw;


    for( int i=0; i<mnw.number_of_node_groups(); i++ )
    {

      for( int j=0; j<mnw.number_of_nodes_in_group(i); j++ )
      {
        Node nd = mnw.node(i,j);
        Neighbour_iterator nit1, nit2;
        int count1=0, count2=0;

        for( nit1 = nd.neighbours_begin(); nit1 != nd.neighbours_end(); nit1++ )
        {
          for( nit2 = nd.neighbours_begin(); nit2 != nit1; nit2++ )
          {          
            const Neighbour_element nd1 = *nit1, 
                                    nd2 = *nit2;
            Edge_link el;

            el.edge_index_1 = nd2.edge_index();
            el.edge_index_2 = nd1.edge_index();

            Node &neigh1 = mnw.node(nd1.group(), nd1.index()),
                 &neigh2 = mnw.node(nd2.group(), nd2.index());

            el.similarity = similarity(neigh1, neigh2);            
  
            links.push_back(el);            

          } // for( nit2 = nd.neihgbours.begin(); ... )

        } // for( nit1 = nd.neighbours.begin(); ... )

      } // for( int j=0; ... )

    } // for( int i=0; ... )

    std::sort(links.begin(), links.end(), Is_greater_edge_link() );

    for( int i=0; i<mnw.number_of_edges(); i++ )
      _union_find_ds.push_back( Union_find_node(i,mnw.edge(i),mnw.number_of_node_groups()) ); 

    for(int i=0; i<mnw.number_of_edges(); i++ )
    {
      Merge_tree_node mnd;

      std::stringstream edge_name;

      std::pair<int, int>  n11 = mnw.edge(i).source_node(),
                           n12 = mnw.edge(i).target_node();

      edge_name << n11.first << "." << n11.second << "." << mnw.node(n11.first,n11.second).name();
      edge_name << "||";
      edge_name << n12.first << "." << n12.second << "." << mnw.node(n12.first,n12.second).name(); 

      mnd.name = edge_name.str();
      mnd.edge_weight = mnw.number_of_edges();
      mnd.total_depth = mnw.number_of_edges();
      mnd.parent = -1;
      mnd.edges_in_subtree = 1;

      _merge_nodes.push_back(mnd); 
    }

    std::vector<int> latest_edge_merge_tree_node;

    for(int i=0; i<mnw.number_of_edges(); i++ )
      latest_edge_merge_tree_node.push_back(i);
    

    int count=0, added_links = 0;

    _cluster_graph_x.push_back(mnw.number_of_edges());
    _cluster_graph_y.push_back(partition_density);
     
    while( added_links < mnw.number_of_edges()-1 && count < links.size() )
    {
      std::vector<Edge_link> links_to_add, links_that_were_added;

      links_to_add.push_back(links[count]);
      count++;

      while( count < links.size() && 
             links[count].similarity == links_to_add.back().similarity )
      {
        links_to_add.push_back(links[count]);
        count++;
      } 

      std::map<int, std::list<int> > list_map;

      for( int i=0; i< links_to_add.size(); i++ )
      {
        Edge_link el = links_to_add[i];

        Union_find_node ufn1 = find(el.edge_index_1),
                        ufn2 = find(el.edge_index_2);

        if( ufn1.cluster_index() != ufn2.cluster_index() )
        {
          latest_edge_merge_tree_node[el.edge_index_1]=
          latest_edge_merge_tree_node[ufn1.cluster_index()];

          latest_edge_merge_tree_node[el.edge_index_2]=
          latest_edge_merge_tree_node[ufn2.cluster_index()];

          if( list_map.find(ufn1.cluster_index()) == list_map.end() )
          {
            std::list<int> l_tmp;
            l_tmp.push_back(latest_edge_merge_tree_node[ufn1.cluster_index()]);
            list_map[ufn1.cluster_index()] = l_tmp;
          }

          if( list_map.find(ufn2.cluster_index()) == list_map.end() )
          {
            std::list<int> l_tmp;
            l_tmp.push_back(latest_edge_merge_tree_node[ufn2.cluster_index()]);
            list_map[ufn2.cluster_index()] = l_tmp;
          }


          double c1_pd = ufn1.cluster_partition_density(),
                 c2_pd = ufn2.cluster_partition_density();

          partition_density = partition_density - double(ufn1.number_of_edges())*c1_pd 
                                                - double(ufn2.number_of_edges())*c2_pd;

          int ci1 = ufn1.cluster_index(),
              ci2 = ufn2.cluster_index();

          std::list<int> &l1 = list_map[ufn1.cluster_index()], 
                         &l2 = list_map[ufn2.cluster_index()];

          unionn(ufn1.cluster_index(),ufn2.cluster_index());
          
          Union_find_node ufn = find(ufn1.cluster_index());

          if(ufn.cluster_index() == ci1)
            l1.splice(l1.end(),l2);
          else
            l2.splice(l2.end(),l1);

          partition_density = partition_density + double(ufn.number_of_edges())*ufn.cluster_partition_density();

          Edge_link event;

          event.edge_index_1 = ufn1.cluster_index();
          event.edge_index_2 = ufn2.cluster_index();


          links_that_were_added.push_back(event);

          added_links++;

        } // if( ufn1.cluster_index() == ufn2.cluster_index() )       

      } // for( int i=0; i< links_to_add.size(); i++ )

      if(_max_partition_density < partition_density )
        _max_partition_density = partition_density;

      Merge_event merged_to;

      std::map<int, std::set<int> > merge_tree_node_groupings;

      for(int i=0; i<links_that_were_added.size(); i++ )
      {
        Edge_link cl = links_that_were_added[i]; 

        Union_find_node ufn = find(cl.edge_index_1);
   
        if( merged_to.groupings.find(ufn.cluster_index()) == merged_to.groupings.end() )
          merged_to.groupings[ufn.cluster_index()] = std::set<int>();

        if( !_union_find_ds[cl.edge_index_1].is_representative() )
          merged_to.groupings[ufn.cluster_index()].insert(cl.edge_index_1);

        if( !_union_find_ds[cl.edge_index_2].is_representative() )
          merged_to.groupings[ufn.cluster_index()].insert(cl.edge_index_2);

      } // for(int i=0; i<links_that_were_added.size(); i++ )

      typename std::map< int, std::list<int> >::iterator git;

      for(git = list_map.begin(); git != list_map.end(); git++ )
        if(git->second.size() > 0)
        {
          Merge_tree_node mtn;

          std::set<int> children_set;

          for( typename std::list<int>::iterator lit=git->second.begin(); lit!= git->second.end(); lit++ )
            children_set.insert(*lit);

          for( typename std::set<int>::iterator sit=children_set.begin(); sit!= children_set.end(); sit++ )
            mtn.children.push_back(*sit);

          std::stringstream sstr;

          sstr << "Node ";
          sstr << _merge_nodes.size();

          mtn.name = sstr.str();
          mtn.edge_weight = mnw.number_of_edges()-added_links;
          mtn.total_depth = mnw.number_of_edges()-added_links;
          mtn.parent = -1;
          mtn.edges_in_subtree=0;

          _merge_nodes.push_back(mtn);

          for(int hh = 0; hh < mtn.number_of_children(); hh++ )
          {
            _merge_nodes[mtn.children[hh]].edge_weight = _merge_nodes[mtn.children[hh]].edge_weight - mtn.edge_weight;
            _merge_nodes[mtn.children[hh]].parent = _merge_nodes.size()-1;
            _merge_nodes.back().edges_in_subtree += _merge_nodes[mtn.children[hh]].edges_in_subtree;
          }

          latest_edge_merge_tree_node[git->first] = _merge_nodes.size()-1;

      } // if(git->second.size() > 0)

      if(links_that_were_added.size() != 0 )
      {
        _cluster_graph_x.push_back(mnw.number_of_edges()-added_links);
        _cluster_graph_y.push_back(partition_density);
      }
      
      merged_to.similarity = links_to_add[0].similarity;
      merged_to.partition_density = partition_density;
      merged_to.added_links = added_links;

      _merge_events.push_back(merged_to);

    } // while( added_links < mnw.number_of_edges()-1 )

    _merge_nodes.back().edge_weight=0;
    _arrange_edges_according_to_cluster_hierarchy();

  } // void operator()( const Multipartite_network &mnw )

  void print_merge_tree(char* filename)
  {
    std::ofstream outf(filename);

    _print_merge_tree(_merge_nodes[_merge_nodes.size()-1], outf);
    outf << ";" << std::endl;
  }

  void compute_max_partition_density_clustering(void)
  {
    if(_mnw.number_of_nodes() == 0)
      return;

    _union_find_ds.clear();

    for( int i=0; i<_mnw.number_of_edges(); i++ )
      _union_find_ds.push_back( Union_find_node(i,_mnw.edge(i),_mnw.number_of_node_groups()) );


    // Find the event after which the 
    // max partition density is achieved.

    int lm=0;

    while(lm<_merge_events.size())
      if(_merge_events[lm].partition_density == _max_partition_density)
        break;
      else
        lm++;

    if( lm == _merge_events.size() )
      return;

    for( int i=0; i<=lm; i++ )
    {
      for( std::map< int, std::set<int> >::iterator mit = _merge_events[i].groupings.begin();
           mit != _merge_events[i].groupings.end(); mit++ )
      {
        int father = mit->first;

        typename std::set<int>::iterator pr_it;

        for( pr_it = mit->second.begin(); pr_it != mit->second.end(); pr_it++ )
          unionn(father,*pr_it);
 
      } // for( std::map< int, std::set<int> >::iterator mit = ... )

    } // for( int i=0; i<=lm; i++ )

  } // void print_max_partition_density_clustering(void) const

  void print_current_cluster_status( std::ostream &oss = std::cout )
  {
    std::vector< std::vector<Edge> >  edge_clusters;
    int clusters_found = 0;

    // This is maybe not the most handy way to print the clusters
    // yet for debugging reasons it is imperative to sort the clusters
    // according to the indices of their edges.

    edge_clusters.assign( _union_find_ds.size(), std::vector<Edge>() );

    for( int i=0; i<_union_find_ds.size(); i++ )
    {
      Union_find_node fn = find(i);

      if(_union_find_ds[i].is_representative())
        clusters_found++;

      edge_clusters[fn.father_index()].push_back(_union_find_ds[i].edge());
    }

    oss << " There were " << clusters_found << " clusters found." << std::endl;
    oss << " The maximum partition density that was observed was: " 
        << max_partition_density() << std::endl;
 
    int cluster_no = 1;

    std::sort( edge_clusters.begin(), edge_clusters.end(), Is_greater_cluster() );

    for( int i=0; i <edge_clusters.size(); i++ )
      if( edge_clusters[i].size() > 0 )
      {
 
        oss << std::endl << std::endl << std::endl;
        oss << " *********************************** " << std::endl;
        oss << " Cluster #" << cluster_no << std::endl;
        oss << " *********************************** " << std::endl;
        oss << std::endl << " Edges: " << std::endl;
        oss << " ------ " << std::endl;

        cluster_no++;

        std::vector<Edge> cedges = edge_clusters[i];
        std::map<int, std::set<int>*> * nd_ptr;
        int num_of_nodes;

        for( int j=0; j<cedges.size(); j++ )
        {


          std::pair<int,int> source_nd = cedges[j].source_node(),
                           target_nd = cedges[j].target_node();

          if( _union_find_ds[cedges[j].index()].is_representative() )
          {
            nd_ptr = _union_find_ds[cedges[j].index()].node_container_pointer();
            num_of_nodes = _union_find_ds[cedges[j].index()].number_of_nodes();
          }

          oss << _mnw.node(source_nd.first, source_nd.second).name() << " - " 
              << _mnw.node(target_nd.first, target_nd.second).name() << std::endl;

        } // for( int j=0; j<cedges.size(); j++ )

        oss << std::endl << std::endl << " Nodes: " << std::endl;
        oss << " ------ " << std::endl;

        for( typename std::map<int, std::set<int>* >::iterator mit  = nd_ptr->begin();
                                                               mit != nd_ptr->end(); mit++ )
        {
          oss << std::endl << " Group #" << mit->first << std::endl;
          oss << " --------" << std::endl;

          for( typename std::set<int>::iterator sit  = mit->second->begin(); 
                                                sit != mit->second->end(); sit++ )
            oss << _mnw.node(mit->first, *sit).name() << std::endl;   
        }
      
      } // if( edge_clusters[i] > 0 )

  } // void print_cluster_status(void) const

  double max_partition_density(void) const
  { return double(2.0)*_max_partition_density/double(_mnw.number_of_edges()); }

  template< class OutputIterator >
  void extract_cluster_structure( OutputIterator o_edges, OutputIterator p_nodes, OutputIterator depths, 
                                  OutputIterator edge_ranges_begin, OutputIterator edge_ranges_end)
  {
    for(int i=0; i<_hierarchy_sorted_edges.size(); i++ )
      *o_edges++ = _hierarchy_sorted_edges[i]+1;

    for(int i=0; i<_merge_nodes.size(); i++ )
    {
      *p_nodes++ = _merge_nodes[i].parent+1;
      *depths++ = _merge_nodes[i].total_depth;
      *edge_ranges_begin++ = _merge_nodes[i].edge_range_begin+1;
      *edge_ranges_end++ = _merge_nodes[i].edge_range_end+1;
    }

  }

  template< class OutputIteratorInt, class OutputIteratorDouble >
  void extract_cluster_graph( OutputIteratorInt cg_x, OutputIteratorDouble cg_y)
  {
    for(int i=0; i<_cluster_graph_x.size(); i++ )
    {
      *cg_x++ = _cluster_graph_x[i];
      *cg_y++ = _cluster_graph_y[i];
    }
  }

 private:

    Multipartite_network _mnw;
    double _max_partition_density;
    std::vector<Union_find_node> _union_find_ds;
    std::vector< Merge_event > _merge_events; // This is a vector that maintains the history of 
                                              // the groups of edges that were merged together
                                              // at each step of the clustering. More than one pair of edges
                                              // may get merged in a single step of the clustering
                                              // process if they have the same similarity value.

    std::vector< Merge_tree_node > _merge_nodes;
    std::vector<int> _cluster_graph_x;
    std::vector<double> _cluster_graph_y;
    std::vector<int> _hierarchy_sorted_edges;


}; // class Single_linkage_clustering

} // MultipartiteNetwork

///////////////////////////////////////////////////////////////////////////
//////////////// Multipartite Network Types Implementation ////////////////
///////////////////////////////////////////////////////////////////////////


namespace MultipartiteNetwork{

template< class TRAITS >
void Single_linkage_clustering<TRAITS>::_arrange_edges_in_subtree(int index, int e1, int e2)
{
  _merge_nodes[index].edge_range_begin = e1;
  _merge_nodes[index].edge_range_end = e2;

  if(e1==e2)
  {
    _hierarchy_sorted_edges[e1] = index;

    _mnw.edge(index).set_mark(true);
    return;
  } 

  int r1 = e1, r2;

  Merge_tree_node mtn = _merge_nodes[index];

  for(int i=0; i<mtn.number_of_children(); i++)
  {
    r2 = r1+_merge_nodes[mtn.children[i]].edges_in_subtree;
    _arrange_edges_in_subtree(mtn.children[i], r1, r2-1);
    r1 = r2;
  }
}

template< class TRAITS >
void Single_linkage_clustering<TRAITS>::_arrange_edges_according_to_cluster_hierarchy()
{
  _hierarchy_sorted_edges.assign(_mnw.number_of_edges(),-1);

  int range_begin = 0;

  for( int i = 0; i < _merge_nodes.size(); i++ )
    if( _merge_nodes[i].parent == -1 )
    {
      int tree_root = i;

      _arrange_edges_in_subtree(tree_root,range_begin,range_begin + _merge_nodes[tree_root].edges_in_subtree-1);
      range_begin += _merge_nodes[tree_root].edges_in_subtree;
    }

} // void _arrange_edges_according_to_cluster_hierarchy()


template< class TRAITS >
void Multipartite_network<TRAITS>::construct_bipartite_network_from_edge_list
( std::vector<std::string> nodes_1, std::vector<std::string> nodes_2,
  std::vector<int> endpoints_1, std::vector<int> endpoints_2, std::vector<double> weights)  
{
  int _number_of_node_groups=2;

  _nodes.push_back(std::vector<Node>());
  _nodes.push_back(std::vector<Node>());

  for( int i = 0; i < nodes_1.size(); i++ )
  {
    Node nd;
    nd.set_name(nodes_1[i]);
    nd.set_index(i);
    nd.set_group(0);
    _nodes[0].push_back(nd);
  }

  for( int i = 0; i < nodes_2.size(); i++ )
  {
    Node nd;
    nd.set_name(nodes_2[i]);
    nd.set_index(i);
    nd.set_group(1);
    _nodes[1].push_back(nd);
  }

  for(int i = 0; i < weights.size(); i++)
  {
    Edge e;

    e.set_adjacent_nodes(1,endpoints_2[i]-1,0,endpoints_1[i]-1);
    e.set_weight(weights[i]);
    e.set_index(i);

    _edges.push_back(e);

    Neighbour_element nelm1(0,endpoints_1[i]-1,weights[i],i), 
                      nelm2(1,endpoints_2[i]-1,weights[i],i);
      
    _nodes[0][endpoints_1[i]-1].insert_neighbour(nelm2);
    _nodes[1][endpoints_2[i]-1].insert_neighbour(nelm1);
  }

  _number_of_nodes = _nodes[0].size() + _nodes[1].size();

} // void construct_bipartite_network_from_edge_list( ... )


// Version that uses the old interface with R

template< class TRAITS >
void Multipartite_network<TRAITS>::construct_bipartite_network_from_edge_list
( int *n_1, char **nodes_1, int *n_2, char **nodes_2,
  int *e, int *endpoints_1, int *endpoints_2, double *weights)  
{
  int _number_of_node_groups=2;

  _nodes.push_back(std::vector<Node>());
  _nodes.push_back(std::vector<Node>());

  for( int i = 0; i < *n_1; i++ )
  {
    Node nd;
    nd.set_name(std::string(nodes_1[i]));
    nd.set_index(i);
    nd.set_group(0);
    _nodes[0].push_back(nd);
  }

  for( int i = 0; i < *n_2; i++ )
  {
    Node nd;
    nd.set_name(std::string(nodes_2[i]));
    nd.set_index(i);
    nd.set_group(1);
    _nodes[1].push_back(nd);
  }

  for(int i = 0; i < *e; i++)
  {
    Edge e;

    e.set_adjacent_nodes(1,endpoints_2[i]-1,0,endpoints_1[i]-1);
    e.set_weight(weights[i]);
    e.set_index(i);

    _edges.push_back(e);

    Neighbour_element nelm1(0,endpoints_1[i]-1,weights[i],i), 
                      nelm2(1,endpoints_2[i]-1,weights[i],i);
      
    _nodes[0][endpoints_1[i]-1].insert_neighbour(nelm2);
    _nodes[1][endpoints_2[i]-1].insert_neighbour(nelm1);
  }

  _number_of_nodes = _nodes[0].size() + _nodes[1].size();

} // void construct_bipartite_network_from_edge_list( ... )

template< class TRAITS >
void Multipartite_network<TRAITS>::extract_bipartite_network_from_csv( char *filename )
{
  int _number_of_node_groups = 2;

  _nodes.push_back(std::vector<Node>());
  _nodes.push_back(std::vector<Node>());

  std::ifstream in(filename);

  if( !( in.is_open() && in.good() ) )
  {
    std::cerr << " ERROR: There was a problem with opening the file with the name list ... aborting " <<std::endl;
    std::exit(-1);
  }

  // Read the first row, the one that contains the species names

  std::vector<std::string> names;
  std::string line;
  char a;

  std::getline(in,line);

  int c=0;
  int node_index1 = 0, 
      node_index2 = 0;

  while(c < line.size() )
  {
    Node nd;
    std::string str;
    a = line[c]; 

    while(a != ',' && a != '\r' && c<line.size() )
    {
      str.push_back(a);
      c++;
      a = line[c];
    }

    if( str.size() > 0 )
    {

      nd.set_name(str);
      nd.set_index(node_index1);
      nd.set_group(0);

      node_index1++;
      _nodes[0].push_back(nd);
    }

    c++;
  }

  // Read the rest of the matrix, adding the network edges on each line.

  std::string str;
  int edge_index=0;

  while( in.good() )
  {
    line.clear();
    Node nd;
    int count=0, start=0;
    std::getline(in,line);

    str.clear();
      
    // Get the node name for this line

    while( line[start] != ',' && line[start] != '\r' && start<line.size() )
    {
      str.push_back(line[start]);
      start++;
    }

    if( str.size() > 0 )
    {

      nd.set_name(str);
      nd.set_index(node_index2);
      nd.set_group(1);

      _nodes[1].push_back(nd);
    }

    while( (    line[start] == ','  || line[start] == '\r' 
             || line[start] == '"' || isspace(line[start]) ) && start<line.size() )
    {
      start++;
    }


    // Get the edge set for this vertex

    int i=start;

    while(i<line.size())
    {
      std::string alph_num;
      double num;

      a = line[i];
     
      while( a>= '0' && a<= '9' && !isspace(a) )
      {
        alph_num.push_back(a);
        i++;

        a = line[i];
      }

      num = ::atof(alph_num.c_str());

      if( num > double(0.0) )
      {
        Edge e;

        e.set_adjacent_nodes(1,_nodes[1].size()-1,0,count);
        e.set_weight(num);
        e.set_index(edge_index);

        _edges.push_back(e);

        Neighbour_element nelm1(0,count,num,edge_index), 
                          nelm2(1,node_index2,num,edge_index);




        
        _nodes[0][count].insert_neighbour(nelm2);
        _nodes[1].back().insert_neighbour(nelm1);

        edge_index++;
        count++;
      }
      else if( num == double(0.0) )
        count++;
      else
      {
        std::cerr << " ERROR: A negative weight was read for an edge of the network ... aborting " <<std::endl;
        std::exit(-1);
      }        

      

      while( (    line[i] == ','  || line[i] == '\r' 
               || line[i] == '"' || isspace(line[i]) ) && i<line.size() )
      {
        i++;
      }

    } // while(i<line.size())

    node_index2++;

  } // while( in.good() )


  _number_of_nodes = _nodes[0].size() + _nodes[1].size();

} // void extract_bipartite_network_from_csv( char *filename )


template< class TRAITS >
void Multipartite_network<TRAITS>::export_as_txt()
{
  std::cout << std::endl << std::endl;
  std::cout << " The network has " << _number_of_nodes << " nodes and " 
            << _edges.size() << " edges. " << std::endl;

  std::cout << " The first group of nodes has " << _nodes[0].size() 
            << " elements, while the second group has " << _nodes[1].size() << " elements " 
            << std::endl << std::endl; 
            
  std::cout << std::endl << " The nodes of the network are the following: " << std::endl;

  for( int k=0; k< _nodes.size(); k++ )
  {

    std::cout << std::endl << " Group #:" << k << std::endl; 
    std::cout << " *******" << std::endl; 

    for( int i=0; i<_nodes[k].size(); i++ )
      std::cout << " Index:" << i << " || " << _nodes[k][i].name() << std::endl; 
  }

  std::cout << std::endl << " The edges of the network are the following: " << std::endl;
  std::cout << " ******************************************* " << std::endl;

  for( int i=0; i< _edges.size(); i++ )
  {
    std::pair<int, int> source_nd = _edges[i].source_node(),
                        target_nd = _edges[i].target_node();

    // Next "std::cout" to be removed:
    std::cout << "Edge #" << _edges[i].index() << " ";


    std::cout << _nodes[source_nd.first][source_nd.second].name() << " <-> " 
              << _nodes[target_nd.first][target_nd.second].name() 
              << " has weight: " << _edges[i].weight() << std::endl; 
  }

} // void export_network_as_txt( char *filename )

} // MultipartiteNetwork

///////////////////////////////////////////////////////////////////////////
/////////////////////// Multipartite Networks Traits //////////////////////
///////////////////////////////////////////////////////////////////////////


template< class EWT >
struct Multipartite_network_traits
{
  typedef EWT                                                               Edge_weight_traits;
  typedef Multipartite_network_traits<EWT>                                  Self;
  typedef typename MultipartiteNetwork::Multipartite_network_node<Self>     Node_type;
  typedef typename MultipartiteNetwork::Multipartite_network_edge<Self>     Edge_type;
  typedef typename Edge_weight_traits::Neighbour_element                    Neighbour_element; 
  typedef typename MultipartiteNetwork::Multipartite_network<Self>          Multipartite_network; 
  typedef typename MultipartiteNetwork::Single_linkage_clustering<Self>     Single_linkage_clustering; 

  typedef typename Edge_weight_traits::Similarity_function                  Similarity_function;  

}; // Multipartite_network_traits


typedef Multipartite_network_traits<MultipartiteNetwork::Standard_unweighted_edge_traits>  
                                   Standard_unweighted_multipartite_network_traits;

typedef Multipartite_network_traits<MultipartiteNetwork::Standard_weighted_edge_traits>  
                                   Standard_weighted_multipartite_network_traits;

typedef Multipartite_network_traits<MultipartiteNetwork::MSW_unweighted_edge_traits>  
                                   MSW_unweighted_multipartite_network_traits;

typedef Multipartite_network_traits<MultipartiteNetwork::MSW_weighted_edge_traits>  
                                   MSW_weighted_multipartite_network_traits;


///////////////////////////////////////////////////////////////////////////
///////////////////////////// CPP File Code ///////////////////////////////
///////////////////////////////////////////////////////////////////////////



template< class TRAITS >
Rcpp::List construct_community_hierarchy_internal( std::vector<std::string> &nodes_1, std::vector<std::string>& nodes_2, 
                                                   std::vector<int> &endpoints_1, std::vector<int>& endpoints_2, 
                                                   std::vector<double> &weights )
{
  typedef TRAITS  Network_traits;

  typedef typename Network_traits::Multipartite_network       Multipartite_network;
  typedef typename Network_traits::Node_type                  Node;
  typedef typename Node::Neighbour_element                    Neighbour_element;
  typedef typename Node::Neighbour_iterator                   Neighbour_iterator;
  typedef typename Network_traits::Single_linkage_clustering  Single_linkage_clustering;

  Multipartite_network  um_network;
  Single_linkage_clustering      clustering_method;

  um_network.construct_bipartite_network_from_edge_list(nodes_1,nodes_2,endpoints_1,endpoints_2,weights);

  clustering_method(um_network);
  clustering_method.compute_max_partition_density_clustering();

 
  std::vector<int> ordered_edges, parent_nodes, depths, edge_ranges_begin, 
                   edge_ranges_end, cluster_graph_x; 

  std::vector<double> cluster_graph_y;

  clustering_method.extract_cluster_structure(std::back_inserter(ordered_edges),std::back_inserter(parent_nodes),
                                              std::back_inserter(depths),std::back_inserter(edge_ranges_begin), 
                                              std::back_inserter(edge_ranges_end));

  clustering_method.extract_cluster_graph(std::back_inserter(cluster_graph_x),std::back_inserter(cluster_graph_y));

  Rcpp::IntegerVector OrderedEdges, ParentNodes, Depths, EdgeRangesBegin, 
                      EdgeRangesEnd, ClusterGraphX; 

  Rcpp::NumericVector ClusterGraphY;

  for(int i=0; i<ordered_edges.size(); i++)
    OrderedEdges.push_back(ordered_edges[i]);

  for(int i=0; i<parent_nodes.size(); i++)
    ParentNodes.push_back(parent_nodes[i]);

  for(int i=0; i<depths.size(); i++)
    Depths.push_back(depths[i]);

  for(int i=0; i<edge_ranges_begin.size(); i++)
    EdgeRangesBegin.push_back(edge_ranges_begin[i]);

  for(int i=0; i<edge_ranges_end.size(); i++)
    EdgeRangesEnd.push_back(edge_ranges_end[i]);

  for(int i=0; i<cluster_graph_x.size(); i++)
    ClusterGraphX.push_back(cluster_graph_x[i]);

  for(int i=0; i<cluster_graph_y.size(); i++)
    ClusterGraphY.push_back(cluster_graph_y[i]);

  return Rcpp::List::create( Rcpp::Named("A") = OrderedEdges,
                             Rcpp::Named("B") = ParentNodes,
                             Rcpp::Named("C") = Depths,
                             Rcpp::Named("D") = EdgeRangesBegin,
                             Rcpp::Named("E") = EdgeRangesEnd,
                             Rcpp::Named("F") = ClusterGraphX,
                             Rcpp::Named("G") = ClusterGraphY);

} // void cluster_network_internal(...)

// [[Rcpp::export]]
Rcpp::List construct_community_hierarchy( Rcpp::CharacterVector cmp_flags, Rcpp::CharacterVector nds1, 
                                          Rcpp::CharacterVector nds2, Rcpp::IntegerVector endp1, 
                                          Rcpp::IntegerVector endp2, Rcpp::NumericVector wgs )
{
  std::vector<std::string> nodes_1, nodes_2;
  std::vector<int> endpoints_1, endpoints_2;
  std::vector<double> weights;

  for(int i=0; i<nds1.size(); i++)
    nodes_1.push_back(Rcpp::as<std::string>(nds1[i]));

  for(int i=0; i<nds2.size(); i++)
    nodes_2.push_back(Rcpp::as<std::string>(nds2[i]));

  for(int i=0; i<endp1.size(); i++)
    endpoints_1.push_back(endp1[i]);

  for(int i=0; i<endp2.size(); i++)
    endpoints_2.push_back(endp2[i]);

  for(int i=0; i<wgs.size(); i++)
    weights.push_back(wgs[i]);

  std::string weight_flag = Rcpp::as<std::string>(cmp_flags[0]), 
              similarity_flag = Rcpp::as<std::string>(cmp_flags[1]);

  if( weight_flag == std::string("UNWEIGHTED") && similarity_flag == std::string("STANDARD") )
    return construct_community_hierarchy_internal<Standard_unweighted_multipartite_network_traits>
    ( nodes_1, nodes_2, endpoints_1, endpoints_2, weights);
  else if( weight_flag == std::string("UNWEIGHTED") && similarity_flag == std::string("MIN_WEIGHT_SUM") )
    return construct_community_hierarchy_internal<MSW_unweighted_multipartite_network_traits>
    ( nodes_1, nodes_2, endpoints_1, endpoints_2, weights);
  else if( weight_flag == std::string("WEIGHTED") && similarity_flag == std::string("STANDARD") )
    return construct_community_hierarchy_internal<Standard_weighted_multipartite_network_traits>
    ( nodes_1, nodes_2, endpoints_1, endpoints_2, weights);
  else if( weight_flag == std::string("WEIGHTED") && similarity_flag == std::string("MIN_WEIGHT_SUM") )
    return construct_community_hierarchy_internal<MSW_weighted_multipartite_network_traits>
    ( nodes_1, nodes_2, endpoints_1, endpoints_2, weights);
  else
  {
    std::cout << " Wrong method specified ... aborting " << std::endl;
    std::exit(-1);
  }

} // void cluster_network(...)


