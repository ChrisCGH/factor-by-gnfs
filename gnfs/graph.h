#ifndef GRAPH_H
#define GRAPH_H
#include <map>
#include <set>
#include <vector>
#include <algorithm>
#include <functional>
#include <typeinfo>

template <class T>
struct Weight : public std::binary_function<T, T, double>
{
   double operator() (const T& t1, const T& t2) const
   {
      return weight(t1, t2);
   }
};

template <class T, class W = Weight<T> >
class Graph
{
   private:
      class Node
      {
         public:
            Node() : component_(0), next_in_component_(0), in_node_set_(false)
            {
               //std::cerr << typeid(*this).name() << "N+:" << node_id_ << std::endl;
            }
            Node(const T& data) : data_(data), component_(0), next_in_component_(0), in_node_set_(false)
            {
               //std::cerr << typeid(*this).name() << "N+:" << node_id_ << std::endl;
            }
            Node(const Node& n) : data_(n.data_), component_(n.component_), in_node_set_(n.in_node_set_)
            {
               //std::cerr << typeid(*this).name() << "N+:" << node_id_ << std::endl;
            }
            ~Node()
            {
               //std::cerr << typeid(*this).name() << "N-:" << node_id_ << std::endl;
            }
            int operator<(const Node& n) const
            {
               return (data_ < n.data_);
            }
            void set_component(size_t component) const
            {
               component_ = component;
            }
            bool is_in_nodeset() const
            {
                return in_node_set_;
            }

            T data_;
            mutable size_t component_;
            mutable Node* next_in_component_;
            mutable bool in_node_set_;
            //unsigned long int node_id_;
            friend std::ostream& operator<< (std::ostream& os, const Node& node)
            {
               os << "Node(" << node.data_ << "," << node.component_ << ")";
               return os;
            }
#if 0
            static unsigned long int get_next_node_id()
            {
               static unsigned long int next_node_id_ = 0;
               //std::cerr << typeid << "N+:" << next_node_id_ << std::endl;
               return next_node_id_++;
            }
#endif
      };

   public:
      Graph() : component_count_(0)
      {}
      void display()
      {
         //std::cout << "Nodes in graph: " << std::endl;
         for (node_set_iterator iter = node_set_.begin();
               iter != node_set_.end();
               ++iter)
         {
            std::cout << *(*iter) << std::endl;
         }

         //std::cout << "Connections: " << std::endl;
         for (graph_iterator iter = graph_.begin();
               iter != graph_.end();
               ++iter)
         {
            std::cout << *(iter->first) << "->" << *(iter->second) << std::endl;
         }
      }
      void connected_components()
      {
         size_t component = 0;
         for (graph_iterator iter = graph_.begin();
               iter != graph_.end();
               ++iter)
         {
            if (iter->first->component_ == 0)
            {
               ++component;
               set_component(iter->first, component);
            }
         }
         component_count_ = component;
         //node_set_.sort();
      }
      void set_component(Node* node, size_t component)
      {
         node->set_component(component);
         node_set_.insert(node);
         for (graph_iterator iter = graph_.lower_bound(node);
               iter != graph_.upper_bound(node);
               ++iter)
         {
            if (iter->second->component_ == 0)
            {
               set_component(iter->second, component);
            }
         }
      }
      void connect(const T& d1, const T& d2)
      {
         //std::cout << "Connect " << *d1 << " <-> " << *d2 << std::endl;
         Node* node1 = new Node(d1);
         Node* node2 = new Node(d2);
         graph_iterator found1 = graph_.find(node1);
         if (found1 != graph_.end())
         {
             delete node1;
             node1 = found1->first;
         }
         graph_iterator found2 = graph_.find(node2);
         if (found2 != graph_.end())
         {
            delete node2;
            node2 = found2->first;
         }
         graph_.insert(std::make_pair(node1, node2));
         graph_.insert(std::make_pair(node2, node1)); // to simulate a non-directed graph
      }

      int get_component_count() const
      {
         return component_count_;
      }
      void clear_components()
      {
         for (graph_iterator iter = graph_.begin();
               iter != graph_.end();
               ++iter)
         {
            iter->first.component_ = 0;
         }
         component_count_ = 0;
      }

      void minimum_spanning_tree(Graph<T>& tree)
      {
         // By Prim's algorithm (originally Jarnik)
         // Assumes that the Graph is connected, otherwise you
         // just get an MST of a connected sub-graph
         // Each connection must have a weight, which is given by
         // template parameter.
         // Make a heap of values (node, connection, weight)
         connected_components();
         if (get_component_count() != 1)
         {
            std::cerr << "Graph::minimum_spanning_tree() : graph is not connected" << std::endl;
            return;
         }
         typedef std::map<Node*, std::pair<Node*, double> > heap_type;
         typedef typename heap_type::iterator heap_type_iterator;
         heap_type node_heap;

         // Step 1.
         std::map<Node*, double> L;
         size_t i = 0;
         node_set_iterator it = node_set_.begin();
         Node* u = *it;
         L[u] = 0.0;
         ++it;

         for (; it != node_set_.end(); ++it)
         {
            L[*it] = 1.0e250;
            node_heap[*it] = std::make_pair<Node*, double>(0, 1.0e250);
         }

         bool done = false;
         while (!done)
         {
            // Step 2.
            double min_w = 1.0e250;
            heap_type_iterator min_it;
            for (heap_type_iterator it = node_heap.begin();
                  it != node_heap.end();
                  ++it)
            {
               Node* v = it->first;
               double w = W()(v->data_, u->data_);
               if (w < L[v])
               {
                  L[v] = w;
                  node_heap[v] = std::make_pair(u, w);
               }
               if (L[v] < min_w)
               {
                  min_w = L[v];
                  min_it = it;
               }
            }
            // Step 3.
            Node* next_u = min_it->first;
            // Step 4.
            node_heap.erase(min_it);
            tree.connect(u->data_, next_u->data_);
            u = next_u;
            ++i;
            if (i == node_set_.size() - 1) done = true;
         }

         tree.connected_components();
         //std::cerr << "tree has " << tree.get_component_count() << " components" << std::endl;
         //tree.display();
      }

   private:
   struct Node_less : public std::binary_function <Node*, Node*, bool>
      {
         bool operator() (const Node* _Left, const Node* _Right) const
         {
            return (*_Left < *_Right);
         }
      };
   struct Node_less_set : public std::binary_function <Node*, Node*, bool>
      {
         bool operator()(const Node* _Left, const Node* _Right) const
         {
            //std::cout << "Node_less_set::operator() : _Left=" << *_Left << " < _Right=" << *_Right << " -> ";
            if (_Left->component_ < _Right->component_)
            {
               //std::cout << "true" << std::endl;
               return true;
            }
            if (_Left->component_ == _Right->component_ && *_Left < *_Right)
            {
               //std::cout << "true" << std::endl;
               return true;
            }
            //std::cout << "false" << std::endl;
            return false;
         }
      };
      typedef std::multimap<Node*, Node*, Node_less> graph_type;
      typedef typename graph_type::iterator graph_iterator;
      graph_type graph_;
#if 0
      typedef std::set<Node*, Node_less_set> node_set_type;
      typedef typename node_set_type::iterator node_set_iterator;
#else
      class NodeSet
      {
         public:
            class component_iterator
            {
               public:
                  component_iterator() : curr_(0)
                  {}
                  component_iterator(Node* n) : curr_(n)
                  {}
                  component_iterator& operator++()
                  {
                     if (curr_)
                     {
                        curr_ = curr_->next_in_component_;
                     }
                     return *this;
                  }
                  Node* operator*()
                  {
                     return curr_;
                  }
                  bool operator!=(const component_iterator& rhs) const
                  {
                     return (curr_ != rhs.curr_);
                  }
                  bool operator==(const component_iterator& rhs) const
                  {
                     return (curr_ == rhs.curr_);
                  }

               private:
                  Node* curr_;
            };
            NodeSet()
            {
                heads_.reserve(4000000);
                nodes_.reserve(4000000);
            }
            ~NodeSet()
            {
               for (typename std::vector<Node*>::iterator it = nodes_.begin();
                    it != nodes_.end();
                    ++it)
               {
                   delete *it;
               }
            }
            component_iterator begin(size_t component)
            {
               if (component && component < heads_.size())
               {
                  Node* head = heads_[component].first;
                  return component_iterator(head);
               }
               return component_iterator();
            }
            component_iterator end(size_t component)
            {
               return component_iterator();
            }
            long int component_length(size_t component)
            {
               if (component && component < heads_.size())
               {
                  return heads_[component].second;
               }
               return 0;
            }
            void insert(Node* n)
            {
#if 0
               std::pair<typename std::set<Node*, Node_less_set>::iterator, bool> res = nodes_.insert(n);
               if (res.second)
               {
#else
               if (!n->is_in_nodeset())                
               { 
                  n->in_node_set_ = true;
                  nodes_.push_back(n); 
#endif
                  // a new Node in nodes_
                  size_t component = n->component_;
                  if (component >= heads_.size())
                  {
                     heads_.resize(component + 1);
                     heads_[component].first = 0;
                     heads_[component].second = 0;
                  }
                  if (heads_[component].first)
                  {
                     n->next_in_component_ = heads_[component].first;
                  }
                  heads_[component].first = n;
                  heads_[component].second++;
               }
            }
            size_t size()
            {
               return nodes_.size();
            }
            //typedef typename std::set<Node*, Node_less_set>::iterator iterator;
            typedef typename std::vector<Node*>::iterator iterator;
            iterator begin()
            {
               return nodes_.begin();
            }
            iterator end()
            {
               return nodes_.end();
            }
         private:
            std::vector<std::pair<Node*, long int> > heads_;
            //std::set<Node*, Node_less_set> nodes_;
            std::vector<Node*> nodes_;
            static bool greater ( const std::pair<Node*, long int>& elem1, const std::pair<Node*, long int>& elem2 )
            {
               return (elem1.second > elem2.second);
            }
         public:
            void sort()
            {
               std::sort(heads_.begin(), heads_.end(), greater);
            }

      };
      typedef NodeSet node_set_type;
      typedef typename node_set_type::iterator node_set_iterator;
      typedef typename node_set_type::component_iterator node_set_component_iterator;
#endif
      node_set_type node_set_;
      size_t component_count_;
   public:
      class connection_iterator
      {
         public:
            connection_iterator& operator++()
            {
               ++gi_;
               return *this;
            }
            std::pair<T, T> operator*()
            {
               return std::make_pair(gi_->first->data_, gi_->second->data_);
            }
            int operator==(const connection_iterator& ci) const
            {
               return (gi_ == ci.gi_);
            }
            int operator!=(const connection_iterator& ci) const
            {
               return !(*this == ci);
            }
            connection_iterator(graph_iterator gi) : gi_(gi)
            {}
         private:
            graph_iterator gi_;
      };
      class component_iterator
      {
         public:
            component_iterator& operator++()
            {
               ++curr_;
               return *this;
            }
            T operator*()
            {
               if (*curr_) return (*curr_)->data_;
               return 0;
            }
            bool operator==(const component_iterator& ci) const
            {
               return (curr_ == ci.curr_);
            }
            int operator!=(const component_iterator& ci) const
            {
               return (curr_ != ci.curr_);
            }
         public:
            component_iterator(size_t component, node_set_type& node_set) : component_(component), node_set_(node_set)
            {
               curr_ = node_set_.begin(component);
            }
            component_iterator(size_t component, node_set_type& node_set, bool the_end) : component_(component), node_set_(node_set)
            {
               curr_ = node_set_.end(component);
            }
            size_t component_;
            node_set_component_iterator curr_;
            node_set_type& node_set_;
      };
   public:
      component_iterator begin(size_t component)
      {
         component_iterator ci(component, node_set_);
         return ci;
      }
      component_iterator end(size_t component)
      {
         component_iterator ci(component, node_set_, true);
         return ci;
      }
      long int component_length(size_t component)
      {
         return node_set_.component_length(component);
      }
      connection_iterator begin()
      {
         connection_iterator ci(graph_.begin());
         return ci;
      }
      connection_iterator end()
      {
         connection_iterator ci(graph_.end());
         return ci;
      }
   public:
      ~Graph()
      {
         std::set<Node*> all_nodes;
         for (typename graph_type::iterator iter = graph_.begin();
               iter != graph_.end();
               ++iter)
         {
            //all_nodes.insert(iter->first);
            //delete iter->first;
            //iter->first = 0;
         }

#if 0
         for (typename std::set<Node*>::const_iterator iter = all_nodes.begin();
               iter != all_nodes.end();
               ++iter)
         {
            delete *iter;
         }
#endif
      }

};
#endif
