#include <iostream>
#include <sys/mman.h>


#include "ligra.h"
//#include "graphfileshared.hpp"

#include <fstream>
#include <iostream>
#include <istream>
#include <ostream>
#include <sstream>


intT nb_visited;

#include "edgelist.hpp"
#include "adjlist.hpp"

namespace pasl{
  namespace graph {

    template <class Adjlist>
void print_adjlist_summary(const Adjlist& graph) {
  std::cout << "nb_vertices\t" << graph.get_nb_vertices() << std::endl;
  std::cout << "nb_edges\t" << graph.nb_edges << std::endl;
}

    
static constexpr uint64_t GRAPH_TYPE_ADJLIST = 0xdeadbeef;
static constexpr uint64_t GRAPH_TYPE_EDGELIST = 0xba5eba11;

static const int bits_per_byte = 8;
static const int graph_file_header_sz = 5;

template <class Vertex_id>
void read_adjlist_from_file(std::string fname, adjlist<flat_adjlist_seq<Vertex_id>>& graph) {
  using vtxid_type = Vertex_id;
  std::ifstream in(fname, std::ifstream::binary);
  uint64_t graph_type;
  int nbbits;
  vtxid_type nb_vertices;
  edgeid_type nb_edges;
  bool is_symmetric;
  uint64_t header[5];
  in.read((char*)header, sizeof(header));
  graph_type = header[0];
  nbbits = int(header[1]);
  nb_vertices = vtxid_type(header[2]);
  nb_edges = edgeid_type(header[3]);
  is_symmetric = bool(header[4]);
  if (graph_type != GRAPH_TYPE_ADJLIST)
    util::atomic::die("read_adjlist_from_file"); //! \todo print the values expected and the value read
  if (sizeof(vtxid_type) * 8 < nbbits)
    util::atomic::die("read_adjlist_from_file: incompatible graph file");
  edgeid_type contents_szb;
  char* bytes;
  in.seekg (0, in.end);
  contents_szb = edgeid_type(in.tellg()) - sizeof(header);
  in.seekg (sizeof(header), in.beg);
  bytes = data::mynew_array<char>(contents_szb);
  if (bytes == NULL)
    util::atomic::die("failed to allocate space for graph");
  in.read (bytes, contents_szb);
  in.close();
  vtxid_type nb_offsets = nb_vertices + 1;
  if (contents_szb != sizeof(vtxid_type) * (nb_offsets + nb_edges))
    util::atomic::die("bogus file");
  graph.adjlists.init(bytes, nb_vertices, nb_edges);
  graph.nb_edges = nb_edges;
}

  } } //end namespace


struct BFS_F {
  intT* Parents;
  BFS_F(intT* _Parents) : Parents(_Parents) {}
  inline bool update (intT s, intT d) { //Update
    if(Parents[d] == -1) { Parents[d] = s; return 1; }
    else return 0;
  }
  inline bool updateAtomic (intT s, intT d){ //atomic version of Update
    return (CAS(&Parents[d],(intT)-1,s));
  }
  //cond function checks if vertex has been visited yet
  inline bool cond (intT d) { return (Parents[d] == -1); }
};

template <class vertex>
intT* ComputeBFS(graph<vertex>& GA, intT start) {
  intT n = GA.n;
  //creates Parents array, initialized to all -1, except for start
  intT* Parents = newA(intT,GA.n);
  par::parallel_for(0l, GA.n, [&] (long i) {
    Parents[i] = -1;
    }); 
  //  memset(Parents, -1, GA.n);
  Parents[start] = start;
  
  vertexSubset Frontier(n,start); //creates initial frontier
  
  while(!Frontier.isEmpty()){ //loop until frontier is empty
    vertexSubset output = edgeMap(GA, Frontier, BFS_F(Parents),GA.m/20);
    Frontier.del();
    Frontier = output; //set new frontier
  }
  Frontier.del();
  return Parents;
}

// **************************************************************
//    Non-DETERMINISTIC BREADTH FIRST SEARCH
// **************************************************************

// **************************************************************
//    THE NON-DETERMINISTIC BSF
//    Updates the graph so that it is the BFS tree (i.e. the neighbors
//      in the new graph are the children in the bfs tree)
// **************************************************************

template <class vertex>
intT* BFS_orig(graph<vertex> GA, intT start) {
  intT numVertices = GA.n;
  intT numEdges = GA.m;
  vertex *G = GA.V;
  intT* Frontier = newA(intT,numEdges);
  intT* Visited = newA(intT,numVertices);
  intT* FrontierNext = newA(intT,numEdges);
  intT* Counts = newA(intT,numVertices);
  par::parallel_for ((intT)0, (intT)numVertices, [&] (long i) {
  Visited[i] = 0;});
  //  memset(Visited, 0, numVertices); 

  Frontier[0] = start;
  intT frontierSize = 1;
  Visited[start] = 1;

  intT totalVisited = 0;
  int round = 0;

  while (frontierSize > 0) {
    round++;
    totalVisited += frontierSize;

    par::parallel_for ((intT)0, (intT)frontierSize, [&] (long i) {
        Counts[i] = G[Frontier[i]].getOutDegree();});
    intT nr = sequence::scan(Counts,Counts,frontierSize,addF<intT>(),(intT)0);

    // For each vertexB in the frontier try to "hook" unvisited neighbors.
    par::parallel_for ((intT)0, (intT)frontierSize, [&] (long i) {
      intT k= 0;
      intT v = Frontier[i];
      intT o = Counts[i];
      for (intT j=0; j < G[v].getOutDegree(); j++) {
        intT ngh = G[v].getOutNeighbor(j); //Neighbors[j];
	if (Visited[ngh] == 0 && CAS(&Visited[ngh],(intT)0,(intT)1)) {
	  FrontierNext[o+j] = ngh;
          k++;
	}
	else FrontierNext[o+j] = -1;
      }
      //      G[v].degree = k;
      });

    // Filter out the empty slots (marked with -1)
    frontierSize = sequence::filter(FrontierNext,Frontier,nr,nonNegF());
  }
  free(FrontierNext); free(Frontier); free(Counts); //free(Visited);
  //  return pair<intT,intT>(totalVisited,round);
  return Visited;
}


/*---------------------------------------------------------------------*/

namespace pasl {
namespace graph {
  
  bool should_disable_random_permutation_of_vertices;

  
template <class Adjlist, class Ligra_graph>
void convert(const Adjlist& adj, Ligra_graph& lig) {
  intT nb_vertices = adj.get_nb_vertices();
  intT nb_edges = adj.nb_edges;
  intT nn = nb_vertices;
  intT mm = nb_edges;
  using vertex_type = typename Ligra_graph::vertex_type;
  vertex_type* VV = newA(vertex_type, nn);
  intE* ai = newA(intE, 2 * mm);
  intE* in_degrees = newA(intE, nn);
  for (intE i = 0; i < nn; i++)
    in_degrees[i] = 0;
  for (intE i = 0; i < nb_vertices; i++) {
    intE n = adj.adjlists[i].get_out_degree();
    intE* neighbors = adj.adjlists[i].get_out_neighbors();
    for (intE j = 0; j < n; j++) {
      intE ngh = neighbors[j];
      in_degrees[ngh]++;
    }
  }
  intE* offsets = newA(intE, nn);
  for (intE i = 0; i < nn; i++) {
    offsets[i] = 0;
    offsets[i] += adj.adjlists[i].get_out_degree();
    offsets[i] += in_degrees[i];
  }
  intE acc = 0;
  for (intE i = 0; i < nn; i++) {
    intE newacc = acc + offsets[i];
    offsets[i] = acc;
    acc = newacc;
  }
  intE** outps = newA(intE*, nn);
  intE** inps = newA(intE*, nn);
  for (intE i = 0; i < nn; i++) {
    intE indeg = in_degrees[i];
    intE outdeg = adj.adjlists[i].get_out_degree();
    intE* outp = &ai[offsets[i]];
    intE* inp = outp + outdeg;
    outps[i] = outp;
    inps[i] = inp;
    VV[i].setOutNeighbors(outp);
    VV[i].setInNeighbors(inp);
    VV[i].setInDegree(indeg);
    VV[i].setOutDegree(outdeg);
    intE* outs = adj.adjlists[i].get_out_neighbors();
    for (intE k = 0; k < outdeg; k++) {
      outp[k] = outs[k];
    }
  }
  intE* incounts = newA(intE, nn);
  for (intE i = 0; i < nn; i++)
    incounts[i] = 0;
  for (intE i = 0; i < nn; i++) {
    intE n = adj.adjlists[i].get_out_degree();
    intE* neighbors = adj.adjlists[i].get_out_neighbors();
    for (intE j = 0; j < n; j++) {
      intE ngh = neighbors[j];
      intE c = incounts[ngh]++;
      intE* inp = inps[ngh];
      inp[c] = i;
    }
  }
  lig = Ligra_graph(VV, nn, mm, ai);
  free(in_degrees);
  free(offsets);
  free(outps);
  free(inps);
  free(incounts);
}


  
}
}


/*---------------------------------------------------------------------*/



int main(int argc, char** argv) {

  using vtxid_type = intE;
  using adjlist_seq_type = pasl::graph::flat_adjlist_seq<vtxid_type>;
  using adjlist_type = pasl::graph::adjlist<adjlist_seq_type>;
  
  using ligra_type = graph<asymmetricVertex>;
  
  using vtxid_type = typename adjlist_type::vtxid_type;
  ligra_type lig;
  intT source;
  intT* Parents;
  bool use_pbbs = false;
  auto init = [&] {
    pasl::graph::should_disable_random_permutation_of_vertices = pasl::util::cmdline::parse_or_default_bool("should_disable_random_permutation_of_vertices", false, false);
    adjlist_type graph;
    source = (intT) pasl::util::cmdline::parse_or_default_long("source", 0);
    std::string infile = pasl::util::cmdline::parse_or_default_string("infile", "");
    read_adjlist_from_file(infile, graph);
    /*
    pasl::util::cmdline::argmap_dispatch tmg;
    tmg.add("from_file",          [&] { load_graph_from_file(graph); });
    tmg.add("by_generator",       [&] { generate_graph(graph); }); 
    pasl::util::cmdline::dispatch_by_argmap(tmg, "load"); */
    convert(graph, lig);
    print_adjlist_summary(graph);
    mlockall(0);
    std::string algo = pasl::util::cmdline::parse_or_default_string("algo", "ligra");
    use_pbbs = algo != "ligra";
    if (! (algo == "ligra" || algo == "pbbs_pbfs_cilk")) {
      pasl::util::atomic::die("bogus -algo");
    }
  };
  auto run = [&] (bool sequential) {
    if (use_pbbs) {
      Parents = BFS_orig(lig, source);
    } else { 
      Parents = ComputeBFS(lig, source);
    }
  };
  auto output = [&] {
    nb_visited = 0;
    for (intT i = 0; i < lig.n; i++)
      if (Parents[i] >= 0)
        nb_visited++;
    free(Parents);
    std::cout << "nb_visited\t" << nb_visited << std::endl;
#ifdef COUNT_SHORTCUT
    std::cout << "nb_times_shortcut\t" << nbTimesShortcut << std::endl;
#endif    
  };
  auto destroy = [&] {
    
  };
  pasl::sched::launch(argc, argv, init, run, output, destroy);
  
  return 0;
}
