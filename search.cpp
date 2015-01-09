#include <iostream>
#include <sys/mman.h>


#include "ligra.h"
#include "graphfileshared.hpp"

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
void Compute(graph<vertex>& GA, intT start) {
  intT n = GA.n;
  //creates Parents array, initialized to all -1, except for start
  intT* Parents = newA(intT,GA.n);
  par::parallel_for(0l, GA.n, [&] (long i) {
    Parents[i] = -1;
  });
  Parents[start] = start;
  
  vertexSubset Frontier(n,start); //creates initial frontier
  
  while(!Frontier.isEmpty()){ //loop until frontier is empty
    vertexSubset output = edgeMap(GA, Frontier, BFS_F(Parents),GA.m/20);
    Frontier.del();
    Frontier = output; //set new frontier
  }
  Frontier.del();
  free(Parents);
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
  intE* ai = newA(intE, mm);
  intE* in_degrees = newA(intE, nn);
  for (intE i = 0; i < nn; i++)
    in_degrees = 0;
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
  lig = Ligra_graph(VV, nn, mm, ai);
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
  free(in_degrees);
  free(offsets);
  free(outps);
  free(inps);
  free(incounts);
}
  
template <class Adjlist, class Ligra_graph>
void search_benchmark_select_input_graph() {

}
  
}
}


/*---------------------------------------------------------------------*/



int main(int argc, char** argv) {

  using vtxid_type = intE;
  using adjlist_seq_type = pasl::graph::flat_adjlist_seq<vtxid_type>;
  using adjlist_type = pasl::graph::adjlist<adjlist_seq_type>;
  
  using ligra_type = graph<asymmetricVertex>;
  
  pasl::graph::search_benchmark_select_input_graph<adjlist_type, ligra_type>();
  
  using vtxid_type = typename adjlist_type::vtxid_type;
  adjlist_type graph;
  ligra_type lig;
  intT source;
  auto init = [&] {
    source = (intT) pasl::util::cmdline::parse_or_default_long("source", 0);
    pasl::util::cmdline::argmap_dispatch tmg;
    tmg.add("from_file",          [&] { load_graph_from_file(graph); });
    tmg.add("by_generator",       [&] { generate_graph(graph); });
    pasl::util::cmdline::dispatch_by_argmap(tmg, "load");
    convert(graph, lig);
    mlockall(0);
  };
  auto run = [&] (bool sequential) {
    Compute(lig, source);
  };
  auto output = [&] {
    //report(graph);
    print_adjlist_summary(graph);
  };
  auto destroy = [&] {
    
  };
  pasl::sched::launch(argc, argv, init, run, output, destroy);
  
  return 0;
}
