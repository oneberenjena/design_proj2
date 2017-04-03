import java.util.*;
import java.lang.Math;


public class Graph {

	private Set<Integer> V;
	private Set<Edge> E;
	private Map<Integer, List<Integer>> adj;
	private Map<Edge, Integer> cost;	
	private Map<Edge, Integer> benefit;
	private Map<Edge, Integer> phie;
	private Map<Edge, Integer> psie;
	private Set<Edge> R;
	private Set<Edge> P;
	private Set<Edge> Q;
	private List<Edge> factibleCicle;
	private List<Edge> sol_parcial;
	private List<Edge> mejor_sol;
	private int beneficio_disponible;
	private Map<Edge, Integer> suc_cost;
	private Map<Edge, Integer> suc_benefit;

	public Graph(){
		this.V = new HashSet<Integer>();
		this.E = new HashSet<Edge>();
		this.cost = new HashMap<Edge, Integer>();
		this.benefit = new HashMap<Edge, Integer>();
		this.adj = new HashMap<Integer, List<Integer>>();
		this.phie = new HashMap<Edge, Integer>();
		this.psie = new HashMap<Edge, Integer>();
		this.R = new HashSet<Edge>();
		this.P = new HashSet<Edge>();
		this.Q = new HashSet<Edge>();
		this.factibleCicle = new ArrayList<Edge>();
		this.sol_parcial = new ArrayList<Edge>();
		this.mejor_sol = new ArrayList<Edge>();
		this.beneficio_disponible = 0;
		this.suc_cost = new HashMap<Edge, Integer>();
		this.suc_benefit = new HashMap<Edge, Integer>();
	}

	public void addEdge(int u, int v, int c, int b){
		if (contains(u,v)){return;}

		if(!V.contains(u)) {V.add(u);}
		if(!V.contains(v)) {V.add(v);}
		
		Edge e = new Edge(u,v,c,b);
		// Edge ee = new Edge(v,u,);
		E.add(e);
		// E.add(ee);

		if (adj.get(u)==null) {
			List<Integer> adjs1 = new ArrayList<Integer>();
			adjs1.add(v);
			adj.put(u, adjs1);
			if (adj.get(v) == null){
				List<Integer> adjs2 = new ArrayList<Integer>();
				adjs2.add(u);
				adj.put(v, adjs2);
			} else {
				adj.get(v).add(u);
			}
		}else {
			adj.get(u).add(v);
			if (adj.get(v) == null){
				List<Integer> adjs2 = new ArrayList<Integer>();
				adjs2.add(u);
				adj.put(v, adjs2);
			} else {
				adj.get(v).add(u);
			}
		}

		cost.put(e, c);
		// cost.put(ee, c); //
		benefit.put(e, b);
		// benefit.put(ee, b); //

		calculate_phie(e, c, b);
		// calculate_phie(ee, c, b); //
		calculate_psie(e, c, b);
		// calculate_psie(ee, c, b); //
	}

	public void addEdge(Edge e){
		if (!E.contains(e)) {E.add(e);}
		if (!V.contains(e.u())) {V.add(e.u());}
		if (!V.contains(e.v())) {V.add(e.v());}
	}

	public void calculate_phie(Edge e, int c, int b){
		int phie_e = b - c;
		phie.put(e, phie_e);
	}

	public void calculate_psie(Edge e, int c, int b){
		int psie_e = b - 2*c;
		psie.put(e, psie_e);
	}

	public Integer get_phie(Edge e){
		return phie.get(e);
	}

	public Integer get_psie(Edge e){
		return psie.get(e);
	}

	public Boolean contains(int u, int v){
		Boolean containsUV = false;
		for (Edge e : this.E) {
			if (e.u() == u && e.v() == v) {containsUV = true;}
		}
		return containsUV;
	}

	public void calculate_P(){
		for (Edge e : this.E) {
			if (get_phie(e) < 0) {P.add(e);}
		}
	}

	public void calculate_R(){
		for (Edge e : this.E) {
			if (get_psie(e) >= 0) {R.add(e);}
		}
	}

	public void calculate_Q(){
		Set<Edge> EnotR = new HashSet<>();
		for (Edge e : E) {
			if (!R.contains(e)) {EnotR.add(e);}
		}
		for (Edge e : EnotR) {
			if (get_phie(e) >= 0) {Q.add(e);}
		}
	}

	public void print_graph(){
		for (Edge e : E) {
			System.out.println(e + " phie: " + phie.get(e) + " psie: " + psie.get(e) + "\n");
		}
	}

	public List<Edge> heuristic(int d){
		calculate_P();
		calculate_R();
		calculate_Q();

		Set<Edge> T = new HashSet<Edge>();
		T.addAll(R);
		T.addAll(Q);

		Boolean estaEnT = false;

		for (Edge e : T) {
			if (e.u() == d || e.v() == d){estaEnT = true;} 
		}

		if (!estaEnT){
			int maxphi = -Integer.MAX_VALUE;
			Edge maxedge = new Edge();
			for (Edge e : E) {
				if (e.u() == d){
					maxphi = Math.max(maxphi, phie.get(e));
					maxedge = e;
				}
			}
			factibleCicle.add(maxedge);
		}

		int b = d;

		Edge ebu = new Edge();
		Boolean ebu_exists = false;
		int enrl = 0;

		while (!T.isEmpty()){
			// Existe un lado cuya fuente sea b y su destino u? 
			for (Edge e : T) {
				for (int u : V) {
					if (e.u() == b && e.v() == u || e.u() == u && e.v() == b){
						ebu_exists = true;
						enrl = u;
					}
				}
			}

			// Verificacion de la busqueda anterior
			if (ebu_exists){
				ebu_exists = false;
				ebu = obtener_lado(T, b);
				if (!T.contains(ebu)) {
					int u = ebu.v();
					int v = ebu.u();
					Edge ebu_right = getEdgeFromE(u,v);
					T.remove(ebu_right);
				} else {
					T.remove(ebu);
				}
				factibleCicle.add(ebu);
				b = enrl;
			} else {
				List<ArrayList<Edge>> CCM = new ArrayList<ArrayList<Edge>>();
				for (Edge e : T) {
					for (int i : e.toSet()) {
						if (i == b) {continue;}
						ArrayList<Edge> CMib = minimumCostPath(i,b);  
						CCM.add(CMib);
					}
				}

				ArrayList<Edge> CMib = obtener_camino(CCM);  
				// Para guardar el i
				Edge ix = CMib.get(0);
				Edge xy = CMib.get(1);
				int first_u = ix.u();
				int first_v = ix.v();
				int sec_u = xy.u();
				int sec_v = xy.v();
				int i = 0;

				if (first_u != sec_u && first_u != sec_v){
					i = first_u;
				} else if (first_v != sec_u && first_v != sec_v){
					i = first_v;
				}
				// Anadir el camino al ciclo
				factibleCicle.addAll(CMib);
				// Borrar lados de CMib en T
				T = removePathEdgesFrom(T, CMib);
				b = i; 
			}
		}

		// Si el ultimo vertice del ciclo no es el deposito:
		int last_i = getLast_i();
		if (last_i != d){
			List<Edge> CMid = minimumCostPath(last_i,d);
			factibleCicle.addAll(CMid); 
		}

		return factibleCicle;
	}

	public int getLast_i(){
		Edge last_edge = factibleCicle.get(factibleCicle.size()-1);
		Edge last_lastEdge = factibleCicle.get(factibleCicle.size()-2);

		int last_u = last_edge.u();
		int last_v = last_edge.v();

		if (last_lastEdge.u() != last_u && last_lastEdge.v() != last_u){
			return last_u;
		}

		if (last_lastEdge.u() != last_v && last_lastEdge.v() != last_v){
			return last_v;
		}

		return 0;
	}

	public Set<Edge> removePathEdgesFrom(Set<Edge> set, List<Edge> CM){
		for (Edge e : CM) {
			if (set.contains(e)) {set.remove(e);}
			e = getEdgeFromE(e.v(), e.u());
			if (set.contains(e)) {set.remove(e);}
		}

		return set;
	}

	public Edge obtener_lado(Set<Edge> T, int b){
		int maxPhi = -Integer.MAX_VALUE;
		int prevMax = maxPhi;
		Edge beneficEdge = new Edge();
		for (Edge e : T) {
			if (e.u() == b || e.v() == b){
				prevMax =  maxPhi;
				maxPhi = Math.max(maxPhi, phie.get(e));
				if (maxPhi > prevMax){
					beneficEdge = e;
				}
			}	
		}

		if (beneficEdge.v() == b) {beneficEdge = beneficEdge.swap();}
		return beneficEdge;
	}

	public ArrayList<Edge> obtener_camino(List<ArrayList<Edge>> CCM){
		int maxPhi = -Integer.MAX_VALUE;
		int prevMax = 0;
		int phiPathPosition = 0;
		int maxPhiPathPosition = -1;
		int actPhi = 0;
		
		for (ArrayList<Edge> path : CCM) {
			for (Edge e : path) {
				if (!E.contains(e)){
					e = getEdgeFromE(e.v(), e.u());
				}
				actPhi += phie.get(e);
			}

			prevMax = maxPhi;
			maxPhi = Math.max(maxPhi, actPhi);

			if (maxPhi > prevMax){
				maxPhiPathPosition = phiPathPosition;
			}

			actPhi = 0;
			phiPathPosition++;
		}

		return CCM.get(maxPhiPathPosition);
	}

	public int getMinDistVertex(Set<Integer> dijkstra_Q, Map<Integer, Integer> dist){
		int min = Integer.MAX_VALUE;
		int prevmin = min;
		int tmp_minVertex = -1;

		for (int v : dijkstra_Q) {
			min = Math.min(min, dist.get(v));

			if (min < prevmin){
				tmp_minVertex = v;
				prevmin = min;
			}
		}

		return tmp_minVertex;
	}

	public Edge getEdgeFromE(int u, int v){
		for (Edge e : E) {
			if (e.u() == u && e.v() == v) {return e;}
			if (e.u() == v && e.v() == u) {return e;}
		}
		return null;
	}

	public int cost(Edge e){
		if (factibleCicle.contains(e)){
			return cost.get(e);
		} else if (P.contains(e)){
			return -phie.get(e);
		} else {
			return 0;
		}
	}

	// Algoritmo de dijkstra extraido del pseudocodigo de wikipedia
	// https://en.wikipedia.org/wiki/Dijkstra%27s_algorithm
	// Calcula el camino minimo desde el nodo source al target
	public ArrayList<Edge> minimumCostPath(int source, int target){
		Set<Integer> dijkstra_Q = new HashSet<Integer>();
		Map<Integer, Integer> dist = new HashMap<Integer, Integer>();
		Map<Integer, Integer> prev = new HashMap<Integer, Integer>();

		for (Integer v : V) {
		 	dist.put(v, Integer.MAX_VALUE);
		 	prev.put(v, -1);
		 	dijkstra_Q.add(v);
		} 

		dist.put(source, 0);

		while (!dijkstra_Q.isEmpty()){
			int u = getMinDistVertex(dijkstra_Q, dist);
			if (u == target) {break;}
			dijkstra_Q.remove(u);

			for (int v : adj.get(u)) {
				if (!dijkstra_Q.contains(v)) {continue;}
				Edge e = getEdgeFromE(u,v);
				int alt = dist.get(u) + cost(e);
					
				if (alt < dist.get(v)){
					dist.put(v,alt);
					prev.put(v,u);
				}
			}
		}
		List<Integer> S = new ArrayList<>();
		int u = target;

		while (prev.get(u) != -1){
			S.add(u);
			u = prev.get(u);
		}

		S.add(u);

		ArrayList<Edge> path = new ArrayList<>();
		for (int i = S.size()-1; i >= 0 ; i--) {
			int v = S.get(i);
			if (i-1 >= 0){
				int w = S.get(i-1);
				Edge e = getEdgeFromE(v,w);

				if (e.u() == w && e.v() == v){
					path.add(e.swap());
				} else {
					path.add(e);
				}
			}
		}

		return path;
	}

	/********************** SEGUNDA PARTE DEL PROYECTO ***********************/

	public List<Edge> branchAndBound(List<Edge> sol_inicial, int d){
		sol_parcial.add(d);
		mejor_sol = sol_inicial;
		beneficio_disponible = getMaxBenefit();

		DFS();
	}

	public int calcBenefit(Edge e){
		if (!E.contains(e)) {e = getEdgeFromE(e.v(), e.u());}
		int benefit_e = benefit.get(e);
		benefit_e -= cost.get(e);
		if (benefit_e == psie.get(e)){
			benefit.put(e, 0);
			benefit_e = -cost.get(e);
		} else {
			benefit.put(e, benefit_e);
		}
		return benefit_e;
	}

	public int getMaxBenefit(){
		int total_benefit = 0;
		for (Edge e : factibleCicle) {
			if (!E.contains(e)) {
				e = getEdgeFromE(e.v(), e.u()); 
			}
			
			total_benefit += calcBenefit(e); 
		}
		return total_benefit;
	}

	public List<Edge> getListSucessors(int v){
		List<Edge> L_v = new ArrayList<Edge>();

		for (Edge e : E) {
			if (e.u() == v || e.v() == v){
				Edge e1 = new Edge(e.u(), e.v(), e.benefit(), e.cost());
				Edge e2 = new Edge(e.u(), e.v(), 0, e.cost());
				L_v.add(e1);
				L_v.add(e2);
				for (Edge e_i : L_v) {
					for (Edge e_j : L_v) {
						if (e_i == e_j) {continue;}
						int i = L_v.indexOf(e_i);
						int j = L_v.indexOf(e_j);
						if (e_i.benefit() < e_j.benefit()){
							Edge etmp = e_i;
							L_v.set(i, e_j);
							L_v.set(j, etmp);
						}
					}
				}
			}
		}

		return L_v;
	}
}

