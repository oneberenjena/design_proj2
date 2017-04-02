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
	}

	public void addEdge(int u, int v, int c, int b){
		if (contains(u,v)){return;}

		if(!V.contains(u)) {V.add(u);}
		if(!V.contains(v)) {V.add(v);}
		
		Edge e = new Edge(u,v);
		Edge ee = new Edge(v,u);
		E.add(e);
		E.add(ee);

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
		benefit.put(e, b);

		calculate_phie(e, c, b);
		calculate_psie(e, c, b);
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
		Set<Edge> Ecopy = E;
		Ecopy.removeAll(R);
		for (Edge e : Ecopy) {
			if (get_phie(e) >= 0) {Q.add(e);}
		}
	}

	public Edge max_phie_neighboor(int d){
		int maxPhie = -Integer.MAX_VALUE;
		int psieEd = 0; 
		Edge maxEdge =  new Edge();
		for (int v : adj.get(d)) {
			Edge ed = new Edge(d,v);
			psieEd = get_psie(ed);
			if (psieEd > maxPhie) {
				maxPhie = psieEd;
				maxEdge = ed;
			}
		}
		return maxEdge;
	}

	public void print_graph(){
		for (Edge e : E) {
			System.out.println(e + " phie: " + phie.get(e) + " psie: " + psie.get(e) + "\n");
		}
	}

	public void print_adjacencies(){
		for (int u : V) {
			System.out.println("Adyacentes del nodo " + u + "\n");
			for (int v : adj.get(u)) {
				System.out.println(v);
			}
		}
	}


	public void heuristic(int d){
		List<Integer> factibleCicle = new ArrayList<Integer>();

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
		}

		int b = d;

		Edge ebu = new Edge();
		Boolean ebu_exists = false;
		int enrl = 0;

		while (!T.isEmpty()){
			// Existe un lado cuya fuente sea b y su destino u? 
			for (Edge e : T) {
				for (int u : V) {
					if (e.u() == b && e.v() == u){
						ebu_exists = true;
						enrl = u;
					}
				}
			}

			// Verificacion de la busqueda anterior
			if (ebu_exists){
				ebu = obtener_lado(T, b);
				T.remove(ebu);
				///////////////////////////
				System.out.println("Adding edge (" + ebu.u() + "," + ebu.v() + ") to cicle");
				///////////////////////////
				factibleCicle.add(ebu.u());
				factibleCicle.add(ebu.v());
				b = enrl;
			} else {
				List<ArrayList<Integer>> CCM = new ArrayList<ArrayList<Integer>>();
				for (Edge e : T) {
					for (int i : e.toSet()) {
						ArrayList<Integer> CMib = minimumCostPath(i,b);  
						CCM.add(CMib);
					}
				}

				List<Integer> CMib = obtener_camino(CCM);  
				int i = CMib.get(0); // Esta es la unica parte que no se
				factibleCicle.addAll(CMib);
				// Borrar lados de CMib en T
				T = removePathEdgesFrom(T, CMib);
				b = i; // Entonces no estoy claro de esto
			}
		}

		// Si el ultimo vertice del ciclo no es el deposito:
		int last_i = factibleCicle.get(factibleCicle.size()-1);
		if (last_i != d){
			List<Integer> CMid = minimumCostPath(last_i,d);
			factibleCicle.addAll(CMid); 
		}
	}

	public Set<Edge> removePathEdgesFrom(Set<Edge> set, List<Integer> CM){
		for (int i = 0; i<CM.size(); i++) {
			int u = CM.get(0);
			if (CM.get(i+1) != null){
				int v = CM.get(i+1);
				Edge e = new Edge(u,v);
				if (set.contains(e)) {set.remove(e);}
			}
		}
		return set;
	}

	public Edge obtener_lado(Set<Edge> T, int b){
		int maxPhi = -Integer.MAX_VALUE;
		int prevMax = 0;
		Edge beneficEdge = new Edge();
		for (Edge e : T) {
			if (e.u() == b){
				prevMax =  maxPhi;
				maxPhi = Math.max(maxPhi, phie.get(e));
				if (maxPhi > prevMax){
					beneficEdge = e;
				}
			}	
		}

		return beneficEdge;
	}

	public ArrayList<Integer> obtener_camino(List<ArrayList<Integer>> CCM){
		int maxPhi = -Integer.MAX_VALUE;
		int prevMax = 0;
		int phiPathPosition = 0;
		int maxPhiPathPosition = -1;
		int actPhi = 0;
		
		for (ArrayList<Integer> path : CCM) {
			for (int v : path) {
				actPhi += phie.get(v);
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

	// Algoritmo de dijkstra extraido del pseudocodigo de wikipedia
	// https://en.wikipedia.org/wiki/Dijkstra%27s_algorithm
	// Calcula el camino minimo desde el nodo source al target
	public ArrayList<Integer> minimumCostPath(int source, int target){
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
			Q.remove(u);

			for (int v : adj.get(u)) {
				Edge e = new Edge(u,v);
				int alt = dist.get(u) + phie.get(e);

				if (alt < dist.get(v)){
					dist.put(v,alt);
					prev.put(v,u);
				}
			}
		}

		// ArrayList<Integer> S = new ArrayList<Integer>();
		Stack<Integer> S = new Stack<>();
		int u = target;

		while (prev.get(u) != -1){
			S.push(u);
			u = prev.get(u);
		}

		S.push(u);

		ArrayList<Integer> path = new ArrayList<>();
		for (int i = 0; i<S.size(); i++ ) {
			path.add(S.pop());
		}

		return path;
	}


}

