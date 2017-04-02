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
		Integer phie_e = b - c;
		phie.put(e, phie_e);
	}

	public void calculate_psie(Edge e, int c, int b){
		Integer psie_e = b - 2*c;
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
		Integer maxPhie = -Integer.MAX_VALUE;
		Integer psieEd = 0; 
		Edge maxEdge =  new Edge();
		for (Integer v : adj.get(d)) {
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
		for (Integer u : V) {
			System.out.println("Adyacentes del nodo " + u + "\n");
			for (Integer v : adj.get(u)) {
				System.out.println(v);
			}
		}
	}


	public void heuristic(int d){
		Graph factibleCicle = new Graph();

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
			Integer maxphi = -Integer.MAX_VALUE;
			Edge maxedge = new Edge();
			for (Edge e : E) {
				if (e.u() == d){
					maxphi = Math.max(maxphi, phie.get(e));
					maxedge = e;
				}
			}
		}

		Integer b = d;

		Edge ebu = new Edge();
		Boolean ebu_exists = false;
		Integer enrl = 0;

		while (!T.isEmpty()){
			// Existe un lado cuya fuente sea b y su destino u? 
			for (Edge e : T) {
				for (Integer u : V) {
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
				factibleCicle.addEdge(ebu);
				b = enrl;
			} else {
				Set<LinkedHashSet<Integer>> CCM = new HashSet<LinkedHashSet<Integer>>();
				for (Edge e : T) {
					break;
				}
			}

		}
	}
}

