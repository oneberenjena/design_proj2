import java.util.Set;
import java.util.HashSet;

public class Edge {

	private int u;
	private int v;
	private int cost;
	private int benefit;

	public Edge(int u, int v, int cost, int benefit){
		this.u = u;
		this.v = v;
		this.cost = cost;
		this.benefit = benefit;
	}

	public Edge(){
		this.u = -1;
		this.v = -1;
	}

	public int u() {return u;}
	public int v() {return v;}
	public int cost() {return cost;}
	public int benefit() {return benefit;}

	public String toString() {
		StringBuilder s = new StringBuilder();
		s.append("(" + u + ", " + v + ")");
		return s.toString();
	}

	public Set<Integer> toSet(){
		Set<Integer> edge = new HashSet<Integer>();
		edge.add(u);
		edge.add(v);
		return edge;
	}	

	public Edge swap(){
		Edge e = new Edge(v,u);
		return e;
	}

}