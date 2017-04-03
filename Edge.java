import java.util.Set;
import java.util.HashSet;

public class Edge {

	private Integer u;
	private Integer v;

	public Edge(int u,int v){
		this.u = u;
		this.v = v;
	}

	public Edge(){
		this.u = -1;
		this.v = -1;
	}

	public Integer u() {return u;}
	public Integer v() {return v;}

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