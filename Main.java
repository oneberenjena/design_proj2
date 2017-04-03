import java.io.*;
import java.util.Scanner;
import java.util.List;
import java.util.ArrayList;

public class Main {
	public static void main(String[] args){	
		Graph G = buildGraph(args[0], Integer.parseInt(args[1]));
		// double start = System.currentTimeMillis();
		// String output_solution = G.runResolvePath(Integer.parseInt(args[1]));
		// try{
		// 	printFile(output_solution, args[0]);
		// } catch (IOException e){
		// 	e.printStackTrace();
		// }
	}

	// public static void printFile(String output, String file) throws IOException{
	// 	File outFile = new File(file + "-salida.txt");
	// 	FileWriter writer = new FileWriter(outFile);
	// 	writer.write(output);
	// 	writer.flush();
	// 	writer.close();
	// }

	public static Graph buildGraph(String arg, int deposit){
		File instance = new File(arg);
		try{
	        Scanner sc = new Scanner(instance);
	        // int V = sc.nextInt();
	        //int E = sc.nextInt();
	        Graph G = new Graph();
	        int w, v, p, c;
	        while(sc.hasNextInt()){
	        	v = sc.nextInt();
	        	w = sc.nextInt();
	        	c = sc.nextInt();
	        	p = sc.nextInt();
	        	G.addEdge(v, w, c, p);
	        }  
	        sc.close();

	        // G.print_graph();
	        List<Edge> solucionFactible = G.heuristic(deposit);
	        // System.out.println(solucionFactible);
	        return G;

		}catch (FileNotFoundException e) {
	        e.printStackTrace();
	        return null;
	    }
	}
}