/**
 *
 */
package de.helmholtz_muenchen.ibis.utils.ngs.frost;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Random;



/**
 * @author tanzeem.haque
 *
 */
public class InputScanner {

	/**
	 * @author tanzeem.haque
	 * This nested class is just an object for (int, boolean) type
	 * int: saving the position for mutation/recombination
	 * boolean: if true recombine, if false unaltered
	 * in case of de novos and mutations, the boolean will always be false
	 * this step was necessary because of the random seed =>
	 * although the positions for recombination were constant, 
	 * the recombination event itself was not constant
	 */
	private static class Position_Action {
		private int i;
		private boolean b;
		
		@SuppressWarnings("unused")
		/**
		 * Not used yet, but soon!
		 * @param i
		 * @param b
		 */
		public Position_Action(int i , boolean b) {
			// TODO Auto-generated constructor stub
			setI(i);
			setB(b);
		}

		public int getI() {
			return i;
		}

		public void setI(int i) {
			this.i = i;
		}

		public boolean isB() {
			return b;
		}

		public void setB(boolean b) {
			this.b = b;
		}

	}

	/**
	 * Attributes
	 */
	private double mutRate;
	private int recombination;
	private int seed;
	private int chunk;


	private InputData iData_parents;
	private InputData iData_recombination_child;
	private InputData iData_deNovo_child;


	public InputScanner(double mutRate, int recombination, int seed, int chunk) {
		this.mutRate = mutRate;
		this.recombination = recombination;
		this.seed = seed;
		this.chunk = chunk;
		this.iData_parents = new InputData("", new ArrayList<ArrayList<Integer>>());
		this.iData_deNovo_child = new InputData("", new ArrayList<ArrayList<Integer>>(10));
		this.iData_recombination_child = new InputData("", new ArrayList<ArrayList<Integer>>(this.recombination));;

	}

/*
	protected double getMutRate() {
		return mutRate;
	}

	protected void setMutRate(double mutRate) {
		this.mutRate = mutRate;
	}

	protected int getRecombination() {
		return recombination;
	}

	protected void setRecombination(int recombination) {
		this.recombination = recombination;
	}

	protected int getSeed() {
		return seed;
	}

	protected void setSeed(int seed) {
		this.seed = seed;
	}
*/
	
	public int getChunk() {
		return this.chunk;
	}

	/*
	public void setChunk(int chunk) {
		this.chunk = chunk;
	}
	*/

	/**
	 * @return the iData_arrList_parents
	 */
	protected InputData getiData_parents() {
		return this.iData_parents;
	}

	/**
	 * @return the iData_arrList_child
	 */
	protected InputData getiData_recombination_child() {
		return this.iData_recombination_child;
	}


	/**
	 * @return the iData_arrList_deNovo_child
	 */
	protected InputData getiData_deNovo_child() {
		return iData_deNovo_child;
	}


	@SuppressWarnings("unused")
	private int getMutationRange() {
		// TODO Auto-generated method stub
		//7994+1 for all chr
		int range = (int)((Math.pow(10, 8)/5300)/this.mutRate)+1;
		System.out.println("Range limit for each variant: " + range);
		return range;
	}

	private int getVariants(int length) {
		// TODO Auto-generated method stub
		//5878.76+1=5879: chr21
		//31019.84+1=31020: chr1
		double a = length/Math.pow(10, 8);// * getMutRate());
		double b = a * this.mutRate*5300;
		int variants =  (int)(b/1);
//		int variants = 10;
		System.out.println("Variants: " + variants);
		return variants;

	}

	private int getDeNovo(int length) {
		// TODO Auto-generated method stub
		//5878/5300=1: chr21
		//31020/5300=5.85=5: chr1
//		System.out.println("Length of chr: " + length);
//		System.out.println("Number of variants: " + getVariants(length));
//		System.out.println("De Novo: " + (int)(getVariants(length)/5300));
		
		int denovo = (int)(getVariants(length)/5300);
//		int denovo = 5;
		System.out.println("denovo: " + denovo);
		return denovo;
	}

	/**
	 * @param fs
	 * @param variants
	 * @param recombination
	 * @throws InterruptedException
	 */
	protected void prepare(String newId, int chr_length) throws InterruptedException {
		
		ArrayList<ArrayList<Integer>> parents;
		ArrayList<ArrayList<ArrayList<Integer>>> child;

		//Prepare parental genotypic data before the mutation
		//for each sequence read by fastareade
		//fs.getLength(i) returns the length of each sequence
		//from this length we will prepare random positions.
		//the number of positions is provided by variants
		parents = prepare_parental_input(chr_length);
		this.iData_parents = new InputData(newId, parents);
		child = prepare_child_input(chr_length);
		this.iData_deNovo_child = new InputData(newId, child.get(0));
		this.iData_recombination_child = new InputData(newId, child.get(1));
		
	}

	/**
	 * @param range: Length of the sequence
	 * @param recombination << variants
	 * @return the sorted positions for each sequence and for each parent to mutate
	 * @throws InterruptedException
	 * Also creates the position for denovo mutation
	 */
	private ArrayList<ArrayList<Integer>> prepare_parental_input(int length) throws InterruptedException {
		// TODO Auto-generated method stub

		int max = length, min = 1, n = getVariants(length);
		/**
		 * SEED ENTRY 1
		 */
		Random rd = new Random(); //new Random();
		ArrayList<ArrayList<Integer>> m_List = new ArrayList<ArrayList<Integer>>(this.chunk);
		ArrayList<Integer> m_tmpList = generatePosition(max, min, n, rd, false); //false means no recombination
		Collections.sort(m_tmpList);
		
		m_List = splitList(m_tmpList);
		for (int i = 0; i < m_List.size(); i++){
			for (int j = 0; j < m_List.get(i).size(); j++) {
//				System.out.println("M: " + m_List.size() + "\t" + m_List.get(i).size() + "\t" + i  + "\t" + j + "\t" + m_List.get(i).get(j)); 	
			}
		}
		return m_List;
	}

	
	private ArrayList<ArrayList<ArrayList<Integer>>> prepare_child_input(int length) throws InterruptedException {
		// TODO Auto-generated method stub
		// n must be < recombination
		ArrayList<ArrayList<ArrayList<Integer>>> child = new ArrayList<>(2);
		ArrayList<ArrayList<Integer>> r_List = new ArrayList<>(this.chunk);
		ArrayList<ArrayList<Integer>> d_List = new ArrayList<>(this.chunk);

		int n = getDeNovo(length);
		int m = this.recombination;
		ArrayList<Integer> r_tmpList = new ArrayList<Integer>(m);
		/**
		 * SEED ENTRY 2
		 */
		Random rd1 = new Random(); //new Random();
		r_tmpList = generatePosition(length, 1, m, rd1, true); // true coz recombination

		ArrayList<Integer> d_tmpList = new ArrayList<Integer>(n);
		/**
		 * SEED ENTRY 3
		 */
		Random rd2 = new Random(); //new Random(getSeed());
		d_tmpList = generatePosition(length, 1, n, rd2, false); // false coz no recombination

		Collections.sort(d_tmpList);
		d_List = splitList(d_tmpList);
		child.add(d_List);
		Collections.sort(r_tmpList);
		r_List = splitList(r_tmpList);
		child.add(r_List);

		for (int i = 0; i < d_List.size(); i++) {
			for (int j = 0; j < d_List.get(i).size(); j++) {
//				System.out.println("D :" + d_List.size() + "\t" + d_List.get(i).size() + "\t" + i + "\t" + j + "\t" + d_List.get(i).get(j)); 	
			}
		}
		for (int i = 0; i < r_List.size(); i++) {
			for (int j = 0; j < r_List.get(i).size(); j++) {
//				System.out.println("R :" + r_List.size() + "\t" + r_List.get(i).size() + "\t" + i + "\t" + j + "\t" + r_List.get(i).get(j)); 	
			}
		}
		
		return child;

	}

	/**
	 * @param arrList1 the base arraylist with the ordered position
	 * @param i = # variants
	 */
	private ArrayList<Integer> generatePosition(int max, int min, int n, Random rand, boolean rec) {
		// TODO Auto-generated method stub
		ArrayList<Integer> shuffled;
		HashSet<Integer> hs = new HashSet<Integer>();
		//genereate i random positions
		for (int i = 0; i < n; i++) {
//			if (parent){
////				System.out.println("max: " + ((i*max)+max) + " " + "min: " + ((i*max)+min));
//
//				generatePosition_help((i*max)+max, (i*max)+min, hs);
//			}
//			else
				generatePosition_help(max, min, hs, rand);
		}

		shuffled = new ArrayList<Integer>(hs);

		return shuffled;
	}

	/**
	 * @param range
	 * @param hs
	 * @param rd
	 */
	protected void generatePosition_help(int max, int min, HashSet<Integer> hs, Random rand) {
//		Random rand;
//		if (getSeed() != 0) {
//		rand = new Random(getSeed());
//			System.out.println("Seed: " + getSeed());
//		}
//		else
//			rand = new Random();
		boolean unique = true;
		int n = rand.nextInt((max+1)-min)+min;
		boolean b = hs.add(n);
		unique = b;
		while (!unique) {
			n = rand.nextInt((max+1)-min)+min;
			b = hs.add(n);
			unique = b;
		}
	}

	private String prettyPrint(final InputData inData) {
		String s = "";
	
		for (int i = 0; i < inData.getPositions().size(); i++) {
			s += inData.getId() + "\t" + inData.getPositions().get(i) + "\n";
		}
		
		return s;
	}

	private ArrayList<ArrayList<Integer>> splitList(ArrayList<Integer> p_tmpList) {
		// TODO Auto-generated method stub

		ArrayList<ArrayList<Integer>> list = new ArrayList<ArrayList<Integer>>(this.chunk);
		ArrayList<Integer> subList = new ArrayList<Integer>();

		int list_idx = 0;//p_tmpList_idx = 0,
		
		for (int j = 0; j < p_tmpList.size(); j++) {

			if (p_tmpList.get(j) <= (list_idx + 1) * FrostRunner.chunk_length){
//				System.out.println(p_tmpList.get(j) + "\t" + (list_idx + 1) * FrostRunner.chunk_length);
				subList.add(p_tmpList.get(j));
//				System.out.println(subList.size());
//				System.out.println(p_tmpList.get(j) + "\t" + list_idx);
			}
			else {
				list.add(subList);
				j--;				
				list_idx++;
				subList = new ArrayList<>();
//				System.out.println("ELSE " + list_idx + "\t" + j);
				continue;
			}
			/**
			 * the empty chunks: case 1
			 */
//			int x = 1;
//			while (p_tmpList.get(j) > (list_idx + x) * FrostRunner.chunk_length) {
//				System.out.println("first while " + p_tmpList.get(j) + "\t" + (list_idx + x) * FrostRunner.chunk_length);
//				subList = new ArrayList<>();
//				list.add(subList);
//				x++;
//				continue;
//			}
			if (j == p_tmpList.size()-1) {
				
				list.add(subList);
				list_idx++;
				subList = new ArrayList<>();
//				int x = 0;
				while(list_idx  < this.chunk ){
					if ((list_idx) * FrostRunner.chunk_length > p_tmpList.get(j)) {
//						System.out.println("second while " + p_tmpList.get(j) + "\t" + (list_idx) * FrostRunner.chunk_length);
						subList = new ArrayList<>();
						list.add(subList);
						list_idx++;
						continue;
					}
			}
//				System.out.println("ELSE " + list_idx + "\t" + j);
				break;
			}
		}	
		
		for (int n = 0; n < list.size(); n++)
			for(int m = 0; m < list.get(n).size(); m++) {
//				System.out.println(n + "\t" + list.size() + "\t" + m + "\t" + list.get(n).size() + "\t" + list.get(n).get(m));
			}
		return list;
	}
	
	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		// TODO Auto-generated method stub
//		return "Parental Mutation" + "\n" + prettyPrint(iData_arrList_parents) + "\n";
		return	"de Novo Mutation" + "\n" + prettyPrint(iData_deNovo_child) + "\n";
//				+ "Recombination" + "\n" +prettyPrint(iData_arrList_recombination_child);

	}

}