/**
 *
 */
package de.helmholtz_muenchen.ibis.utils.ngs.frost;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Random;
import java.util.Set;



/**
 * @author tanzeem.haque
 *
 */
public class InputScanner {

	/**
	 * @author tanzeem.haque
	 * 
	 */


	/**
	 * Attributes
	 */
	private String chrID;
	private double mutRate;
	private int seed;
	private boolean mutVary;
	private boolean recVary;
	private boolean deNovoVary;
//	private int chunk;
	
	
	public int inDelCounter_F = 0;
	public int inDelCounter_M = 0;
	public int variants = 0;
	public int deNovo = 0;
	/**
	 * Extra global variable for checking the chr length to correct the reco position
	 */
	private int chrLength;


	private InputData iData_parents;
	private InputData iData_recombination_child;
	private InputData iData_deNovo_child;
	private ArrayList<Strand> nMap;
//	private ArrayList<String> strandMap = new ArrayList<>();



	public InputScanner(String chrID, double mutRate, int seed, 
			boolean mutVary, boolean recVary, boolean deNovoVary, /*int chunk,*/ GenomeMap gm) {
		this.chrID = chrID;
		this.mutRate = mutRate;
		this.seed = seed;
		this.mutVary = mutVary;
		this.recVary = recVary;
		this.deNovoVary = deNovoVary;
		this.iData_parents = new InputData("", new ArrayList<Integer>());
		this.iData_deNovo_child = new InputData("", new ArrayList<Integer>(10));
		this.iData_recombination_child = new InputData("", new ArrayList<Integer>());;
		/**
		 * get the N region map at first
		 */
		this.nMap = gm.getGenomeMap().get(this.chrID);
//		for (int i = 0; i < this.nMap.size(); i++)
//			System.out.println(this.chrID + "\t" + this.nMap.get(i).get(0) + "\t" + this.nMap.get(i).get(1));
	
	}

	public ArrayList<Strand> getStrandMap() {
		return this.nMap;
	}

	private Random [] getSeed() {
		/**
		 * 0: mut seed
		 * 1: rec seed
		 * 2: denovo seed
		 */
		Random[] seeds4all = new Random [3];
		if(mutVary) {
			seeds4all[0] = new Random();
			seeds4all[1] = new Random(this.seed);
			seeds4all[2] = new Random(this.seed+1);
		}
		else if(recVary) {
			seeds4all[0] = new Random(this.seed);
			seeds4all[1] = new Random();
			seeds4all[2] = new Random(this.seed+1);
		}
		else if(deNovoVary) {
			seeds4all[0] = new Random(this.seed);
			seeds4all[1] = new Random(this.seed+1);
			seeds4all[2] = new Random();
		}
		return seeds4all;
	}
	

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
//		double a = length/Math.pow(10, 8);// * getMutRate());
//		double b = a * this.mutRate*this.generation;
//		int variants = (int)(b/1);
//		variants = 240650;
		int variants = this.chrLength/200;// (this.chrLength is calculated at prepare_parental_input method)
//		System.out.println("Variants: " + variants);
		this.variants = variants;
		return variants;

	}

	private int getDeNovo(int length) {
		// TODO Auto-generated method stub
		//5878/5300=1: chr21
		//31020/5300=5.85=5: chr1	
//		int deNovo = (int)(getVariants(length)/this.generation);

		double a = length/Math.pow(10, 8);// * getMutRate());
		double b = a * this.mutRate*1; //this.generation;
		int deNovo = (int)(b/1);
		deNovo = (deNovo == 0)? 1: deNovo;
//		deNovo = 5;
//		System.out.println("DeNovo: " + deNovo);
		this.deNovo = deNovo;
		return deNovo;
	}

	/**
	 * @param fs
	 * @param variants
	 * @param recombination
	 * @throws InterruptedException
	 */
	protected void prepare(int chr_length) throws InterruptedException {
		
		this.chrLength = chr_length;
		ArrayList<Integer> parents;
		ArrayList<ArrayList<Integer>> child;

		//Prepare parental genotypic data before the mutation
		//for each sequence read by fastareade
		//fs.getLength(i) returns the length of each sequence
		//from this length we will prepare random positions.
		//the number of positions is provided by variants
		parents = prepare_parental_input(chr_length);
		this.iData_parents = new InputData(this.chrID, parents);
		child = prepare_child_input(chr_length);
		this.iData_deNovo_child = new InputData(this.chrID, child.get(0));
		this.iData_recombination_child = new InputData(this.chrID, child.get(1));
//		prepareStrandFile(parents, child.get(0));
		
	}

//	private void prepareStrandFile(ArrayList<Integer> parents,
//			ArrayList<Integer> deNovoChild){
//		// TODO Auto-generated method stub
//		ArrayList<Integer> allMutations = new ArrayList<>(parents.size()+deNovoChild.size());
//		for (int i : parents)
//			allMutations.add(i);
//		for (int i : deNovoChild)
//			allMutations.add(i);
//		
//		Collections.sort(allMutations);
//		
//		for (int n : allMutations) {
//			for(int i = 0; i < this.nMap.size(); i++) {
//				if (n >= this.nMap.get(i).getStart() && n <= this.nMap.get(i).getEnd()){
//					String strand = (n+1) + "\t" + this.nMap.get(i).getStream();
//					this.strandMap.add(strand);
//					break;
//				}			
//			}
//		}
//	}

	/**
	 * @param range: Length of the sequence
	 * @param recombination << variants
	 * @return the sorted positions for each sequence and for each parent to mutate
	 * @throws InterruptedException
	 * Also creates the position for denovo mutation
	 */
	private ArrayList<Integer> prepare_parental_input(int length) throws InterruptedException {
		// TODO Auto-generated method stub
		
		//int max = length-1, min = 0, n = getVariants(length);

		int max = 0, min = 0, n = 0;
		ArrayList<Integer> m_tmpList = new ArrayList<>();
		
		for (int i = 0; i < this.nMap.size(); i++) {
			Random rand = new Random();
			min = this.nMap.get(i).getStart();
			max = this.nMap.get(i).getEnd();
			/**
			 * one can vary this n!!! just by changing 100 and 300
			 */
			int range = rand.nextInt((300-100)+1)+100;
			n = (max - min)/range; //one SNP every 100-300th base
//			System.out.println("#Mut: " + range);
			if (n == 0)
				continue;
			Random rd = getSeed()[0];
			ArrayList<Integer> m_tmp = generatePosition(max, min, n, rd, false);
			for(int tmp_i : m_tmp)
				m_tmpList.add(tmp_i);
		}
	
		/**
		 * If mutVary == false 
		 * => mutation will not vary
		 * => mutation is fixed
		 * => use of seed for the fixed positions
		 */
		Set<Integer> m_tmpSet = new HashSet<Integer>(m_tmpList); 
		ArrayList<Integer> m_List = new ArrayList<Integer>(m_tmpSet);
		Collections.sort(m_List);
		FrostRunner.createLog(FrostRunner.bw_log, "Number of SNPs in the exome of parents: " + m_List.size());
//		System.out.println("Number of SNPs in the exome of parents: " + m_List.size());
//		m_List = splitList(m_tmpList);
//		for (int i = 0; i < m_tmpList.size(); i++){
//			System.out.println(m_tmpList.get(i)); 	
//		}
		return m_List;
	}

	
	private ArrayList<ArrayList<Integer>> prepare_child_input(int length) throws InterruptedException {
		// TODO Auto-generated method stub
		// n must be < recombination
		ArrayList<ArrayList<Integer>> child = new ArrayList<>(2);

		int denovo = getDeNovo(length);
		/**
		 * Recombination module
		 */
		int max = 0, min = 0, reco = 0;
		ArrayList<Integer> r_tmpList = new ArrayList<>();
		
		for(int i = 0; i < this.nMap.size(); i++) {
			min = (i == 0) ? 0 : this.nMap.get(i-1).getEnd();
			max = this.nMap.get(i).getStart();
			if (i == this.nMap.size()-1) {
				min = this.nMap.get(i).getEnd();
				max = length-1;
			}
			
			if (max - min < FrostRunner.reco_gap) {
//				System.out.println("Skipping reco");
				continue;
			}
			
			Random rand = new Random();
			reco = (max - min)/(rand.nextInt(((max-min)-FrostRunner.reco_gap)+1)+FrostRunner.reco_gap); //one SNP every 100-300th base
			if (reco == 0)
				continue;
			Random rd1 = getSeed()[1]; //new Random();	
			ArrayList<Integer> r_tmp = generatePosition(max, min, reco, rd1, false); // false coz no denovo
			
			for (int tmp_i : r_tmp)
				r_tmpList.add(tmp_i);
			
		}
		/**
		 * SEED ENTRY 2
		 * Random rd1 = getSeed()[1]; //new Random();
		 */
		
		/**
		 * If recVary == false 
		 * => crossovers will not vary
		 * => crossovers are fixed
		 * => use of seed for the fixed positions
		 */
		ArrayList<Integer> d_tmpList = new ArrayList<Integer>(denovo);
		/**
		 * SEED ENTRY 3
		 */
		Random rd2 = getSeed()[2]; //new Random(); //new Random(getSeed());
		/**
		 * If deNovoVary == false 
		 * => denovo positions will not vary
		 * => denovos are fixed
		 * => use of seed for the fixed positions
		 */	
		d_tmpList = generatePosition(length-1, 0, denovo, rd2, true); // false coz no recombination
		
		/**
		 * if the position for denovo is an N, it will be skipped and we might not get any denovo at all
		 * so we force the program to generate perfect denovos
		 */
//		System.out.println("Denovo size " +n + "\t" + d_tmpList.size());
		while (d_tmpList.size() < denovo) {
			d_tmpList = generatePosition(length, 1, denovo, rd2, false);
//			System.out.println("denovo N " + d_tmpList.size() );
		}
//
//		System.out.println("Denovos: " + "\t" + "actually " + denovo + "\t"+ "got " + d_tmpList.size());
//		for (int i = 0; i < d_tmpList.size(); i++)
//			System.out.println(d_tmpList.get(i));
		/**
		 * denovo pos
		 */
		Collections.sort(d_tmpList);
		child.add(d_tmpList);
		/**
		 * recombination pos
		 */
		Set<Integer> r_tmpSet = new HashSet<Integer>(r_tmpList); 
		ArrayList<Integer> r_List = new ArrayList<Integer>(r_tmpSet);
		Collections.sort(r_List);
		/**
		 * check whether the rec pos are at least 
		 * 18kbp away from each other or else correction
		 * 		checkRecoGap(r_tmpList);

		 */
		child.add(r_List);
		FrostRunner.createLog(FrostRunner.bw_log, "Number of crossover point in child: " + r_List.size());
		FrostRunner.createLog(FrostRunner.bw_log, "Number of total de novo: " + d_tmpList.size());
		
		return child;

	}

	/**
	 * 
	 * @param r_tmpList
	 * check whether the rec pos are at least 
	 * 18kbp away from each other or else correction
	 */
	@SuppressWarnings("unused")
	private void checkRecoGap(ArrayList<Integer> r_tmpList) {
		// TODO Auto-generated method stub
//		int j = 0;
		boolean cut = false;
		for (int i = 0; i < r_tmpList.size()-1; i++) {
//			System.out.println(i + ". " + r_tmpList.get(i));
			if(r_tmpList.get(i+1)-r_tmpList.get(i) < FrostRunner.reco_gap){
//				j++;
//				System.out.println("Before: " + (i+1) + ". " + r_tmpList.get(i+1) + "\t" + (r_tmpList.get(i+1)-r_tmpList.get(i)));
				if (r_tmpList.get(i)+FrostRunner.reco_gap <= this.chrLength) 
					r_tmpList.set(i+1, r_tmpList.get(i)+FrostRunner.reco_gap);
				else {
					r_tmpList.remove(i+1);
					cut = true;
				}
//				System.out.println("After: " + (i+1) + ". " + r_tmpList.get(i+1));
			}
		}
		if (cut) {
			FrostRunner.createLog(FrostRunner.bw_log,"Total corrected number of recombination: " + r_tmpList.size());
//			System.out.println("Total corrected number of recombination: " + r_tmpList.size());
		}
//		System.out.println(j);
		
	}

	/**
	 * @param arrList1 the base arraylist with the ordered position
	 * @param i = # variants
	 */
	private ArrayList<Integer> generatePosition(int max, int min, int n, Random rand, boolean denovo) {
		// TODO Auto-generated method stub
		ArrayList<Integer> shuffled;
		HashSet<Integer> hs = new HashSet<Integer>();
		//genereate i random positions

		for (int i = 0; i < n; i++) {
			generatePosition_help(max, min, hs, rand, denovo);
		}

		shuffled = new ArrayList<Integer>(hs);
//		System.out.println("shuffled: " + shuffled.size());
		return shuffled;
	}

	/**
	 * @param range
	 * @param hs
	 * @param rd
	 */
	protected void generatePosition_help(int max, int min, HashSet<Integer> hs, Random rand, boolean denovo) {
//		Random rand;
//		if (getSeed() != 0) {
//		rand = new Random(getSeed());
//			System.out.println("Seed: " + getSeed());
//		}
//		else
//			rand = new Random();
		int n;// = rand.nextInt((max+1)-min)+min;

		do {
			n = rand.nextInt((max-min)+1)+min;
			/**
			 * only if its deNovo then skipN check
			 * 
			 */
			if (denovo && skipN(n)) {
//				FrostRunner.skipped_N++;
				return;
			}
		} while (!hs.add(n));
		
	}

	/**
	 * 
	 * @param n only for denovo coz this number is quite random and cannot not fixed-> depends on mutation rate
	 * @return
	 */
	private boolean skipN(int n) {
		boolean a = false;		
		for(int i = 0; i < this.nMap.size(); i++) {
			int prevEnd = (i == 0) ? 0 : this.nMap.get(i-1).getEnd();
			if ((n >= prevEnd && n <= this.nMap.get(i).getStart()) || 
			(i == this.nMap.size()-1 && n > this.nMap.get(i).getEnd())){
//				System.out.println("skipping");
				a = true;
				break;
			}
		
			
		}
//		if (a) 
//			FrostRunner.skipped_N++;
		return a;
	}


	private String prettyPrint(final InputData inData) {
		String s = "";
	
		for (int i = 0; i < inData.getPositions().size(); i++) {
			s += inData.getId() + "\t" + inData.getPositions().get(i) + "\n";
		}
		
		return s;
	}

/*
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
*/	
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