/**
 *
 */
package de.helmholtz_muenchen.ibis.utils.ngs.frost;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.Random;



/**
 * @author tanzeem.haque
 *
 */
public class InputScanner {

	/**
	 * Attributes
	 */
	private FastaReader fs;
	private double mutRate;
	private int recombination;
	private int seed;


	private ArrayList<InputData> iData_arrList_parents;
	private ArrayList<InputData> iData_arrList_recombination_child;
	private ArrayList<InputData> iData_arrList_deNovo_child;


	public InputScanner(FastaReader fs, double mutRate, int recombination, int seed) {
		setFs(fs);
		setMutRate(mutRate);
		setRecombination(recombination);
		setSeed(seed);
	}

	/**
	 * @return the fs
	 */
	protected FastaReader getFs() {
		return fs;
	}

	/**
	 * @param fs the fs to set
	 */
	protected void setFs(FastaReader fs) {
		this.fs = fs;
	}


	/**
	 * @return the deNovo_mutRate
	 */
	protected double getMutRate() {
		return mutRate;
	}

	/**
	 * @param mutRate the deNovo_mutRate to set
	 * 2.36e-8 default
	 */
	protected void setMutRate(double mutRate) {
		this.mutRate = mutRate;
	}

	/**
	 * @return the recombination
	 */
	protected int getRecombination() {
		return recombination;
	}

	/**
	 * @param recombination the recombination to set
	 */
	protected void setRecombination(int recombination) {
		this.recombination = recombination;
	}

	/**
	 * @return the seed
	 */
	protected int getSeed() {
		return seed;
	}

	/**
	 * @param seed the seed to set
	 */
	protected void setSeed(int seed) {
		this.seed = seed;
	}

	/**
	 * @return the iData_arrList_parents
	 */
	protected ArrayList<InputData> getiData_arrList_parents() {
		return this.iData_arrList_parents;
	}

	/**
	 * @return the iData_arrList_child
	 */
	protected ArrayList<InputData> getiData_arrList_recombination_child() {
		return this.iData_arrList_recombination_child;
	}


	/**
	 * @return the iData_arrList_deNovo_child
	 */
	protected ArrayList<InputData> getiData_arrList_deNovo_child() {
		return iData_arrList_deNovo_child;
	}

	/**
	 * @param fs
	 * @param variants
	 * @param recombination
	 * @throws InterruptedException
	 */
	protected void prepare() throws InterruptedException {
		iData_arrList_parents = new ArrayList<InputData>();
		//the number of de novo mutations for child is <= 10. really 7
		iData_arrList_deNovo_child = new ArrayList<InputData>(10);
		iData_arrList_recombination_child = new ArrayList<InputData>(getRecombination());

		ArrayList<Integer> mutation;
		ArrayList<ArrayList<Integer>> child;

		//Prepare parental genotypic data before the mutation
		//for each sequence read by fastareade
		//fs.getLength(i) returns the length of each sequence
		//from this length we will prepare random positions.
		//the number of positions is provided by variants

		for (int i = 0; i < getFs().size(); i++) {
			mutation = prepare_parental_input(getFs().getLength(i));
			this.iData_arrList_parents.add(new InputData(getFs().getIdentifier(i), mutation));
			child = prepare_child_input(getFs().getLength(i));
			this.iData_arrList_deNovo_child.add(new InputData(getFs().getIdentifier(i), child.get(0)));
			this.iData_arrList_recombination_child.add(new InputData(getFs().getIdentifier(i), child.get(1)));
		}
	}

	@SuppressWarnings("unused")
	private int getMutationRange() {
		// TODO Auto-generated method stub
		//7994+1 for all chr
		int range = (int)((Math.pow(10, 8)/5300)/getMutRate())+1;
		System.out.println("Range limit for each variant: " + range);
		return range;
	}

	private int getVariants(int length) {
		// TODO Auto-generated method stub
		//5878.76+1=5879: chr21
		//31019.84+1=31020: chr1
		double a = length/Math.pow(10, 8);// * getMutRate());
		double b = a * getMutRate() * 5300;
		int variants =  (int)(b/1);
		return variants;
	}

	private int getDeNovo(int length) {
		// TODO Auto-generated method stub
		//5878/5300=1: chr21
		//31020/5300=5.85=5: chr1
		System.out.println("Length of chr: " + length);
		System.out.println("Number of variants: " + getVariants(length));
		System.out.println("De Novo: " + (int)(getVariants(length)/5300));
		return (int)(getVariants(length)/5300);
	}

	/**
	 * @param range: Length of the sequence
	 * @param recombination << variants
	 * @return the sorted positions for each sequence and for each parent to mutate
	 * @throws InterruptedException
	 * Also creates the position for denovo mutation
	 */
	private ArrayList<Integer> prepare_parental_input(int length) throws InterruptedException {
		// TODO Auto-generated method stub

		int max = length, min = 1, n = getVariants(length);
//		ArrayList<Integer> p_arrList = generatePosition(max, min, n, length, true);
		Random rd = new Random(getSeed());
		ArrayList<Integer> p_arrList = generatePosition(max, min, n, rd);
		Collections.sort(p_arrList);

		return p_arrList;
	}


	private ArrayList<ArrayList<Integer>> prepare_child_input(int length) throws InterruptedException {
		// TODO Auto-generated method stub
		// n must be < recombination
		ArrayList<ArrayList<Integer>> child = new ArrayList<ArrayList<Integer>>(2);
		int n = getDeNovo(length);
		int m = getRecombination();
		ArrayList<Integer> pos_recombination = new ArrayList<Integer>(m);
//		pos_recombination = generatePosition(length, 1, m, -99, false);
		Random rd1 = new Random(getSeed());
		pos_recombination = generatePosition(length, 1, m, rd1);

		ArrayList<Integer> pos_deNovo = new ArrayList<Integer>(n);
//		pos_deNovo = generatePosition(length, 1, n, -99, false);
		Random rd2 = new Random(getSeed());
		pos_deNovo = generatePosition(length, 1, n, rd2);


//		Thread.sleep(1000);

		Collections.sort(pos_deNovo);
		child.add(pos_deNovo);
		Collections.sort(pos_recombination);
		child.add(pos_recombination);

		return child;

	}


	/**
	 * @param arrList1 the base arraylist with the ordered position
	 * @param i = # variants
	 */
	private ArrayList<Integer> generatePosition(int max, int min, int n, Random rand) {
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

	private String prettyPrint(final ArrayList<InputData> in_arrList) {
		String s = "";
		for (InputData inData : in_arrList) {
			for (int i = 0; i < inData.getPositions().size(); i++) {
				s += inData.getId() + "\t" + inData.getPositions().get(i) + "\n";
			}
		}
		return s;
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		// TODO Auto-generated method stub
//		return "Parental Mutation" + "\n" + prettyPrint(iData_arrList_parents) + "\n";
		return	"de Novo Mutation" + "\n" + prettyPrint(iData_arrList_deNovo_child) + "\n";
//				+ "Recombination" + "\n" +prettyPrint(iData_arrList_recombination_child);

	}

}