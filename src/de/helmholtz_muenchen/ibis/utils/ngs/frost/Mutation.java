/**
 * 
 */
package de.helmholtz_muenchen.ibis.utils.ngs.frost;

import java.util.ArrayList;
import java.util.Random;





/**
 * @author tanzeem.haque
 *
 */
public class Mutation {
	/**
	 * Attribute
	 */
	private char[] bases;
	protected int M;
	protected int F;
	
	private String individuum = "";
	private String type = "";
	private ArrayList<ParentalGenotype> pg_arrList = new ArrayList<ParentalGenotype>();
	private ArrayList<Genotype> g_arrList = new ArrayList<Genotype>();//only for child
	private ArrayList<Allelic_Change> parent_ac_arrList = new ArrayList<Allelic_Change>();
	private ArrayList<Allelic_Change> child_ac_arrList = new ArrayList<Allelic_Change>();


	private final double IN_RATE = FrostRunner.insertion_rate;
	private final double DEL_RATE = FrostRunner.deletion_rate;

	public Mutation(char[] bases) {
		// TODO Auto-generated constructor stub
		setBases(bases);
	}

	public Mutation(char[] basesToMutate, int M, int F) {
		// TODO Auto-generated constructor stub
		setBases(basesToMutate);
		this.M = M;
		this.F = F;
	}
	/**
	 * @return the original base 
	 */
	protected char[] getBases() {
		return bases;
	}

	/**
	 * @param bases2
	 */
	protected void setBases(char[] bases) {
		this.bases = bases;
	}

	/**
	 * @return an arraylist of the 4 bases
	 */
	private static ArrayList<Character> base_list() {
		ArrayList<Character> arrList = new ArrayList<Character>();
		arrList.add('A');
		arrList.add('T');
		arrList.add('G');
		arrList.add('C');
		return arrList;
	}

	/**
	 * @return a random base from a list of bases
	 */
	private static char randomBase(ArrayList<Character> arrList) {
		Random base = new Random();
		return arrList.get(base.nextInt(arrList.size()));
	}

	/**
	 * @param toMutate
	 * @param f : true if father : GF ; for child false
	 * @param m : true if mother : GM ; for child false
	 * @return
	 */
	protected String[] mutationType(char[] toMutate, String s) {
		String[] allele = new String[2];
		Random rand = new Random();
		boolean insert , delete;
		//probability of insert 1/50
		insert = (rand.nextInt((int)(1/IN_RATE)) == 0) ? true : false; 
		// same probability for deletion
		delete = (rand.nextInt((int)(1/DEL_RATE)) == 0) ? true : false;

		//T && T = T, F && F = T
		if ((insert==true && delete==true)||(insert==false && delete==false))	{	
			allele = substitution(toMutate, s);
		}
		if (insert == true && delete == false) {
			allele = insertion(toMutate, s);
		}
		if (delete == true && insert == false) {
			allele = deletion(toMutate, s);
		}
		this.individuum = s;
		return allele;
	}

	/**
	 * @param m 
	 * @param f 
	 * @param all_ll an arraylist of size 2 with thwo chromatid sequences as linked liste
	 * @param i position for mutation
	 * @return chr_ar_ll => the changed ll arraylist
	 */
	private String[] insertion(char[] c, String s) { 
		this.type = "IN";

		String[] allele = new String [2];
		char[] chr = new char[2];
		c = getBases();

		// gets 2 random bases for two different chromatids
		char[] tmp = { randomBase(base_list()), randomBase(base_list()) }; 

		char chrd1_1 = tmp[0];
		char chrd2_1 = tmp[1];
		// 2nd possibility that one of the both gets insertion : second chromatid
		// 3rd possibility that one of the both gets insertion : first chromatid			
		char[][] chr_combi = { { c[0], chrd2_1 }, { chrd1_1, c[1] } };
		
//		else {
		Random rand = new Random();
		int x = rand.nextInt(2);
		chr[0] = chr_combi[x][0];
		chr[1] = chr_combi[x][1];
		//x=0: 0/1, x=1: 1/0
		if (x == 0) {
			allele[0] = chr[0]+"";
			allele[1] = c[1]+""+chr[1]+"";
			heterozygous_genotype(s, 0, c[1]+"", allele[1]);
		}
		else if (x == 1) {
			allele[0] = c[0]+""+chr[0];
			allele[1] = chr[1]+"";
			heterozygous_genotype(s, 1, c[0]+"", allele[0]);
		}
		return allele;
	}


	private String[] deletion(char[] c, String s) {
		// TODO Auto-generated method stub
		this.type = "DEL";

		String[] allele = new String [2];
		char[] chr = new char[2];
		c = getBases();
		
		char chrd1_1 = 42;//*
		char chrd2_1 = 42;//*

		char[][] chr_combi = { { c[0], chrd2_1 },{ chrd1_1, c[1] } };

		Random rand = new Random();
		int x = rand.nextInt(2); // either 1 or 2
		chr[0] = chr_combi[x][0];
		chr[1] = chr_combi[x][1];
		allele[0] = chr[0]+"";
		allele[1] = chr[1]+"";
		//x=0: 0/1, x=1: 1/0
		if (x == 0) 		
			heterozygous_genotype(s, 0, c[1]+"", chr[1]+"");	
		else if (x == 1) 
			heterozygous_genotype(s, 1, c[0]+"", chr[0]+"");	
		return allele;

	}

	private String[] substitution(char[] c, String s) {
		// TODO Auto-generated method stub
		String[] allele = new String [2];
		// 1st possibility one chromatid gets mutated
	
		char[] chr = new char[2];
		c = getBases();
		
		char snp = getSnp(c[0]);
		// The first chromatid gets mutated, the second one doesnt
		char chrd1_1 = snp;
		char chrd2_1 = c[1];
		// The second chromatid get mutated and the first one doesnt
		char chrd1_2 = c[0];
		char chrd2_2 = snp;

		// or set own probability for homozygous= rand.nextInt(50) == 0;
		Random rand = new Random();
		boolean homo = rand.nextBoolean(); 

		// homozygous substitution
		if (homo) {
			this.type = "HOMO";

			chr[0] = chrd1_1;
			chr[1] = chrd2_2;
			allele[0] = snp+"";
			allele[1] = allele[0];
			homozygous_genotype(s, c[0]+"", snp+"");//chr[0]=chr[1]
		} 
		else {
			this.type = "HETERO";

			char[][] chr_combi = { { chrd1_1, chrd2_1 }, { chrd1_2, chrd2_2 } };
			// chrd1_1.charAt(check)
			if (snp == chrd2_1) {
				chr[0] = chrd1_2;
				chr[1] = chrd2_2;
				//"0/1"
				allele[0] = c[1]+"";
				allele[1] = snp+"";
				heterozygous_genotype(s, 0, allele[0], allele[1]);
			} 
			else if (snp == chrd2_2) {
				chr[0] = chrd1_1;
				chr[1] = chrd2_1;
				//"1/0";
				allele[0] = snp+"";
				allele[1] = c[0]+"";
				heterozygous_genotype(s, 1, allele[1], allele[0]);
			} 
			else {
				int x = rand.nextInt(2);
				chr[0] = chr_combi[x][0];
				chr[1] = chr_combi[x][1];
				//x=0: 1/0, x=1: 0/1
				if (x == 0) {
					allele[0] = snp+"";
					allele[1] = c[0]+"";
					heterozygous_genotype(s, 1, allele[0], allele[1]);
				}
				else if (x == 1) {
					allele[0] = c[1]+"";
					allele[1] = snp+"";
					heterozygous_genotype(s, 0, allele[0], allele[1]);
				}
			}
		}
		return allele;

	}
	
	private char getSnp(char c) {
		// TODO Auto-generated method stub
		// array list for the other 3 possible bases
//		checkBiallelic();
		ArrayList<Character> possibleBases = new ArrayList<Character>(3);
		// base_list().size() = 4
		for (int k = 0; k < base_list().size(); k++) {
			// if the base from (ATCG) is different!
			// the bases are saved in uppercase
			//c: Upper, c-32: it was lower case, after -32, its uppercase
			if (c != base_list().get(k) && (c-32) != base_list().get(k)) {
				possibleBases.add((char)(base_list().get(k)));
			}
		}
		return randomBase(possibleBases);

	}

	/**
	 * @param s
	 * @param a 
	 */
	private void homozygous_genotype(String s, String ref, String alt) {
		ParentalGenotype pg = new ParentalGenotype(new Genotype(0,0), new Genotype(0,0));
		Allelic_Change ac = new Allelic_Change(ref, alt);
		Genotype child = new Genotype(-1,-1);
		if (s.equals("M"))
			pg.setMother(new Genotype(1,1));
		else if (s.equals("F"))
			pg.setFather(new Genotype(1,1));
		else if (s.equals("MF")) {
			pg.setMother(new Genotype(1,1));
			pg.setFather(new Genotype(1,1));
		}
		else if (s.equals("C")){
			child.setAllele1(1);
			child.setAllele2(1);
			child_ac_arrList.add(ac);
		}
		pg_arrList.add(pg);
		parent_ac_arrList.add(ac);
		g_arrList.add(child);
	}

	/**
	 * @param s
	 * @param x
	 * @param a 
	 */
	private void heterozygous_genotype(String s, int x, String ref, String alt) {
		ParentalGenotype pg = new ParentalGenotype(new Genotype(0,0), new Genotype(0,0));
		Allelic_Change ac = new Allelic_Change(ref, alt);
		Genotype child = new Genotype(-1,-1);
		if (x == 0) {
			if (s.equals("M")) 
				pg.setMother(new Genotype(0,1));		
			else if (s.equals("F")) 
				pg.setFather(new Genotype(0,1));
			
			else if (s.equals("MF")) {
				pg.setMother(new Genotype(0,1));
				pg.setFather(new Genotype(0,1));
			}
			else if (s.equals("C")) {
				child.setAllele1(0);
				child.setAllele2(1);
				child_ac_arrList.add(ac);
			}
		}
		else if (x == 1) {
			if (s.equals("M")) 
				pg.setMother(new Genotype(1,0));
			else if (s.equals("F")) 
				pg.setFather(new Genotype(1,0));
			else if (s.equals("MF")) {
				pg.setMother(new Genotype(1,0));
				pg.setFather(new Genotype(1,0));
			}
			else if (s.equals("C")) {
				child.setAllele1(1);
				child.setAllele2(0);
				child_ac_arrList.add(ac);
			}
		}
		if (!this.type.equals("HETERO"))
			inDelCounter(x,s);

		pg_arrList.add(pg);
		parent_ac_arrList.add(ac);
		g_arrList.add(child);

	}
	
	private void inDelCounter(int x, String s) {
		// TODO Auto-generated method stub
		String a = (x == 0)? 1+"" : 0+"";
//		System.out.println("POP: " + F);
//		System.out.println("MOM: " +M);

		if (s.equals("M")) {
//			System.out.println("---------"+ FrostRunner.parental_chromatids[0] + "\t" + "M"+a + "----------");

			if (FrostRunner.parental_chromatids[0].equals("M"+a)) {
//				System.out.println("Prev MOM: " + this.M);

				if(this.type.equals("IN"))
					this.M++;
				if(this.type.equals("DEL"))
					this.M--;
				
//				System.out.println(this.type);
//				System.out.println("Now MOM: " + this.M);
			}
		}
		if (s.equals("F")) {
//			System.out.println("---------"+ FrostRunner.parental_chromatids[1] + "\t" + "F"+a + "----------");

			if (FrostRunner.parental_chromatids[1].equals("F"+a)) {
//				System.out.println("Prev POP: " + this.F);

				if(this.type.equals("IN"))
					this.F++;
				if(this.type.equals("DEL"))
					this.F--;
//				System.out.println(this.type);
//				System.out.println("Now POP: " + this.F);
			}
		}
			
		if (s.equals("MF")) {
//			System.out.println("---------"+ FrostRunner.parental_chromatids[0] + "\t" + "M"+a + "----------");
//			System.out.println("---------"+ FrostRunner.parental_chromatids[1] + "\t" + "F"+a + "----------");

			if (FrostRunner.parental_chromatids[0].equals("M"+a)) {
//				System.out.println("Prev MOM: " + this.M);

				if(this.type.equals("IN"))
					this.M++;
				if(this.type.equals("DEL"))
					this.M--;
				
//				System.out.println(this.type);
//				System.out.println("Now MOM: " + this.M);
			}
			if (FrostRunner.parental_chromatids[1].equals("F"+a)) {
//				System.out.println("Prev POP: " + this.F);

				if(this.type.equals("IN"))
					this.F++;
				if(this.type.equals("DEL"))
					this.F--;
				
//				System.out.println(this.type);
//				System.out.println("Now POP: " + this.F);
			}
		}
	}
	
	protected String getParentString() {
		String data = "";
		if (this.parent_ac_arrList.size() != 0) {
			for (int i = 0; i < this.parent_ac_arrList.size(); i++) {
				data += this.parent_ac_arrList.get(i).getReference() + "\t" + this.parent_ac_arrList.get(i).getAltered()
						+ "\t" + this.pg_arrList.get(i).getMother().getAllele1() + "/" + this.pg_arrList.get(i).getMother().getAllele2()
						+ "\t" + this.pg_arrList.get(i).getFather().getAllele1() + "/" + this.pg_arrList.get(i).getFather().getAllele2()
						+ "\t" + this.individuum;
			}
		}
		return data;
	}
	
	protected String getChildString() {
		String data = "";
		if (this.child_ac_arrList.size() != 0) {
			for (int i = 0; i < this.child_ac_arrList.size(); i++) {
				data += this.child_ac_arrList.get(i).getReference() + "\t" + this.child_ac_arrList.get(i).getAltered()
						+ "\t" + this.g_arrList.get(i).getAllele1()
						+ "/" + this.g_arrList.get(i).getAllele2()
					 + "\t" + this.individuum;
			}
		}
		return data;
	}

	// private String prettyPrint (final String s) {
	// return "";
	// }

}