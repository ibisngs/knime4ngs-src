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
	protected int F;
	protected int M;
	
	private String individuum = "";
//	@SuppressWarnings("unused")
	private String type = "";
	private ArrayList<ParentalGenotype> pg_arrList = new ArrayList<ParentalGenotype>();
	private ArrayList<Genotype> g_arrList = new ArrayList<Genotype>();//only for child
	private ArrayList<Allelic_Change> parent_ac_arrList = new ArrayList<Allelic_Change>();
	private ArrayList<Allelic_Change> child_ac_arrList = new ArrayList<Allelic_Change>();


	private final double IN_RATE = FrostRunner.insertion_rate;
	private final double DEL_RATE = FrostRunner.deletion_rate;

	public Mutation(char[] basesToMutate, char[] basesNotToMutate, int F, int M) {
		// TODO Auto-generated constructor stub
		setBases(basesToMutate);
		setBases(basesNotToMutate);
		this.F = F;
		this.M = M;
	}
	
	public Mutation(char[] basesToMutate, char[] basesNotToMutate) {
		// TODO Auto-generated constructor stub
		setBases(basesToMutate);
		setBases(basesNotToMutate);
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
	 * @return the pg_arrList
	 */
	protected ArrayList<ParentalGenotype> getPg_arrList() {
		return pg_arrList;
	}

	/**
	 * @return the g_arrList
	 */
	protected ArrayList<Genotype> getG_arrList() {
		return g_arrList;
	}

	/**
	 * @return the ac_arrList
	 */
	protected ArrayList<Allelic_Change> getParent_ac_arrList() {
		return parent_ac_arrList;
	}
	
	/**
	 * @return the ac_arrList
	 */
	protected ArrayList<Allelic_Change> getChild_ac_arrList() {
		return child_ac_arrList;
	}

	/**
	 * @return an arraylist of the 4 bases
	 */
	private final ArrayList<Character> base_list() {
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
	private char randomBase(ArrayList<Character> arrList) {
		Random base = new Random();
		return arrList.get(base.nextInt(arrList.size()));
	}

	/**
	 * @param toMutate
	 * @param f : true if father : GF ; for child false
	 * @param m : true if mother : GM ; for child false
	 * @return
	 */
	protected String[] mutationType(char[] toMutate, char[] toDelete, String s) {

		String[] allele = new String[2];
		Random rand = new Random();
		boolean insert , delete;
		//probability of insert 1/50
		insert = (rand.nextInt((int)(1/IN_RATE)) == 0) ? true : false; 
		// same probability for deletion
		delete = (rand.nextInt((int)(1/DEL_RATE)) == 0) ? true : false;

		//T && T = T, F && F = T
		if ((insert==true && delete==true)||(insert==false && delete==false))		
			allele = substitution(toMutate, s);
		
		if (insert == true && delete == false) 
			allele = insertion(toMutate, s);
			
		
		if (delete == true && insert == false) 
			allele = deletion(toMutate, toDelete, s);
	
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


	private String[] deletion(char[] c, char[] d, String s) {
		// TODO Auto-generated method stub
		this.type = "DEL";

		String[] allele = new String [2];
		String[] chr = new String[2];
		c = getBases();
		
//		char chrd1_1 = 42;//*
//		char chrd2_1 = 42;//*

//		char[][] chr_combi = { { c[0], chrd2_1 },{ chrd1_1, c[1] } };
		String[][] chr_combi = { { (c[0]+"")+""+(d[0]+""), c[1]+"" },{ c[0]+"", (c[1]+"")+""+(d[1]+"") } };

		Random rand = new Random();
		int x = rand.nextInt(2); // either 1 or 2
		chr[0] = chr_combi[x][0];
		chr[1] = chr_combi[x][1];
		allele[0] = chr[0];//+"";
		allele[1] = chr[1];//+"";

		//x=0: 0/1, x=1: 1/0
//		if (x == 0) 		
//			heterozygous_genotype(s, 0, c[1]+"", chr[1]+"");	
//		else if (x == 1) 
//			heterozygous_genotype(s, 1, c[0]+"", chr[0]+"");
		if (x == 0) 		
			heterozygous_genotype(s, 0, (c[1]+"")+""+(d[1]+""), c[1]+"");	
		else if (x == 1) 
			heterozygous_genotype(s, 1, (c[0]+"")+""+(d[0]+""), c[0]+"");	
		return allele;

	}

	private String[] substitution(char[] c, String s) {
		// TODO Auto-generated method stub
		String[] allele = new String [2];
		c = getBases();
		
		char snp = getSnp(c[0]);

		// or set own probability for homozygous= rand.nextInt(50) == 0;
		Random rand = new Random();
		boolean homo = rand.nextBoolean(); 

		// homozygous substitution
		if (homo) {
			this.type = "HOMO";

			allele[0] = snp+"";
			allele[1] = allele[0];
			homozygous_genotype(s, c[0]+"", snp+"");//chr[0]=chr[1]
		} 
		else {
			this.type = "HETERO";

			int x = rand.nextInt(2);
			if (x == 0) {
				allele[0] = c[0]+"";
				allele[1] = getSnp(c[0])+"";
				heterozygous_genotype(s, 0, allele[0], allele[1]);
			}
			else if (x == 1) {
				allele[0] = getSnp(c[1])+"";
				allele[1] = c[1]+"";
				heterozygous_genotype(s, 1, allele[0], allele[1]);
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
		if (s.equals("F")) {
			pg.setFather(new Genotype(1,1));
		}
		else if (s.equals("M"))
			pg.setMother(new Genotype(1,1));
		else if (s.equals("FM")) {
			pg.setFather(new Genotype(1,1));
			pg.setMother(new Genotype(1,1));
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

			if (s.equals("F"))  
				pg.setFather(new Genotype(0,1));

			else if (s.equals("M")) 
				pg.setMother(new Genotype(0,1));
			
			else if (s.equals("FM")) {
				pg.setFather(new Genotype(0,1));
				pg.setMother(new Genotype(0,1));
			
			}
			else if (s.equals("C")) {
				child.setAllele1(0);
				child.setAllele2(1);
				child_ac_arrList.add(ac);
			}
		}
		else if (x == 1) {

			if (s.equals("F")) 
				pg.setFather(new Genotype(1,0));
			
			else if (s.equals("M")) 
				pg.setMother(new Genotype(1,0));
			
			else if (s.equals("FM")) {
				pg.setFather(new Genotype(1,0));
				pg.setMother(new Genotype(1,0));
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

		switch(s) {
		
		case "F":
//			System.out.println("---------"+ FrostRunner.parental_chromatids[0] + "\t" + "F"+a + "----------");

			if (FrostRunner.parental_chromatids[0].equals("F"+a)) {
//				System.out.println("Prev POP: " + this.F);

				if(this.type.equals("IN"))
					this.F++;
				if(this.type.equals("DEL"))
					this.F--;
				
//				System.out.println(this.type);
//				System.out.println("Now POP: " + this.F);
			}
		
		case "M":
//			System.out.println("---------"+ FrostRunner.parental_chromatids[1] + "\t" + "M"+a + "----------");

			if (FrostRunner.parental_chromatids[1].equals("M"+a)) {
//				System.out.println("Prev MOM: " + this.M);

				if(this.type.equals("IN"))
					this.M++;
				if(this.type.equals("DEL"))
					this.M--;
//				System.out.println(this.type);
//				System.out.println("Now MOM: " + this.M);
			}
			
		case "FM":
//			System.out.println("---------"+ FrostRunner.parental_chromatids[0] + "\t" + "F"+a + "----------");
//			System.out.println("---------"+ FrostRunner.parental_chromatids[1] + "\t" + "M"+a + "----------");

			if (FrostRunner.parental_chromatids[0].equals("F"+a)) {
//				System.out.println("Prev POP: " + this.F);

				if(this.type.equals("IN"))
					this.F++;
				if(this.type.equals("DEL"))
					this.F--;
				
//				System.out.println(this.type);
//				System.out.println("Now POP: " + this.F);
			}
			if (FrostRunner.parental_chromatids[1].equals("M"+a)) {
//				System.out.println("Prev MOM: " + this.M);

				if(this.type.equals("IN"))
					this.M++;
				if(this.type.equals("DEL"))
					this.M--;
				
//				System.out.println(this.type);
//				System.out.println("Now MOM: " + this.M);
			}
		}
//		System.out.println(FrostRunner.parental_chromatids[0] + "\t" + "F"+a + "\t" + F);
//		System.out.println(FrostRunner.parental_chromatids[1] + "\t" + "M"+a + "\t" + M);


		
	}

	protected String getParentString() {
		String data = "";
		if (getParent_ac_arrList().size() != 0) {
			for (int i = 0; i < getParent_ac_arrList().size(); i++) {
				data += getParent_ac_arrList().get(i).getReference() + "\t" + getParent_ac_arrList().get(i).getAltered()
						+ "\t" + getPg_arrList().get(i).getMother().getAllele1() + "/" + getPg_arrList().get(i).getMother().getAllele2()
						+ "\t" + getPg_arrList().get(i).getFather().getAllele1() + "/" + getPg_arrList().get(i).getFather().getAllele2()
						+ "\t" + this.individuum;
			}
		}
		return data;
	}
	
	protected String getChildString() {
		String data = "";
		if (getChild_ac_arrList().size() != 0) {
			for (int i = 0; i < getChild_ac_arrList().size(); i++) {
				data += getChild_ac_arrList().get(i).getReference() + "\t" + getChild_ac_arrList().get(i).getAltered()
						+ "\t" + getG_arrList().get(i).getAllele1() + "/" + getG_arrList().get(i).getAllele2()
					 + "\t" + this.individuum;
			}
		}
		return data;
	}

	// private String prettyPrint (final String s) {
	// return "";
	// }

}
