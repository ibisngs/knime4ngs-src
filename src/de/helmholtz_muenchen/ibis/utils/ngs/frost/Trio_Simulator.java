package de.helmholtz_muenchen.ibis.utils.ngs.frost;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Random;


/**
 * @author tanzeem.haque
 *
 * invokeMutation(), invokeMutation_help() and mutation_help()
 * all are only for parents 
 * unchanged() is for both -> improvised as in parents will get 1 input and 2 outputs, 
 * child will get 2 inputs and 1 output
 * 
 * Mental Note: Put the parameters carefully
 *
 *
 */

public class Trio_Simulator {

	//there cannot be more than 8 deNovo mutation
//	private final double recombination_rate = 0.01;

	/**
	 * Attribute
	 */
	private InputScanner in;
	/**
	 * ID_List is needed for Knime in order to submit the output files to art
	 */
	public ArrayList<String> ID_List = new ArrayList<String>();
	
	private ArrayList<String> ID1 = new ArrayList<String>();
	private ArrayList<Integer> POS1 = new ArrayList<Integer>();
	private ArrayList<String> prettyMutation1 = new ArrayList<String>();
	private ArrayList<String> ID2 = new ArrayList<String>();
	private ArrayList<Integer> POS2 = new ArrayList<Integer>();
	private ArrayList<String> prettyMutation2 = new ArrayList<String>();
	private ArrayList<String> recombined = new ArrayList<String>();
	private String[] parentalAlleles;
	private ArrayList<String[]> id_p_allele = new ArrayList<String[]>();


	public Trio_Simulator(InputScanner in) {
		// TODO Auto-generated constructor stub
		setIn(in);
		for(int i = 0; i < getIn().getFs().size(); i++) {
			/**
			 * ID_List will be called in the main method to write the output files.
			 */
			ID_List.add(getIn().getFs().getIdentifier(i));
		}
	}

	/**
	 * @return the in
	 */
	protected InputScanner getIn() {
		return in;
	}

	/**
	 * @param in the in to set
	 */
	protected void setIn(InputScanner in) {
		this.in = in;
	}


	/**
	 * @return the iD
	 */
	protected ArrayList<String> getID1() {
		return ID1;
	}
	/**
	 * @return the pOS
	 */
	protected ArrayList<Integer> getPOS1() {
		return POS1;
	}
	/**
	 * @return the prettyMutation
	 */
	protected ArrayList<String> getPrettyMutation1() {
		return prettyMutation1;
	}

	/**
	 * @return the iD
	 */
	protected ArrayList<String> getID2() {
		return ID2;
	}
	/**
	 * @return the pOS
	 */
	protected ArrayList<Integer> getPOS2() {
		return POS2;
	}
	/**
	 * @return the prettyMutation
	 */
	protected ArrayList<String> getPrettyMutation2() {
		return prettyMutation2;
	}

	/**
	 * @return the recombines
	 */
	protected ArrayList<String> getRecombined() {
		return recombined;
	}
	//a1, a2, id
	protected String[] getParentalAlleles() {
		return parentalAlleles;
	}

	/**
	 * @param al1
	 * @param al2
	 * Id basically not necessary since the combination is already decided from the very beginning
	 */
	protected void setParentalAllele (String al1, String al2){ //, String id) {
		this.parentalAlleles = new String[2];
		this.parentalAlleles[0] = al1;
		this.parentalAlleles[1] = al2;
//		this.parentalAlleles[2] = id;

	}
	/**
	 * @return the id_p_allele
	 */
	protected ArrayList<String[]> getId_p_allele() {
		return id_p_allele;
	}

	/**
	 * @param path
	 * @param mutationRate: this one is the mutation rate of child : deNovo
	 * Basically this method creates the whole trio
	 * @throws InterruptedException
	 */
	protected void createTrio () throws InterruptedException {

		ArrayList<InputData> iData_parents = getIn().getiData_arrList_parents();
		ArrayList<InputData> iData_deNovo_child = getIn().getiData_arrList_deNovo_child();
		ArrayList<InputData> iData_recombination_child = getIn().getiData_arrList_recombination_child();

		long startTime = System.currentTimeMillis();
		generateParents(iData_parents);
		long endTime = System.currentTimeMillis();
		NumberFormat formatter = new DecimalFormat("#0.00000");
		System.out.println("Execution time for MUTATION (parents) "+ formatter.format((endTime - startTime) / 1000d) + " seconds");

		startTime = System.currentTimeMillis();
		generateChild(choose_parentalChromatids(), iData_deNovo_child, iData_recombination_child);
		endTime = System.currentTimeMillis();
		formatter = new DecimalFormat("#0.00000");
		System.out.println("Execution time for DENOVO+RECOMBINATION (child) "+ formatter.format((endTime - startTime) / 1000d) + " seconds");
		
	}
	/**
	 * @param path
	 * @param iData_parents
	 * @param i
	 * creates the parents by inserting mutations
	 */
	private void generateParents(ArrayList<InputData> iData_parents) {

		String id = "";
		String path = FrostRunner.INTERNAL_PATH;
		// parent
		for (int i = 0; i < getIn().getFs().size(); i++) {

			if (iData_parents.get(i).getPositions().size() != 0) {
				id = iData_parents.get(i).getId();

				String[] parent_file = { path + id + "_F_0.fa",
						path + id + "_F_1.fa", path + id + "_M_0.fa",
						path + id + "_M_1.fa" };
				/**
				 * check if file exist, else carry on
				 */
				delete_existing(parent_file);
				manipulate(getIn().getFs().getSequence(i), null, iData_parents.get(i), parent_file, id);

			}
		}
	}
	/**
	 * @param path
	 * @param iData_deNovo_child
	 * @param iData_recombination_child
	 * @param i
	 * Create the child, first denovo then recombination
	 */
	private void generateChild(String[] parental_chromatids, ArrayList<InputData> iData_deNovo_child,
			ArrayList<InputData> iData_recombination_child) {
		
		String id = "";
		String path = FrostRunner.INTERNAL_PATH;
		for (int i = 0; i < getIn().getFs().size(); i++) {

			if (iData_deNovo_child.get(i).getPositions().size() != 0) {

				id = iData_deNovo_child.get(i).getId();
				String[] parent_file = {path+id+parental_chromatids[0], path+id+parental_chromatids[1]};
				String[] child_file = { path + id + "_C_0_no_rec.fa",
						path + id + "_C_1_no_rec.fa", "", "" };

				delete_existing(child_file);

				FastaReader fr = new FastaReader();
				fr.readSequenceFromFile(parent_file[0]);
				String pop = fr.getSequence();

				fr = new FastaReader();
				fr.readSequenceFromFile(parent_file[1]);
				String mom = fr.getSequence();
				/**
				 * create child sequence from mom and dad and insert deNovo
				 */
				manipulate(pop, mom, iData_deNovo_child.get(i), child_file, id);
				// Work for child_file is done

				fr = new FastaReader();
				fr.readSequenceFromFile(child_file[0]);
				String child_0 = fr.getSequence();

				fr = new FastaReader();
				fr.readSequenceFromFile(child_file[1]);
				String child_1 = fr.getSequence();

				String[] child_file_fin = { path + id + "_C_0.fa",
						path + id + "_C_1.fa"};

				delete_existing(child_file_fin);
				/**
				 * recombine child sequence from the non-recombined ones
				 */
				recombine(child_0, child_1, iData_recombination_child.get(i),child_file_fin, id);
				delete_existing(child_file); 
				// deletes the *C_0_no_rec.fa and *C_1_no_rec.fa since not necessary anymore


			}
		}
	}


	/**
	 * @param seq1 : the reference chromosome sequence
	 * @param seq2
	 * @param inputData: (id, ArrayList<Integer>)->id with positions subjected to muattion
	 * @param output00: F0
	 * @param output01: F1
	 * @param output10: M0
	 * @param output11: M1
	 * @param id
	 * @param sequence_length
	 */
	private void manipulate(String seq1, String seq2, InputData inputData, String[] output, String id) {

		String output00 = output[0], output01 = output[1], output10 = output[2], output11 = output[3];
		boolean parent = false;
		if (seq2 == null)
			parent = true;

		appendFile(output00, ">"+id+"\n");
		appendFile(output01, ">"+id+"\n");
		if (parent) {
			appendFile(output10, ">"+id+"\n");
			appendFile(output11, ">"+id+"\n");
		}

		ArrayList<Integer> input_pos = inputData.getPositions();

		//the first position
//		if (input_pos.get(0) == 1) {
//			invokeMutation(char_arrList1, char_arrList2, output00,
//					output01, output10, output11, 0, id, parent);
//		}
		
		for (int i = 0; i < input_pos.size(); i++) {
//			System.out.println("Mutating ... " + (i+1) + ". " + input_pos.get(i) );
			if (i == 0 && input_pos.get(0) == 1) {
				invokeMutation(seq1, seq2, output00,
						output01, output10, output11, i, id, parent);
				continue;
			}
			//if j > 0 => input_pos.get(j) !=1
//			else {
				//extracting the previous positions
			int start = 0;
			if (i > 0)
				start = input_pos.get(i-1);
				//stop is also the real index of that position
				//position is always index+1
			int stop = input_pos.get(i)-1;

			unaltered_sequence(seq1, output00, output10, start, stop, parent, id+"_m_"+parent);
			if (parent)
				unaltered_sequence(seq1,output01, output11, start, stop, parent, id+"_m_"+parent);
			else
				unaltered_sequence(seq2,output01, output11, start, stop, parent, id+"_m_"+parent);
			invokeMutation(seq1, seq2, output00,
						output01, output10, output11, stop, id, parent);
//			}
			//Still something left to append! check the last index
			//inputData.getPositions().size()-1 = my last index
			if ((i == input_pos.size()-1 && input_pos.get(i) < seq1.length())
					|| (seq2 != null && i == input_pos.size()-1
							&& input_pos.get(i) < seq2.length())){
				start = input_pos.get(i);
				stop = seq1.length();
				unaltered_sequence(seq1, output00, output10, start, stop, parent, id+"_m_"+parent);
				if (parent)
					unaltered_sequence(seq1,output01, output11, start, stop, parent, id+"_m_"+parent);
				else
					unaltered_sequence(seq2,output01, output11, start,
							seq2.length(), parent, id+"_m_"+parent);

			}
		}
	}

	private void recombine(String seq1, String seq2, InputData inputData, String[] output, String id) {

		String output00 = output[0], output01 = output[1];
		
		//possible only for child
		appendFile(output00, ">"+id+"\n");
		appendFile(output01, ">"+id+"\n");

		ArrayList<Integer> input_pos = inputData.getPositions();
		//for the first position, if not recombining
		//write the sequence unchanged
		if (input_pos.get(0) > 1) {
			unaltered_sequence(seq1, output00, "", 0, input_pos.get(0)-1, false, id+"_r");
			unaltered_sequence(seq2, output01, "", 0, input_pos.get(0)-1, false, id+"_r");
//			index = 1;
		}
//		else
//			index = 0;
		for (int i = 0; i < input_pos.size()-1; i++) {
				//extracting the previous positions
				int start = input_pos.get(i)-1;
				int stop = input_pos.get(i+1)-1;
				swap(seq1, seq2, output00, output01, start, stop, stop, id+"_r");
		}

		recombine_last_pos(seq1, seq2, input_pos, output00, output01, id+"_r");
	}

	/**
	 * @param seq1
	 * @param seq2
	 * @param input_pos
	 * @param output00
	 * @param output01
	 * @param id
	 */
	private void recombine_last_pos(String seq1, String seq2,
			ArrayList<Integer> input_pos, String output00, String output01, String id) {
		//Still something left to append! check the last index
		//inputData.getPositions().size()-1 = my last index
		//case 2: the last rec position < length of both allele
		//case 2: the last rec position = the last position of either of the allele=> do nothing
		if (input_pos.get(input_pos.size()-1) <= seq1.length()
				&& input_pos.get(input_pos.size()-1) <= seq2.length()){
			int start = input_pos.get(input_pos.size()-1)-1;
			int stop1 = seq1.length();
			int stop2 = seq2.length();
			swap(seq1, seq2, output00, output01, start, stop1, stop2, id);
		}
		else if (input_pos.get(input_pos.size()-1) == seq1.length()
				&& input_pos.get(input_pos.size()-1) < seq2.length()){
			int start = input_pos.get(input_pos.size()-1)-1;
			int stop1 = seq1.length()-1;
			int stop2 = seq2.length();
			swap(seq1, seq2, output00, output01, start, stop1, stop2, id);
		}
		else if (input_pos.get(input_pos.size()-1) < seq1.length()
				&& input_pos.get(input_pos.size()-1) == seq2.length()){
			int start = input_pos.get(input_pos.size()-1)-1;
			int stop1 = seq1.length();
			int stop2 = seq2.length()-1;
			swap(seq1, seq2, output00, output01, start, stop1, stop2, id);
		}
		else
			return;
	}

	/**
	 * @param seq1
	 * @param seq2
	 * @param output00
	 * @param output01
	 * @param rand
	 * @param start
	 * @param stop
	 */
	protected void swap(String seq1, String seq2, String output00, String output01, 
			int start, int stop1, int stop2, String id) {
		Random rand = new Random();
		boolean flip = false;
		flip = rand.nextBoolean();
		if (flip) {
			unaltered_sequence(seq2, output00, "", start, stop2, false, id);
			unaltered_sequence(seq1, output01, "", start, stop1, false, id);
			getRecombined().add(id + "\t" + (start+1) + "\t" + getParentalAlleles()[1] + " " + getParentalAlleles()[0]);
		}
		else {
			unaltered_sequence(seq1, output00, "", start, stop1, false, id);
			unaltered_sequence(seq2, output01, "", start, stop2, false, id);
			getRecombined().add(id + "\t" + (start+1) + "\t" + getParentalAlleles()[0] + " " + getParentalAlleles()[1]);

		}
	}

	/**
	 * @param seq
	 * @param char_arrList2
	 * @param output00
	 * @param output01
	 * @param output10
	 * @param output11
	 * @param start
	 * @param stop
	 * @param parent
	 */
	private void unaltered_sequence(String seq,
			String output00, String output10, int start, int stop, boolean parent, String id) {

		if (stop > seq.length()){
			System.out.println("bad shit is happening");
			System.out.println(stop + " " + seq.length() + " " + id);
			return;
		}

		char[] tmp00 = new char[stop-start];
		char[] tmp10 = new char[stop-start];
		int x = 0;
		for (int j = start; j < stop; j++) {
			if (parent) {
				tmp00[x] = seq.charAt(j);
				tmp10[x] = seq.charAt(j);
			}
			else {
				tmp00[x] = seq.charAt(j);
			}
			x++;
		}
		//writing the bases till previous position is done
		appendFile(output00, new String(tmp00));
		if (parent) {
			appendFile(output10, new String(tmp10));
		}
		
	}


	/**
	 * @param seq1
	 * @param seq2
	 * @param mutationRate
	 * @param output00
	 * @param output01
	 * @param output10
	 * @param output11
	 * @param stop :the position to mutate
	 * @param id
	 * @param parent
	 */
	private void invokeMutation(String seq1, String seq2,
			String output00, String output01, String output10, String output11, int stop, String id, boolean parent) {
		// TODO Auto-generated method stub
		//Now we invoke mutation
		String[] mutated = new String[2];
		String individuum = "";
		char[] toMutate = new char[2];
		if (parent) {
			toMutate[0] = seq1.charAt(stop);
			toMutate[1] = seq1.charAt(stop);
		}
		else {
			toMutate[0] = seq1.charAt(stop);
			toMutate[1] = seq2.charAt(stop);
		}

		Mutation m = new Mutation(toMutate);

		if (parent) {
			Random rand = new Random();
			boolean mutate0 = rand.nextBoolean();
			boolean mutate1 = rand.nextBoolean();
			if (mutate0 == true && mutate1 == false) {
				individuum = "F";
				mutated = m.mutationType(toMutate, individuum);
			} else if (mutate1 == true && mutate0 == false) {
				individuum = "M";
				mutated = m.mutationType(toMutate, individuum);
			} else if ((mutate0 == true && mutate1 == true)
					|| (mutate0 == false && mutate1 == false)) {
				individuum = "FM";
				mutated = m.mutationType(toMutate, individuum);
			}
		}
		else {
			individuum = "C";
			mutated = m.mutationType(toMutate, individuum);
		}
		if (parent) {
			getID1().add(id);
			//because the real position was at the index of the array which is always x-1
			getPOS1().add(stop+1);
			getPrettyMutation1().add(m.getParentString());
		}
		else {
			getID2().add(id);
			getPOS2().add(stop+1);
			getPrettyMutation2().add(m.getChildString());
		}

		write_ref_alt(individuum, mutated, toMutate, output00, output01, output10, output11);
	}
//		else prettyMutation.add("NO MUTATION");

	/**
	 * @param toMutate
	 * @param output112
	 * @param output00: always for F0 or C0
	 * @param output01: always for F1 or C1
	 * @param output10: always for M0
	 * @param output11: always for M1
	 * @param mutated
	 */
	private void write_ref_alt(String individuum, String[] allele, char[] toMutate, String output00,
			String output01, String output10, String output11) {
		//DEL
		String a1 = allele[0], a2 = allele[1];
		char c1 = toMutate[0], c2 = toMutate[1];

		if (allele[0].equals("*"))
			a1 = "";
		if (allele[1].equals("*"))
			a2 = "";

		if (individuum.equals("F")) {
			appendFile(output00, a1);
			appendFile(output01, a2);
			appendFile(output10, c1+"");
			appendFile(output11, c2+"");
		}
		if (individuum.equals("M")) {
			appendFile(output00, c1+"");
			appendFile(output01, c2+"");
			appendFile(output10, a1);
			appendFile(output11, a2);
		}
		if (individuum.equals("FM")) {
			appendFile(output00, a1);
			appendFile(output01, a2);
			appendFile(output10, a1);
			appendFile(output11, a2);
		}
		if (individuum.equals("C")) {
			appendFile(output00, a1);
			appendFile(output01, a2);
		}
	}

	/**
	 * @param fileName
	 * @return the content of the file as a char array
	 */
	@SuppressWarnings("unused")
	private ArrayList<Character> convert_file_to_char_arr(String fileName) {
		ArrayList<Character> char_arrList = new ArrayList<Character> ();
		File file = new File(fileName);
		if (!file.exists()) {
			System.out.println("The file " + fileName + " does not exist");
			System.exit(0);
		}
		try {
			BufferedReader br = new BufferedReader(new FileReader(file));
			String currentLine;
			// There are only 2 lines
			while ((currentLine = br.readLine()) != null) {
				if (currentLine.startsWith(">"))
					currentLine = br.readLine().trim();
				char[] char_arr = new char[currentLine.length()];
				char_arr = currentLine.toCharArray();
				for (int i = 0; i < char_arr.length; i++) {
					char_arrList.add(char_arr[i]);
				}
			}

		} catch (FileNotFoundException f) {
			System.err.println("Can't find file" + file.getPath());
			f.printStackTrace();
		} catch (IOException io) {
			System.err.println("Can't read from file" + file.getPath());
			io.printStackTrace();
		}

		return char_arrList;
	}

	/**
	 * @param fileID : the 4 files from parents
	 * @param id
	 * @return the two files to get the mutated sequences to inherit
	 */
	private String[] choose_parentalChromatids() {
		// TODO Auto-generated method stub
		String[] readFile_for_child = new String[2];
		int calc_prob = probability();

		switch(calc_prob) {

		case 0: //F0M0
			setParentalAllele("F0", "M0");
			this.id_p_allele.add(getParentalAlleles());
			readFile_for_child[0]="_F_0.fa";
			readFile_for_child[1]="_M_0.fa";
			break;
		case 1:  //F0M1
			setParentalAllele("F0", "M1");
			this.id_p_allele.add(getParentalAlleles());
			readFile_for_child[0]="_F_0.fa";
			readFile_for_child[1]="_M_1.fa";
			break;
		case 2: //F1M0
			setParentalAllele("F1", "M0");
			this.id_p_allele.add(getParentalAlleles());
			readFile_for_child[0]="_F_1.fa";
			readFile_for_child[1]="_M_0.fa";
			break;
		case 3: //F1M1
			setParentalAllele("F1", "M1");
			this.id_p_allele.add(getParentalAlleles());
			readFile_for_child[0]="_F_1.fa";
			readFile_for_child[1]="_M_1.fa";
			break;

		}
//		System.out.println("I CAME IN TO CHOOSE FILE");
//		for (String s : readFile_for_child) 
//			System.out.println(s);

		return readFile_for_child;
	}
	/**
	 * Calculates the probability for a random combination of the alleles
	 * random combination of alleles = 2^23 (for human) = 8388608
	 * probability of a diploid combination: 2^23!/((2^23-2)!*2!) = 3.5184367894528*10^13,  100/... = 2.84*10^(-12)
	 */
	private int probability () {
		Random rand;
//		if (in.getSeed() != 0) {
//			rand = new Random(in.getSeed());
//		}
//		else
			rand = new Random();
		int i = rand.nextInt(4); // 4 combinations of chromosomes
		return i;

	}

	@SuppressWarnings("unused")
	private String readSequence(String fileName) {
		
		String seq = "";
		try {
			BufferedReader in = new BufferedReader(new FileReader(fileName));
			StringBuffer buffer = new StringBuffer();
			String currentLine = "";
			if ((currentLine = in.readLine()) == null)
				throw new IOException(fileName + " File is empty ");

			if (currentLine.charAt(0) != '>')
				throw new IOException("Header of the file " + fileName
						+ " must be '>'");

			for (currentLine = in.readLine().trim(); currentLine != null; currentLine = in
					.readLine()) {
				if (currentLine.length() > 0 && currentLine.charAt(0) == '>') {
					// do nothing
					System.out.println("I CAME IN HEADER");
					continue;
				} else {
//					System.out.println("I CAME IN SEQUENCE");
					buffer.append(currentLine.trim());
				}
			}
			if (buffer.length() != 0)
				seq = buffer.toString();

		} catch (IOException e) {
			System.out.println("Error when reading " + fileName);
			e.printStackTrace();

		}

		// System.out.println(seq);
		return seq;
	}

	/**
	 * @param file
	 * Delete if file exists
	 * 
	 */
	private void delete_existing(String[] file) {
		try {
			for (String s : file) {
				File fileTmp = new File(s);
				if (fileTmp.exists()) {
					fileTmp.delete();
				}
			}

		} catch (Exception e) {
			// if any error occurs
			e.printStackTrace();
		}
	}

	/**
	 * @param fileName
	 * What do you think it does?
	 */
	private void appendFile(String fileName, String s) {
		// TODO Auto-generated method stub
		try (PrintWriter pw = new PrintWriter(new FileOutputStream (fileName, true))) { // new PrintWriter(new BufferedWriter(new FileWriter(fileName, true))))
			pw.println(s);
		} catch (Exception e) {
			System.err.println("Error: " + e.getMessage());
			e.printStackTrace();
		}
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	protected String getParentInfo() {
		String data = "";

		if(getPrettyMutation1().size() != 0) {
			for (int i = 0; i < getPrettyMutation1().size(); i++) {
				data += getID1().get(i) + "\t" + getPOS1().get(i) + "\t" + getPrettyMutation1().get(i) + "\n";

			}
		}
		return data;
	}

	protected String getChildInfo() {

		String data = "";
		if(getPrettyMutation2().size() != 0) {
			for (int i = 0; i < getPrettyMutation2().size(); i++) {
				data += getID2().get(i) + "\t" + getPOS2().get(i) + "\t" + getPrettyMutation2().get(i) + "\t";
				for (int j = 0; j < id_p_allele.size(); j++) {
//					if (getID2().get(i).equals(id_p_allele.get(j)[2])) {
						data += id_p_allele.get(j)[0] + "\t" + id_p_allele.get(j)[1]+ "\n";
//					}
				}

			}
		}
		return data;
	}
}