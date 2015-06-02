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

public class TrioSimulator {

	//there cannot be more than 8 deNovo mutation
//	private final double recombination_rate = 0.01;


	int current_chunk;
	/**
	 * Attribute
	 */
	private FastaCheck fc;
	private FastaReader fr;
	private InputScanner in;
	/**
	 * record for the adjustment of the sequence for denovo mutation later
	 */
	private ArrayList<ArrayList<Integer>> sequenceAdjustment = new ArrayList<>();
	/**
	 * ID_List is needed for Knime in order to submit the output files to art
	 */	
	private ArrayList<String> ID1 = new ArrayList<String>();
	private ArrayList<Integer> POS1 = new ArrayList<Integer>();
	private ArrayList<String> prettyMutation1 = new ArrayList<String>();
	private ArrayList<String> ID2 = new ArrayList<String>();
	private ArrayList<Integer> POS2 = new ArrayList<Integer>();
	private ArrayList<String> prettyMutation2 = new ArrayList<String>();
	private ArrayList<String> recombined = new ArrayList<String>();


	public TrioSimulator(FastaCheck fc, FastaReader fr, InputScanner in) {
		// TODO Auto-generated constructor stub
		this.fc = fc;
		this.fr = fr;
		this.in = in;
	}
/*
	public FastaCheck getFc() {
		return fc;
	}

	protected FastaReader getFr() {
		return fr;
	}

	protected InputScanner getIn() {
		return in;
	}

	protected int getCurrentLength() {
		return currentLength;
	}
*/
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


	/**
	 * @param currentChr 
	 * @param chunk 
	 * @param chunk 
	 * @param parentalChromatids 
	 * @param currentLength 
	 * @param path
	 * @param mutationRate: this one is the mutation rate of child : deNovo
	 * Basically this method creates the whole trio
	 * @throws InterruptedException
	 * @throws IOException 
	 */
	protected void createTrio (int chunk_index, String currentChr) throws InterruptedException, IOException {

		InputData iData_parents = this.in.getiData_parents();
		InputData iData_deNovo_child = this.in.getiData_deNovo_child();
		InputData iData_recombination_child = this.in.getiData_recombination_child();
		/**
		 * DO NOT DELETE THIS
		 */
		this.current_chunk = chunk_index;

//		long startTime = System.currentTimeMillis();
		generateParents(iData_parents);
//		long endTime = System.currentTimeMillis();
//		NumberFormat formatter = new DecimalFormat("#0.00000");
//		System.out.println("Execution time for MUTATION (parents) "+ formatter.format((endTime - startTime) / 1000d) + " seconds");
		
//		startTime = System.currentTimeMillis();
		generateChild(this.fc.choose_parentalChromatids(), iData_deNovo_child, iData_recombination_child, currentChr);
//		endTime = System.currentTimeMillis();
//		formatter = new DecimalFormat("#0.00000");
//		System.out.println("Execution time for DENOVO+RECOMBINATION (child) "+ formatter.format((endTime - startTime) / 1000d) + " seconds");

	}
	/**
	 * @param path
	 * @param iData_parents
	 * @param chunk_index 
	 * creates the parents by inserting mutations
	 * @param chunk 
	 */
	private void generateParents(InputData iData_parents) {

		String id = "";
		String path = FrostRunner.INTERNAL_PATH;
		// parent

		if (iData_parents.getPositions().size() != 0) {
			id = this.current_chunk + "_" + iData_parents.getId();
//			System.out.println("FILE: "+ id);


			String[] parent_file = { path + id + "_F_0.fa",
					path + id + "_F_1.fa", path + id + "_M_0.fa",
					path + id + "_M_1.fa" };
			/**
			 * check if file exist, else carry on
			 */
			delete_existing(parent_file);
			int begin = this.current_chunk*FrostRunner.chunk_length;
			int end = begin + FrostRunner.chunk_length;
			if (end > this.fr.getLength())
				end = this.fr.getLength();
//			System.out.println("BEGIN: " + begin + "\t" + "END: " + end);
			mutate(this.fr.getSequence().substring(begin, end), null, iData_parents, parent_file, id, "m");
//			System.out.println("FASTA LENGTH TOTAL: " + this.fr.getLength());
//			System.out.println("FASTA LENGTH PARENTS: " + this.fr.getSequence().substring(begin, end).length());
		}
	}
	/**
	 * @param path
	 * @param iData_deNovo_child
	 * @param iData_recombination_child
	 * @param currentChr 
	 * @param chunk 
	 * @param i 
	 * @param i
	 * Create the child, first denovo then recombination
	 * @throws IOException 
	 */
	private void generateChild(String[] parental_chromatids, InputData iData_deNovo_child,
			InputData iData_recombination_child, String currentChr) throws IOException {
		
		String id = "";
		String path = FrostRunner.INTERNAL_PATH;

		if (iData_deNovo_child.getPositions().size() != 0) {

			id = this.current_chunk + "_" + iData_deNovo_child.getId();
			String[] parent_file = {path + id + parental_chromatids[0], path+id+parental_chromatids[1]};
			String[] child_file = { path + id + "_C_0_no_rec.fa",
					path + id + "_C_1_no_rec.fa", "", "" };

			delete_existing(child_file);

			FastaReader fr = new FastaReader();
			fr.readSequenceFromFile(parent_file[0], id);
			String pop = fr.getSequence();

//			System.out.println("FASTA LENGTH CHILD (pop): " + pop.length());

			fr = new FastaReader();
			fr.readSequenceFromFile(parent_file[1], id);
			String mom = fr.getSequence();
			
//			System.out.println("FASTA LENGTH CHILD (mom): " + mom.length());

			/**
			 * create child sequence from mom and dad and insert deNovo
			 */
			mutate(pop, mom, iData_deNovo_child, child_file, id, "d");
			// Work for child_file is done

			fr = new FastaReader();
			fr.readSequenceFromFile(child_file[0], id);
			String child_0 = fr.getSequence();

			fr = new FastaReader();
			fr.readSequenceFromFile(child_file[1], id);
			String child_1 = fr.getSequence();
			
			String[] child_file_fin = { path + id + "_C_0.fa",
					path + id + "_C_1.fa"};

//			System.out.println("FASTA LENGTH CHILD (nr1): " + child_0.length());
//			System.out.println("FASTA LENGTH CHILD (nr2): " + child_1.length());


			delete_existing(child_file_fin);
			/**
			 * recombine child sequence from the non-recombined ones
			 * 		if (iData_recombination_child.getPositions().size() != 0) {
			 */
			recombine(child_0, child_1, iData_recombination_child,child_file_fin, id);
			delete_existing(child_file); 
			// deletes the *C_0_no_rec.fa and *C_1_no_rec.fa since not necessary anymore
	
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
	 * @param chunk_idx 
	 * @param sequence_length
	 */
	private void mutate(String seq1, String seq2, InputData inputData, String[] output, String id, String mutation_tag) {

		ArrayList<ArrayList<Integer>> input_pos = inputData.getPositions();
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


		/**
		 * it is ok to check only seq1, (most of the time, this is the chunk length)
		 * 1) for parents, seq2 is null
		 * 2) for child, seq2 is never smaller than de novo mut position, 
		 * since it had insertions inherited from parents
		 
		int mutation_idx = 0; //, tmp_idx = this.current_chunk*FrostRunner.chunk_length;

		if (mutation_tag.equals("m"))
			mutation_idx = FrostRunner.mutation_index_parent;
		else if (mutation_tag.equals("d"))
			mutation_idx = FrostRunner.denovo_index_child;
		*/
		int i = this.current_chunk;
//		System.out.println("chunky: " + i);
		int tmp_idx = i*FrostRunner.chunk_length;
		
		if (input_pos.get(i).size()==0) {
			int start = 0;
			int stop = seq1.length();
//			System.out.println("IS EMPTY: " + start + "\t" + stop);
			unaltered_sequence(seq1, output00, output10, start, stop, parent, id+"_m_"+parent);
			if (parent)
				unaltered_sequence(seq1,output01, output11, start, stop, parent, id+"_m_"+parent);
			else
				unaltered_sequence(seq2,output01, output11, start,
						seq2.length(), parent, id+"_m_"+parent);
			return;
		}

		for (int j = 0; j < input_pos.get(i).size(); j++) {

			if (j == 0 && input_pos.get(i).get(0) == 1) {
				invokeMutation(seq1, seq2, output00, output01, output10, output11, 0, id, parent);
				continue;
			}
			//extracting the previous positions
			int start = 0;
			if (j > 0)
				start = input_pos.get(i).get(j-1)-tmp_idx;
			//stop is also the real index of that position
			//position is always index+1
			int stop = input_pos.get(i).get(j)-tmp_idx-1;
			/**
			 * last chunk
			 */
			if (i == this.in.getChunk()-1 && (j == 0)) {
					start = 0;
			}
			unaltered_sequence(seq1, output00, output10, start, stop, parent, id+"_m_"+parent);
			if (parent)
				unaltered_sequence(seq1,output01, output11, start, stop, parent, id+"_m_"+parent);
			else
				unaltered_sequence(seq2,output01, output11, start, stop, parent, id+"_m_"+parent);
			/**
			 * checking
			 */
			invokeMutation(seq1, seq2, output00, output01, output10, output11, stop, id, parent);
			//Still something left to append! check the last index
			//inputData.getPositions().size()-1 = my last index
			if ((j == input_pos.get(i).size()-1 && input_pos.get(i).get(j)-tmp_idx < seq1.length())
					|| (seq2 != null && j == input_pos.get(i).size()-1 && input_pos.get(i).get(j)-tmp_idx < seq2.length())){
				start = input_pos.get(i).get(j) - tmp_idx;
				stop = seq1.length();
//				System.out.println(start + "\t" + stop);
				unaltered_sequence(seq1, output00, output10, start, stop, parent, id+"_m_"+parent);
				if (parent)
					unaltered_sequence(seq1,output01, output11, start, stop, parent, id+"_m_"+parent);	
				else
					unaltered_sequence(seq2,output01, output11, start, seq2.length(), parent, id+"_m_"+parent);
			}
		}
	}
	
	private void recombine(String seq1, String seq2, InputData inputData, String[] output, String id) {

		String output00 = output[0], output01 = output[1];
		//possible only for child
		appendFile(output00, ">"+id+"\n");
		appendFile(output01, ">"+id+"\n");

		ArrayList<ArrayList<Integer>> input_pos = inputData.getPositions();

		int i = this.current_chunk;
		int tmp_idx = i*FrostRunner.chunk_length;
		
		if (input_pos.get(i).size()==0) {
			int start = 0;
			int stop = seq1.length();
//			System.out.println("IS EMPTY: " + start + "\t" + stop);
			unaltered_sequence(seq1, output00, "", start, stop, false, id+"_r");
			unaltered_sequence(seq2,output01, "", start, seq2.length(), false, id+"_r");
			return;
		}
		//for the first position, if not recombining
		//write the sequence unchanged
		for (int j = 0; j < input_pos.get(i).size(); j++) {
			if (j == 0 && input_pos.get(i).get(j) > 1+tmp_idx) {
				unaltered_sequence(seq1, output00, "", 0, input_pos.get(i).get(0)-tmp_idx-1, false, id+"_r");
				unaltered_sequence(seq2, output01, "", 0, input_pos.get(i).get(0)-tmp_idx-1, false, id+"_r");
//				System.out.println("RECOMBINE first: " + "\t" + 0 + "\t" + (input_pos.get(i).get(0)-tmp_idx-1));
			}	
			if (j < input_pos.get(i).size()-1) {
				//extracting the previous positions
				int start_idx = input_pos.get(i).get(j)-tmp_idx-1;
				/**
				 * last chunk
				 */
//				if (i == this.in.getChunk()-1 && (j == 0)) {
//					start_idx = 0;
//				}
				int stop_idx = input_pos.get(i).get(j+1)-tmp_idx-1;
				swap(seq1, seq2, output00, output01, start_idx, stop_idx, stop_idx, id+"_r", tmp_idx);
//				System.out.println("RECOMBINE: " + "\t" + start_idx + "\t" + stop_idx);
			}
			if (j == input_pos.get(i).size()-1)
				recombine_last_pos(seq1, seq2, input_pos.get(i), i, output00, output01, id+"_r");

		}		
	}

	/**
	 * @param seq1
	 * @param seq2
	 * @param input_pos_sublist
	 * @param output00
	 * @param output01
	 * @param id
	 */
	private void recombine_last_pos(String seq1, String seq2, ArrayList<Integer> input_pos_sublist, int idx, String output00, String output01, String id) {
		//Still something left to append! check the last index
		//inputData.getPositions().size()-1 = my last index
		//case 2: the last rec position < length of both allele
		//case 2: the last rec position = the last position of either of the allele=> do nothing
		int tmp_idx = idx*FrostRunner.chunk_length;
		int start = input_pos_sublist.get(input_pos_sublist.size()-1)-tmp_idx-1, stop1 = 0, stop2 = 0;

		if (input_pos_sublist.get(input_pos_sublist.size()-1)-tmp_idx <= seq1.length()
				&& input_pos_sublist.get(input_pos_sublist.size()-1)-tmp_idx <= seq2.length()){
			stop1 = seq1.length();
			stop2 = seq2.length();
			swap(seq1, seq2, output00, output01, start, stop1, stop2, id, tmp_idx);
//			System.out.println("RECOMBINE last: " + "\t" + start + "\t" + stop1 + "\t" + stop2);

		}
		else if (input_pos_sublist.get(input_pos_sublist.size()-1)-tmp_idx == seq1.length()
				&& input_pos_sublist.get(input_pos_sublist.size()-1)-tmp_idx < seq2.length()){
			stop1 = seq1.length()-1;
			stop2 = seq2.length();
			swap(seq1, seq2, output00, output01, start, stop1, stop2, id, tmp_idx);
//			System.out.println("RECOMBINE last: " + "\t" + start + "\t" + stop1 + "\t" + stop2);

		}
		else if (input_pos_sublist.get(input_pos_sublist.size()-1)-tmp_idx < seq1.length()
				&& input_pos_sublist.get(input_pos_sublist.size()-1)-tmp_idx == seq2.length()){
			stop1 = seq1.length();
			stop2 = seq2.length()-1;
			swap(seq1, seq2, output00, output01, start, stop1, stop2, id, tmp_idx);
//			System.out.println("RECOMBINE last: " + "\t" + start + "\t" + stop1 + "\t" + stop2);

		}

		else {
//			System.out.println("RECOMBINE last: DEI MUDDA " + (idx*FrostRunner.chunk_length) + "\t"+ (input_pos_sublist.get(input_pos_sublist.size()-1)- tmp_idx));
			return;
		}
	}
	/**
	 * @param seq1
	 * @param seq2
	 * @param output00
	 * @param output01
	 * @param rand
	 * @param start
	 * @param tmp_idx2 
	 * @param stop
	 */
	protected void swap(String seq1, String seq2, String output00, String output01, 
			int start, int stop1, int stop2, String id, int tmp_idx) {
		Random rand = new Random();
		boolean flip = false;
		flip = rand.nextBoolean();
		if (flip) {
			unaltered_sequence(seq2, output00, "", start, stop2, false, id);
			unaltered_sequence(seq1, output01, "", start, stop1, false, id);
			getRecombined().add(id + "\t" + (start+tmp_idx+1) + "\t" + FrostRunner.parental_chromatids[1] + " " + FrostRunner.parental_chromatids[0]);
		}
		else {
			unaltered_sequence(seq1, output00, "", start, stop1, false, id);
			unaltered_sequence(seq2, output01, "", start, stop2, false, id);
			getRecombined().add(id + "\t" + (start+tmp_idx+1) + "\t" + FrostRunner.parental_chromatids[0] + " " + FrostRunner.parental_chromatids[1]);

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
	private void unaltered_sequence(String seq, String output00, String output10, 
			int start, int stop, boolean parent, String id) {
		
		if (stop > seq.length()){
			System.out.println("bad shit is happening");
			System.out.println(stop + " " + seq.length() + " " + id);
			return;
		}

//		System.out.println("Start: " + start + "\t" + "Stop: " + stop);
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
			String output00, String output01, String output10, String output11, 
			int stop, String id, boolean parent) {
		// TODO Auto-generated method stub
		ArrayList<Integer> adjust = new ArrayList<>(2);

		//Now we invoke mutation
		int tmp_idx = this.current_chunk*FrostRunner.chunk_length;
		String[] mutated = new String[2];
		String individuum = "";
		char[] toMutate = new char[2], toDelete = new char[2];

		if (parent) {
			
			if (seq1.charAt(stop+1) > tmp_idx) 		
				stop = stop -1;
							
			toMutate[0] = seq1.charAt(stop);
			toMutate[1] = seq1.charAt(stop);
			toDelete[0] = seq1.charAt(stop+1);
			toDelete[1] = seq1.charAt(stop+1);
		}
		else {
			int stopAdjusted = adjustDeNovo(stop);
			while (stopAdjusted >= tmp_idx) {
				int diff = stopAdjusted-tmp_idx;
				stop = stop-diff;
				stopAdjusted = adjustDeNovo(stop);
			}
			if (seq1.charAt(stopAdjusted+1) > tmp_idx || (seq2.charAt(stopAdjusted+1) > tmp_idx)) 
				stopAdjusted = stopAdjusted - 1;
				
			toMutate[0] = seq1.charAt(stopAdjusted);
			toMutate[1] = seq2.charAt(stopAdjusted);
			toDelete[0] = seq1.charAt(stopAdjusted+1);
			toDelete[1] = seq2.charAt(stopAdjusted+1);
		}

		Mutation m = new Mutation(toMutate, toDelete);

		if (parent) {
						
			Random rand = new Random();
			boolean mutate0 = rand.nextBoolean();
			boolean mutate1 = rand.nextBoolean();
			if (mutate0 == true && mutate1 == false) {
				individuum = "F";
				mutated = m.mutationType(toMutate, toDelete, individuum);
			} else if (mutate1 == true && mutate0 == false) {
				individuum = "M";
				mutated = m.mutationType(toMutate, toDelete, individuum);
			} else if ((mutate0 == true && mutate1 == true)
					|| (mutate0 == false && mutate1 == false)) {
				individuum = "FM";
				mutated = m.mutationType(toMutate, toDelete, individuum);
			}
			getID1().add(id);
			//because the real position was at the index of the array which is always x-1
			getPOS1().add(stop+tmp_idx+1);
			getPrettyMutation1().add(m.getParentString());
			/**
			 * record for denovo adjustment
			 */
			adjust.add(stop);
			adjust.add(FrostRunner.mutationCounter);
			this.sequenceAdjustment.add(adjust);
		}
		else {

			individuum = "C";
			mutated = m.mutationType(toMutate, toDelete, individuum);
			getID2().add(id);
			getPOS2().add(stop+tmp_idx+1);
			getPrettyMutation2().add(m.getChildString());
		}		
		write_ref_alt(individuum, mutated, toMutate, output00, output01, output10, output11);
	}
//		else prettyMutation.add("NO MUTATION");

	private int adjustDeNovo(int position) {
		// TODO Auto-generated method stub
		int p = 0;
		//+tmp_idx+1
		for (int k = 0; k < this.sequenceAdjustment.size(); k++) {
			if (this.sequenceAdjustment.get(k).get(0) == position) {
				p = position + this.sequenceAdjustment.get(k).get(1);
				break;
			}
			else if (this.sequenceAdjustment.get(k).get(0) > position) {
				if (k > 0) 
					p = position + this.sequenceAdjustment.get(k-1).get(1);
				else 
					p = position;
				break;
			}
			else
				continue;						
		}
		return p;
	}
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
			pw.write(s);
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

		/**
		 * CHILD ALWAYS GETS M/F AND NOT F/M
		 */
		String data = "";
		if(getPrettyMutation2().size() != 0) {
			for (int i = 0; i < getPrettyMutation2().size(); i++) {
				if (getPrettyMutation2().get(i) != "") {
					data += getID2().get(i) + "\t" + getPOS2().get(i) + "\t" + getPrettyMutation2().get(i) + "\t";
					data += FrostRunner.parental_chromatids[1] + "\t" + FrostRunner.parental_chromatids[0]+ "\n";
				}
				else
					continue;			

			}
		}
		return data;
	}

}