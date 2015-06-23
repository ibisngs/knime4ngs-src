package de.helmholtz_muenchen.ibis.utils.ngs.frost;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
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
	private ArrayList<ArrayList<Integer>> sequenceAdjustment = new ArrayList<>();


//	int current_chunk;
	/**
	 * Attribute
	 */
//	private FastaCheck fc;
	private FastaReader fr;
	private InputScanner in;
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


	public TrioSimulator(/*FastaCheck fc,*/ FastaReader fr, InputScanner in) {
		// TODO Auto-generated constructor stub
//		this.fc = fc;
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
//	protected ArrayList<String> getID1() {
//		return ID1;
//	}
	/**
	 * @return the pOS
	 */
//	protected ArrayList<Integer> getPOS1() {
//		return POS1;
//	}
	/**
	 * @return the prettyMutation
	 */
//	protected ArrayList<String> getPrettyMutation1() {
//		return prettyMutation1;
//	}

	/**
	 * @return the iD
	 */
//	protected ArrayList<String> getID2() {
//		return ID2;
//	}
	/**
	 * @return the pOS
	 */
//	protected ArrayList<Integer> getPOS2() {
//		return POS2;
//	}
	/**
	 * @return the prettyMutation
	 */
//	protected ArrayList<String> getPrettyMutation2() {
//		return prettyMutation2;
//	}

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
	protected void createTrio (/*int chunk_index, */String currentChr) throws InterruptedException, IOException {

		InputData iData_parents = this.in.getiData_parents();
		InputData iData_deNovo_child = this.in.getiData_deNovo_child();
		InputData iData_recombination_child = this.in.getiData_recombination_child();
		/**
		 * DO NOT DELETE THIS
		 */
//		this.current_chunk = chunk_index;

		generateParents(iData_parents);

		/*
		for (int i = 0; i < this.sequenceAdjustment.size(); i++)
			{
//				if (this.sequenceAdjustment.get(i).get(1) != 0 || this.sequenceAdjustment.get(i).get(2) != 0) 
					System.out.println("Adjustment" + "\n" 
							+ this.sequenceAdjustment.get(i).get(0) + "\t" 
							+ this.sequenceAdjustment.get(i).get(1) + "\t"
							+ this.sequenceAdjustment.get(i).get(2));
			}*/
		
		generateChild(FrostRunner.parental_chromatids, iData_deNovo_child, iData_recombination_child, currentChr);

	}
	/**
	 * @param path
	 * @param iData_parents
	 * @param chunk_index 
	 * creates the parents by inserting mutations
	 * @param chunk 
	 * @throws IOException 
	 */
	private void generateParents(InputData iData_parents) throws IOException {

		String id = "";
		String path = FrostRunner.INTERNAL_PATH;
		// parent

		if (iData_parents.getPositions().size() != 0) {
			id = /*this.current_chunk + "_" + */iData_parents.getId();
//			System.out.println("FILE: "+ id);


			String[] parent_file = { path + id + "_M0.fa",
					path + id + "_M1.fa", path + id + "_F0.fa",
					path + id + "_F1.fa" };
			/**
			 * check if file exist, else carry on
			 */
			delete_existing(parent_file);
			/*
			int begin = this.current_chunk*FrostRunner.chunk_length;
			int end = begin + FrostRunner.chunk_length;
			if (end > this.fr.getLength())
				end = this.fr.getLength();
				*/
//			System.out.println("BEGIN: " + begin + "\t" + "END: " + end);
			mutate(this.fr.getSequence()/*.substring(begin, end)*/, null, iData_parents, parent_file, id, "m");
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

			id = /*this.current_chunk + "_" + */iData_deNovo_child.getId();
			String[] parent_file = {path + id + "_"+parental_chromatids[0]+".fa", 
					path+id+"_"+parental_chromatids[1]+".fa"};
			String[] child_file = { path + id + "_C0_no_rec.fa",
					path + id + "_C1_no_rec.fa", "", "" };

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
			
			String[] child_file_fin = { path + id + "_C0.fa",
					path + id + "_C1.fa"};

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
	 * @throws IOException 
	 */
	private void mutate(String seq1, String seq2, InputData inputData, String[] output, String id, String mutation_tag) throws IOException {
		
		boolean parent = false;
		ArrayList/*<ArrayList*/<Integer>/*>*/ input_pos = inputData.getPositions();
		
		BufferedWriter output00 = new BufferedWriter(new FileWriter(output[0], true), 10000000)/*M0*/, 
				output01 = new BufferedWriter(new FileWriter(output[1], true), 10000000)/*M1*/,
				output10 = null, output11 = null;
		
				
		if (seq2 == null) {
			parent = true;
			output10 = new BufferedWriter(new FileWriter(output[2], true), 10000000)/*F0*/;
			output11 = new BufferedWriter(new FileWriter(output[3], true), 10000000)/*F1*/;
		}

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
		/*
		int i = this.current_chunk;
		System.out.println("chunky: " + i);
		int tmp_idx = i*FrostRunner.chunk_length;
		*/
		/**
		 * no mutation here, also no change for child
		 */
		if (input_pos./*get(i).*/size()==0) {
			int start = 0;
			int stop = seq1.length();
//			System.out.println("IS EMPTY: " + start + "\t" + stop);
			unaltered_sequence(seq1, output00, output10, start, stop, parent, id+"_m_"+parent); //M0, F0
			if (parent)
				unaltered_sequence(seq1,output01, output11, start, stop, parent, id+"_m_"+parent); //M1, F1
			else
				unaltered_sequence(seq2,output01, output11, start,
						seq2.length(), parent, id+"_m_"+parent); //
			return;
		}

		for (int j = 0; j < input_pos./*get(i).*/size(); j++) {
//			if(parent)
//				System.out.println("Mut: " + input_pos.get(j));
			if (j == 0 && input_pos./*get(i).*/get(0) == 1) {
				int stop = 0 /*for parents, and if child, mother allele*/,
						stop_2 = 0 /*not for parents, only for child father allele*/;
				if (!parent){
					stop = adjustDeNovo(stop)[0]/*M*/;
					stop_2 = adjustDeNovo(stop)[1]/*F*/; 
				}
				invokeMutation(seq1, seq2, output00, output01, output10, output11, stop, 
						stop_2/*will not be taken for parents*/, id, parent);
				 
					
				continue;
			}
			//extracting the previous positions
			int start = 0, start_2 = 0;
			if (j > 0)
				start = input_pos./*get(i).*/get(j-1)+1/*-tmp_idx*/;
			//stop is also the real index of that position
			//position is always index+1
			int stop = input_pos./*get(i).*/get(j)/*-tmp_idx-1*/, stop_2 = 0;
			
			if(!parent)
			{
				start = adjustDeNovo(start)[0];/*for Child start, mother allele*/
				start_2 = adjustDeNovo(start)[1];/*for Child start, father allele*/
				stop = adjustDeNovo(stop)[0];/*for Child stop, mother allele*/
				stop_2 = adjustDeNovo(stop)[1];/*for Child stop, father allele*/

			}
			/**
			 * Mother or mother allele for child
			 */
			unaltered_sequence(seq1, output00, output10, start, stop, parent, id+"_m_"+parent);//M0, F0
			if (parent)
				/**
				 * father
				 */
				unaltered_sequence(seq1,output01, output11, start, stop, parent, id+"_m_"+parent);
			else
				/**
				 * father allele for child
				 */
				unaltered_sequence(seq2,output01, output11, start_2, stop_2, parent, id+"_m_"+parent);
			/**
			 * checking
			 */
			invokeMutation(seq1, seq2, output00, output01, output10, output11, stop,
					stop_2/*will not be taken into account if parents*/, id, parent);
			//Still something left to append! check the last index
			if ((j == input_pos./*get(i).*/size()-1 && input_pos./*get(i).*/get(j)/*-tmp_idx*/ < seq1.length())
					|| (seq2 != null && j == input_pos./*get(i).*/size()-1 && input_pos./*get(i).*/get(j)/*-tmp_idx*/ < seq2.length())){
				start = input_pos./*get(i).*/get(j)+1 /*- tmp_idx*/;
				stop = seq1.length(); 
				if(!parent)
				{
					start = adjustDeNovo(start)[0];/*for Child start, mother allele*/
					start_2 = adjustDeNovo(start)[1];/*for Child start, father allele*/
				}
				/**
				 * Mother or mother allele for child
				 */
				unaltered_sequence(seq1, output00, output10, start, stop, parent, id+"_m_"+parent); //M0, F0
				if (parent)
					/**
					 * father allele
					 */
					unaltered_sequence(seq1,output01, output11, start, stop, parent, id+"_m_"+parent);//M1, F1	
				else
					/**
					 * father allele for child
					 */
					unaltered_sequence(seq2,output01, output11, start_2, seq2.length(), parent, id+"_m_"+parent);
			}
		}
		output00.close();
		output01.close();
		if (parent) {
			output10.close();
			output11.close();
		}
		
	}
	
	private void recombine(String seq1, String seq2, InputData inputData, String[] output, String id) throws IOException {

		BufferedWriter output00 = new BufferedWriter(new FileWriter(output[0], true), 10000000)/*C0*/, 
			output01 = new BufferedWriter(new FileWriter(output[1], true), 10000000)/*C1*/;
		//possible only for child
		appendFile(output00, ">"+id+"\n");
		appendFile(output01, ">"+id+"\n");

		ArrayList/*<ArrayList*/<Integer>/*>*/ input_pos = inputData.getPositions();
/*
		int i = this.current_chunk;
		int tmp_idx = i*FrostRunner.chunk_length;
*/		
		if (input_pos./*get(i).*/size()==0) {
			int start = 0;
			int stop = seq1.length();
//			System.out.println("IS EMPTY: " + start + "\t" + stop);
			unaltered_sequence(seq1, output00, null, start, stop, false, id+"_r");//C0
			unaltered_sequence(seq2,output01, null, start, seq2.length(), false, id+"_r");//C1
			return;
		}
		//for the first position, if not recombining
		//write the sequence unchanged
		for (int j = 0; j < input_pos./*get(i).*/size(); j++) {
			if (j == 0 && input_pos./*get(i).*/get(j) > 1/*+tmp_idx*/) {
				unaltered_sequence(seq1, output00, null, 0, input_pos./*get(i).*/get(0)/*-tmp_idx*/-1, false, id+"_r");
				unaltered_sequence(seq2, output01, null, 0, input_pos./*get(i).*/get(0)/*-tmp_idx*/-1, false, id+"_r");
			}	
			if (j < input_pos./*get(i).*/size()-1) {
				//extracting the previous positions
				int start_idx = input_pos./*get(i).*/get(j)/*-tmp_idx*/-1;
				/**
				 * last chunk
				 */
//				if (i == this.in.getChunk()-1 && (j == 0)) {
//					start_idx = 0;
//				}
				int stop_idx = input_pos./*get(i).*/get(j+1)/*-tmp_idx*/-1;
				swap(seq1, seq2, output00, output01, start_idx, stop_idx, stop_idx, id+"_r"/*, tmp_idx*/);
			}
			if (j == input_pos./*get(i).*/size()-1)
				recombine_last_pos(seq1, seq2, input_pos/*.get(i)*/, /*i,*/ output00, output01, id+"_r");

		}
		output00.close();
		output01.close();
	}

	/**
	 * @param seq1
	 * @param seq2
	 * @param input_pos_sublist
	 * @param output00
	 * @param output01
	 * @param id
	 */
	private void recombine_last_pos(String seq1, String seq2, ArrayList<Integer> input_pos_sublist, /*int idx,*/ 
			BufferedWriter output00, BufferedWriter output01, String id) {
		//Still something left to append! check the last index
		//inputData.getPositions().size()-1 = my last index
		//case 2: the last rec position < length of both allele
		//case 2: the last rec position = the last position of either of the allele=> do nothing
		
		/*int tmp_idx = idx*FrostRunner.chunk_length;*/
		int start = input_pos_sublist.get(input_pos_sublist.size()-1)/*-tmp_idx*/-1, stop1 = 0, stop2 = 0;

		if (input_pos_sublist.get(input_pos_sublist.size()-1)/*-tmp_idx*/ <= seq1.length()
				&& input_pos_sublist.get(input_pos_sublist.size()-1)/*-tmp_idx*/ <= seq2.length()){
			stop1 = seq1.length();
			stop2 = seq2.length();
			swap(seq1, seq2, output00, output01, start, stop1, stop2, id/*, tmp_idx*/);
//			System.out.println("RECOMBINE last: " + "\t" + start + "\t" + stop1 + "\t" + stop2);

		}
		else if (input_pos_sublist.get(input_pos_sublist.size()-1)/*-tmp_idx*/ == seq1.length()
				&& input_pos_sublist.get(input_pos_sublist.size()-1)/*-tmp_idx*/ < seq2.length()){
			stop1 = seq1.length()-1;
			stop2 = seq2.length();
			swap(seq1, seq2, output00, output01, start, stop1, stop2, id/*, tmp_idx*/);
//			System.out.println("RECOMBINE last: " + "\t" + start + "\t" + stop1 + "\t" + stop2);

		}
		else if (input_pos_sublist.get(input_pos_sublist.size()-1)/*-tmp_idx*/ < seq1.length()
				&& input_pos_sublist.get(input_pos_sublist.size()-1)/*-tmp_idx*/ == seq2.length()){
			stop1 = seq1.length();
			stop2 = seq2.length()-1;
			swap(seq1, seq2, output00, output01, start, stop1, stop2, id/*, tmp_idx*/);
//			System.out.println("RECOMBINE last: " + "\t" + start + "\t" + stop1 + "\t" + stop2);

		}

		else {
			stop1 = seq1.length();
			stop2 = seq2.length();
			swap(seq1, seq2, output00, output01, start, stop1, stop2, id/*, tmp_idx*/);
//			System.out.println("RECOMBINE last: " + "\t" + start + "\t" + stop1 + "\t" + stop2);

			
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
	protected void swap(String seq1, String seq2, BufferedWriter output00, BufferedWriter output01, 
			int start, int stop1, int stop2, String id/*, int tmp_idx*/) {
		Random rand = new Random();
		boolean flip = false;
		flip = rand.nextBoolean();
		if (flip) {
			unaltered_sequence(seq2, output00, null, start, stop2, false, id);
			unaltered_sequence(seq1, output01, null, start, stop1, false, id);
			this.recombined.add(id + "\t" + (start/*+tmp_idx*/+1) + "\t" + FrostRunner.parental_chromatids[1] + " " + FrostRunner.parental_chromatids[0]);
		}
		else {
			unaltered_sequence(seq1, output00, null, start, stop1, false, id);
			unaltered_sequence(seq2, output01, null, start, stop2, false, id);
			this.recombined.add(id + "\t" + (start/*+tmp_idx*/+1) + "\t" + FrostRunner.parental_chromatids[0] + " " + FrostRunner.parental_chromatids[1]);

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
	private void unaltered_sequence(String seq, BufferedWriter output00, BufferedWriter output10, 
			int start, int stop, boolean parent, String id) {
		
		if (stop > seq.length()){
			System.out.println("bad shit is happening");
			System.out.println(stop + " " + seq.length() + " " + id);
			return;
		}

		int char_idx = stop-start;
		if (char_idx <= 0) 
			return; /*nothing to write, consecutive mutation*/
		
//		System.out.println("Start: " + start + "\t" + "Stop: " + stop);
		char[] tmp00 = new char[char_idx];
		char[] tmp10 = new char[char_idx];
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
			BufferedWriter output00, BufferedWriter output01,
			BufferedWriter output10, BufferedWriter output11, 
			int stop, int stop_2, String id, boolean parent) {
		// TODO Auto-generated method stub
		ArrayList<Integer> adjust = new ArrayList<>(3);

		//Now we invoke mutation
		/*int tmp_idx = this.current_chunk*FrostRunner.chunk_length;*/
		String[] mutated = new String[2];
		String individuum = "";
		char[] toMutate = new char[2];
		
		if (parent) {
			toMutate[0] = seq1.charAt(stop);
			toMutate[1] = seq1.charAt(stop);
			
			Mutation m = new Mutation(toMutate, in.inDelCounter_M, in.inDelCounter_F);

			Random rand = new Random();
			boolean mutate0 = rand.nextBoolean();
			boolean mutate1 = rand.nextBoolean();
			if (mutate0 == true && mutate1 == false) {
				individuum = "M";
				mutated = m.mutationType(toMutate, individuum);
			} else if (mutate1 == true && mutate0 == false) {
				individuum = "F";
				mutated = m.mutationType(toMutate, individuum);
			} else if ((mutate0 == true && mutate1 == true)
					|| (mutate0 == false && mutate1 == false)) {
				individuum = "MF";
				mutated = m.mutationType(toMutate, individuum);
			}
			in.inDelCounter_M = m.M;
			in.inDelCounter_F = m.F;
			
//			System.out.println("indelcounter_F: " + in.inDelCounter_F);
//			System.out.println("indelcounter_M: " + in.inDelCounter_M);
			
			adjust.add(stop);
			adjust.add(in.inDelCounter_M);
			adjust.add(in.inDelCounter_F);
			this.sequenceAdjustment.add(adjust);
			
			this.ID1.add(id);
			//because the real position was at the index of the array which is always x-1
			this.POS1.add(stop/*+tmp_idx*/+1);
			this.prettyMutation1.add(m.getParentString());
			/**
			 * record for denovo adjustment
			 */

		}
		else {
//			int stop_F = adjustDeNovo(stop)[0], stop_M = adjustDeNovo(stop)[1];
//			System.out.println("Child: " + "\t"+ "stop_M: "+ stop+ "\t"+"stop_F: "+ stop_2);
//			System.out.println("Child: " + "\t"+ "length_M: "+ seq1.length()+ "\t"+"length_F: "+ seq2.length());


			if (stop < 0 || stop >= seq1.length() 
					|| stop_2 < 0 || stop_2 >= seq2.length()) /*if the position is negative it means 
			deletion happened in previous chunk, since we have no record of that anymore we skip this*/
				/*if the position is positive it means 
				insertion will happen in the next chunk, since we have no record of that we skip this*/
				return;
			toMutate[0] = seq1.charAt(stop);
			toMutate[1] = seq2.charAt(stop_2);
//			System.out.println("CHILD");
//			System.out.println("Previously: " +"\t" + stop + "\t"
//											+ "M: " +seq1.charAt(stop) + "\t"
//											+ "F: " +seq2.charAt(stop_2));
//			System.out.println("Currently M: " +"\t" + stop + "\t" +seq1.charAt(stop));
//			System.out.println("Currently F: " +"\t" + stop_2 + "\t" +seq2.charAt(stop_2));
			
			Mutation m = new Mutation(toMutate);

			individuum = "C";
			mutated = m.mutationType(toMutate, individuum);
//			System.out.println("toMutate: " +"\t" + toMutate[0] + "\t" +toMutate[1]);
//			System.out.println("mutated: " +"\t" + mutated[0] + "\t" +mutated[1]);
			this.ID2.add(id);
			this.POS2.add(stop/*+tmp_idx*/+1);
			this.prettyMutation2.add(m.getChildString());

		}
		write_ref_alt(individuum, mutated, toMutate, output00, output01, output10, output11);
	}
//		else prettyMutation.add("NO MUTATION");

	private int [] adjustDeNovo(int position) { // 0: pop, 1: mom
		// TODO Auto-generated method stub
		int p[] = new int[2]; //0: M, 1: F
		//+tmp_idx+1
//		System.out.println(position);

		for (int k = 0; k < this.sequenceAdjustment.size(); k++) {
			if (this.sequenceAdjustment.get(k).get(0) == position){
				p[0] = position + this.sequenceAdjustment.get(k).get(1);
				p[1] = position + this.sequenceAdjustment.get(k).get(2);

//				System.out.println(this.sequenceAdjustment.get(k).get(0) + "\t" + this.sequenceAdjustment.get(k).get(1));
//				System.out.println(this.sequenceAdjustment.get(k).get(0) + "\t" + this.sequenceAdjustment.get(k).get(2));

				break;
			}
			else if (this.sequenceAdjustment.get(k).get(0) > position) {
				if (k > 0 && k <= this.sequenceAdjustment.size()-1) {
					p[0] = position + this.sequenceAdjustment.get(k-1).get(1);
					p[1] = position + this.sequenceAdjustment.get(k-1).get(2);
					
//						System.out.println(this.sequenceAdjustment.get(k-1).get(0) + "\t" + this.sequenceAdjustment.get(k-1).get(1));
//						System.out.println(this.sequenceAdjustment.get(k-1).get(0) + "\t" + this.sequenceAdjustment.get(k-1).get(2));	
				}
				/*
				 * k = 0
				 */
				else {
					p[0] = position;
					p[1] = position;
				}

				break;
			}
			else if (k == this.sequenceAdjustment.size()-1 && this.sequenceAdjustment.get(k).get(0) < position) {
				p[0] = position + this.sequenceAdjustment.get(k).get(1);
				p[1] = position + this.sequenceAdjustment.get(k).get(2);
			}
			
			else
				continue;
						
		}
//		System.out.println(p);
//		System.out.println();
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
	private void write_ref_alt(String individuum, String[] allele, char[] toMutate, BufferedWriter output00,
			BufferedWriter output01, BufferedWriter output10, BufferedWriter output11) {
		//DEL
		String a1 = allele[0], a2 = allele[1];
		char c1 = toMutate[0], c2 = toMutate[1];

		if (allele[0].equals("*"))
			a1 = "";
		if (allele[1].equals("*"))
			a2 = "";

		if (individuum.equals("M")) {
			appendFile(output00, a1);
			appendFile(output01, a2);
			appendFile(output10, c1+"");
			appendFile(output11, c2+"");
		}
		if (individuum.equals("F")) {
			appendFile(output00, c1+"");
			appendFile(output01, c2+"");
			appendFile(output10, a1);
			appendFile(output11, a2);
		}
		if (individuum.equals("MF")) {
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
	private void appendFile(BufferedWriter bw, String s) {
		// TODO Auto-generated method stub
		try  { // new PrintWriter(new BufferedWriter(new FileWriter(fileName, true))))
			bw.write(s);
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

		if(this.prettyMutation1.size() != 0) {
			for (int i = 0; i < this.prettyMutation1.size(); i++) {
				data += this.ID1.get(i) + "\t" + this.POS1.get(i) + "\t" + this.prettyMutation1.get(i) + "\n";

			}
		}
		return data;
	}

	protected String getChildInfo() {

		/**
		 * CHILD ALWAYS GETS M/F AND NOT F/M
		 */
		String data = "";
		if(this.prettyMutation2.size() != 0) {
			for (int i = 0; i < this.prettyMutation2.size(); i++) {
				if (this.prettyMutation2.get(i) != "") {
					data += this.ID2.get(i) + "\t" + this.POS2.get(i) + "\t" + this.prettyMutation2.get(i) + "\t";
					data += FrostRunner.parental_chromatids[0]
							/*because allele 2 was programmed to be mother allele but conventionally mother allele comes first*/ 
							+ "\t" + FrostRunner.parental_chromatids[1]+ "\n";
				}
				else
					continue;			

			}
		}
		return data;
	}

}