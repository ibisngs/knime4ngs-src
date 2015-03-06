package de.helmholtz_muenchen.ibis.utils.ngs.frost;

import java.io.IOException;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;

/**
 * @author tanzeem.haque
 *
 */
public class FrostRunner {


	public static final String INTERNAL_PATH = "/home/ibis/tanzeem.haque/Documents/Frost_outputs/";
	public static ArrayList<String> ID_List = new ArrayList<String> ();
	public static int total_mutations;
	public static int total_deNovo;
	static String[] parental_chromatids = new String[2];

//	static int mutation_index_parent = 0;
//	static int denovo_index_child = 0;
//	static int recombination_index = 0;
	public final static int chunk_length = 10000000;


	public static void main(String[] args) throws InterruptedException, IOException {
		// TODO Auto-generated method stub

		//if -m not provided the mutation rate will be 2.36e-8;

//		for (int i = 0; i < args.length; i++) {
//			System.out.println((i+1) + "\t" + args[i]);
//		}
		
		String tag_input = "-i", tag_mutRate = "-m", tag_recombination = "-r", tag_seed = "-s";

		/**
		 * Default values are already given in Knime
		 */
		String input = "";
		double mutRate = 2.36;
		int recombination = 1000;
		int seed = 999;

		long startTime = System.currentTimeMillis();

		if (args.length < 2) {
			print_help();
		}

		else if (args.length >= 2 && args.length <= 8) {
			for (int i = 0; i < args.length; i++) {
				//Input-> must
				if (args[i].equals(tag_input)) {
					input = args[i + 1];
					continue;
				}
				//mutationrate -> if not there it will be 2.36e-8
				else if (args[i].equals(tag_mutRate)) {
					mutRate = Double.parseDouble(args[i + 1]);
					if (mutRate < 0) {
						print_help();
						System.out.println("MUT WRONG");
						break;
					}
					continue;
				}
				//recombination -> if not there it will be 1000
				else if (args[i].equals(tag_recombination)) {
					recombination = Integer.parseInt(args[i + 1]);
					if (recombination < 0) {
						print_help();
						System.out.println("REC WRONG");

						break;
					}
					continue;
				}
				else if (args[i].equals(tag_seed)) {
					seed = Integer.parseInt(args[i + 1]);
					continue;
				}
			}
		}
		else {
			System.out.println("I CAME IN ELSE");

			print_help();
		}

		if (input.equals("")) {
			System.out.println("INPUT WRONG");

			print_help();
		}

		else
			run(input, mutRate, recombination, seed);

		long endTime   = System.currentTimeMillis();
		NumberFormat formatter = new DecimalFormat("#0.00000");
		System.out.println("Execution time is (main) " + formatter.format((endTime - startTime) / 1000d) + " seconds");


	}

	/**
	 * @param input
	 * @param mutRate
	 * @param recombination
	 * @param seed
	 * @throws InterruptedException
	 * @throws IOException 
	 */
	protected static void run(String input, double mutRate, int recombination, int seed) throws InterruptedException, IOException {
		// TODO Auto-generated method stub

		/**
		 * Checking input Fasta
		 */	
		long startTime = System.currentTimeMillis();
		FastaCheck fc = new FastaCheck(input);
		FrostRunner.parental_chromatids = fc.getParentalChromatids();
		
		long endTime   = System.currentTimeMillis();
		NumberFormat formatter = new DecimalFormat("#0.00000");
		System.err.println("Execution time to check FASTA " + formatter.format((endTime - startTime) / 1000d) + " seconds" + "\n");

//		System.out.println(fc.input_chr_length.size());
		for (int i = 0; i < fc.input_chr_length.size(); i++)
			System.out.println((i+1) +". " + fc.input_chr_length.get(i));
		
		//Writing the IDs and Chunks as file
		RecordWriters rw = new RecordWriters();
		String ids = "";
		
		System.gc();
		memory();

		for (int i = 0; i < fc.input_chr_length.size(); i++) {
			String currentChr = fc.input_chr_length.get(i).split("\t")[0];
			int currentLength = Integer.parseInt(fc.input_chr_length.get(i).split("\t")[1]);
			/**
			 * Some info printing
			 */
			System.out.println();
			System.out.println((i+1) + ". "+ currentChr + " " + currentLength);
			FrostRunner.total_mutations = (int)((currentLength/Math.pow(10, 8) * mutRate * 5300)/(1));
			FrostRunner.total_deNovo = (int)(FrostRunner.total_mutations/5300);
			System.out.println("Number of total mutations: " + FrostRunner.total_mutations);
			System.out.println("Number of total de novo: " + FrostRunner.total_deNovo);
			
			int chunk = (currentLength/FrostRunner.chunk_length)+1;
			
			ids += fc.input_chr_length.get(i) + "\n";

			/**
			 * Reading reference Fasta
			 */	
			startTime = System.currentTimeMillis();
			FastaReader fr = new FastaReader();
			fr.readSequenceFromFile(input, currentChr);	
			System.out.println("Actual length of " + currentChr + " from file: " + fr.getLength());

			endTime   = System.currentTimeMillis();
			formatter = new DecimalFormat("#0.00000");
			System.err.println("Execution time to read " + currentChr + " FASTA: " + formatter.format((endTime - startTime) / 1000d) + " seconds");
		

			/**
			 * preparing the mutations and recombination positions
			 */	

			startTime = System.currentTimeMillis();
			InputScanner in = new InputScanner(mutRate, recombination, seed, chunk);
			in.prepare(currentChr, fr.getLength());//currentLength
			endTime   = System.currentTimeMillis();
			formatter = new DecimalFormat("#0.00000");
			System.err.println("Execution time to prepare " + currentChr + ": " + formatter.format((endTime - startTime) / 1000d) + " seconds");

				/**
				 * Creating the trio: invoke parental mutation, denovo for child,
				 * recombination for child
				 */	
			
			System.out.println("CHUNK SIZE: " + chunk + "\t" + currentChr);
			startTime = System.currentTimeMillis();
		
			
			for (int j = 0; j < chunk; j++) {

//				System.err.println("P1: " + FrostRunner.parental_chromatids[0] + "\t" + "P2: " + FrostRunner.parental_chromatids[1]);

				/**
				 * ID_List will be called in the main method to write the output files.
				 */
				FrostRunner.ID_List.add(j + "_" + currentChr); //fc.input_chr_length.get(i).split("\t")[0] is the ID itself

				startTime = System.currentTimeMillis();
				System.out.println("CHUNKY #"+j);
				
//				System.out.println("mutation_index_parents: " + mutation_index_parent);
//				System.out.println("denovo_index_child: " + denovo_index_child);
//				System.out.println("recombination_index: " + recombination_index);

				Trio_Simulator trio = new Trio_Simulator(fc, fr, in);
				trio.createTrio(j, currentChr);	
				
				/**
				 * Appending into the record files
				 */
				rw = new RecordWriters(fr, in, trio);
				rw.write_simple_string(FrostRunner.INTERNAL_PATH + "parents_run_" + seed + ".txt", trio.getParentInfo());
				rw.write_simple_string(FrostRunner.INTERNAL_PATH + "child_run_" + seed + ".txt", trio.getChildInfo());

				rw.write_InputData(INTERNAL_PATH + "deNovo_" + seed + ".txt",in.getiData_deNovo_child(), j);			
				rw.write_InputData(INTERNAL_PATH + "recombination_" + seed + ".txt",in.getiData_recombination_child(), j);
				//writing the recombined file
				String rec = "";
				for (int k = 0; k < rw.getTrio().getRecombined().size(); k++) {
					rec += rw.getTrio().getRecombined().get(k) + "\n";
				}
				rw.write_simple_string(FrostRunner.INTERNAL_PATH + "recombined_seq_" + seed + ".txt", rec);


			}
		
			endTime   = System.currentTimeMillis();
			formatter = new DecimalFormat("#0.00000");
			System.err.println("Execution time to create trio of " + currentChr + ": " + formatter.format((endTime - startTime) / 1000d) + " seconds");
			
//			rw.write_InputData(INTERNAL_PATH + "deNovo_" + seed + ".txt",in.getiData_deNovo_child());			
//			mutation_index_parent = 0;
//			denovo_index_child = 0;
//			recombination_index = 0;
			System.gc();		
			memory();
			System.out.println();

		}

		rw.write_simple_string(INTERNAL_PATH + "ids_chunk.txt", ids);
		


	}

	static void memory() {
			
		int mb = 1024*1024;
        
        //Getting the runtime reference from system
        Runtime runtime = Runtime.getRuntime();
         
        System.out.println("##### Heap utilization statistics [MB] #####");
         
        //Print used memory
        System.out.println("Used Memory:" + (runtime.totalMemory() - runtime.freeMemory()) / mb);
 
        //Print free memory
        System.out.println("Free Memory:" + runtime.freeMemory() / mb);
         
        //Print total available memory
        System.out.println("Total Memory:" + runtime.totalMemory() / mb);
 
        //Print Maximum available memory
        System.out.println("Max Memory:" + runtime.maxMemory() / mb);
		
	}
	
	public static void print_help() {
		System.out.println("Usage:");
		System.out.println("    -i <String> path to chromosome file in fasta format");
		System.out.println("    -m <Double>  Mutation rate (e.g. 1.5) " + "\n" +
				                "default: 2.36 (*e-8 per bp per generation)");
		System.out.println("    -r <Integer> number of switch positions for recombination"  + "\n" +
								"default: 1000");
		System.out.println("    -s <Integer>  random seed for same output");
		System.out.println("Compulsory parameter is -i");

		// -verbose:gc -Dsun.rmi.dgc.client.gcInterval=3600000

//		System.exit(0);
	}

}
