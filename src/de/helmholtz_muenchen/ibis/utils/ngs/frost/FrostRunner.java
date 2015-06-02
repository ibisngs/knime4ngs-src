package de.helmholtz_muenchen.ibis.utils.ngs.frost;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;

import de.helmholtz_muenchen.ibis.ngs.frost.FrostNodeModel;

/**
 * @author tanzeem.haque
 *
 */
public class FrostRunner {


	public static final String INTERNAL_PATH = FrostNodeModel.INTERNAL_OUTPUT_PATH;
	public static ArrayList<String> id_list = new ArrayList<String> ();
	public static int total_mutations;
	public static int total_deNovo;
	public static String[] parental_chromatids = new String[2];
	public static String[] records = new String[6];

//	static int mutation_index_parent = 0;
//	static int denovo_index_child = 0;
//	static int recombination_index = 0;
	public static int skipped_N = 0;
	public static int mutationCounter = 0;

	public final static int chunk_length = 10000000;
	public final static double insertion_rate = 0.15;
	public final static double deletion_rate = 0.05;

	public static void main(String[] args) throws InterruptedException, IOException {
		// TODO Auto-generated method stub

		//if -m not provided the mutation rate will be 2.36e-8;

//		for (int i = 0; i < args.length; i++) {
//			System.out.println((i+1) + "\t" + args[i]);
//		}
		
		ArrayList<String> args_al = new ArrayList<>(Arrays.asList(args));

		String tag_input = "-i", tag_mutRate = "-m", tag_recombination = "-r", 
				tag_generation = "-g", tag_seed = "-s",
				tag_mutVary = "--mut", tag_recVary = "--reco", tag_deNovoVary = "--denovo";

		/**
		 * Default values are already given in Knime
		 */
		String input = "";
		double mutRate = 2.36;
		int recombination = 1000;
		int generation = 5300;
		int seed = 999;
		boolean mutVary = false;
		boolean recVary = false;
		boolean deNovoVary = false;
		
		long startTime = System.currentTimeMillis();

		if (args.length < 3) {
			print_help();
			return;
		}
		


		if (args.length >= 3 && args.length <= 11) {
			if (!(args_al.contains(tag_mutVary) /* no booleans*/
					|| args_al.contains(tag_recVary) 
					|| args_al.contains(tag_deNovoVary))) {
					print_help();
					System.out.println("I CAME IN NO BOOLEANS");
					return;
			}
			if (args_al.contains(tag_mutVary)  /* all 3 booleans*/
					&& args_al.contains(tag_recVary) 
					&& args_al.contains(tag_deNovoVary)) {
					print_help();
					System.out.println("I CAME IN 3 BOOLEANS");
					return;
			}
			if (args_al.contains(tag_mutVary) /* two booleans */
					&& args_al.contains(tag_recVary)
					&& !args_al.contains(tag_deNovoVary)) {
					print_help();
					System.out.println("I CAME IN 2 (mut rec) BOOLEANS");
					return;
			}
			if (args_al.contains(tag_mutVary)  /* two booleans */
					&& args_al.contains(tag_deNovoVary)
					&& !args_al.contains(tag_recVary)) {
					print_help();
					System.out.println("I CAME IN 2 (mut denovo) BOOLEANS");
					return;
			}
			if (args_al.contains(tag_recVary) /* two booleans */
					&& args_al.contains(tag_deNovoVary)
					&& !args_al.contains(tag_mutVary)) {
					print_help();
					System.out.println("I CAME IN 2 (rec denovo) BOOLEANS");
					return;
			}
			else {
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
					//generation -> if not there it will be 5300
					else if (args[i].equals(tag_generation)) {
						generation = Integer.parseInt(args[i + 1]);
						if (generation < 0) {
							print_help();
							System.out.println("GENERATION NEXT");

							break;
						}
						continue;
					}
					else if (args[i].equals(tag_seed)) {
						seed = Integer.parseInt(args[i + 1]);
						continue;
					}
					else if (args[i].equals(tag_mutVary)) {
						mutVary = true;
						continue;
					}
					else if (args[i].equals(tag_recVary)) {
						recVary = true;
						continue;
					}
					else if (args[i].equals(tag_deNovoVary)) {
						deNovoVary = true;
						continue;
					}
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

		else {
			/**
			 *
			 * files[0] = FrostNodeModel.INTERNAL_OUTPUT_PATH + "parents_run_" + FrostNodeModel.seed + ".txt";
    	files[1] = FrostNodeModel.INTERNAL_OUTPUT_PATH + "child_run_" + FrostNodeModel.seed + ".txt";
    	files[2] = FrostNodeModel.INTERNAL_OUTPUT_PATH + "deNovo_" + FrostNodeModel.seed + ".txt";
    	files[3] = FrostNodeModel.INTERNAL_OUTPUT_PATH + "recombination_" + FrostNodeModel.seed + ".txt";
    	files[4] = FrostNodeModel.INTERNAL_OUTPUT_PATH + "recombined_seq_" + FrostNodeModel.seed + ".txt";
		files[5] = FrostNodeModel.INTERNAL_OUTPUT_PATH + "trio_" + FrostNodeModel.seed + ".vcf";
    	files[6] = FrostNodeModel.INTERNAL_OUTPUT_PATH + "ids_chunks.txt";    	
    	
			 */
			/**
			 * Handling the existing files first
			 */
//			FrostRunner.records = FrostNodeModel.recordFiles();
			for (String s : records) {
				File f = new File (s);
				if (f.exists())
					f.delete();
			}

//			System.out.println("Mut: " + mutVary + "\t" + "Rec: " + recVary + "\t" + "Denovo: " + deNovoVary);
			run(input, mutRate, recombination, generation, seed, mutVary, recVary, deNovoVary, records);
		}
		long endTime   = System.currentTimeMillis();
		NumberFormat formatter = new DecimalFormat("#0.00000");
		FrostRunner.createLog("Skipped N's: " + FrostRunner.skipped_N);
		FrostRunner.createLog("Execution time is (main) " + formatter.format((endTime - startTime) / 1000d) + " seconds");

	}

	/**
	 * @param input
	 * @param mutRate
	 * @param recombination
	 * @param seed
	 * @throws InterruptedException
	 * @throws IOException 
	 */
	protected static void run(String input, double mutRate, int recombination, int generation, int seed, 
			boolean mutVary, boolean recVary, boolean deNovoVary, String[] recordFiles) throws InterruptedException, IOException {
				// TODO Auto-generated method stub

		/**
		 * Checking input Fasta
		 */	
//		long startTime = System.currentTimeMillis();
		FastaCheck fc = new FastaCheck(input);
		FrostRunner.parental_chromatids = fc.getParentalChromatids();
		
//		long endTime   = System.currentTimeMillis();
//		NumberFormat formatter = new DecimalFormat("#0.00000");
//		System.err.println("Execution time to check FASTA " + formatter.format((endTime - startTime) / 1000d) + " seconds" + "\n");

//		System.out.println(fc.input_chr_length.size());
//		for (int i = 0; i < fc.input_chr_length.size(); i++)
//			System.out.println((i+1) +". " + fc.input_chr_length.get(i));
		
		//Writing the IDs and Chunks as file
		RecordWriters rw = new RecordWriters();
		String ids = "";
		
		System.gc();
//		memory();

		for (int i = 0; i < fc.input_chr_length.size(); i++) {
			String currentChr = fc.input_chr_length.get(i).split("\t")[0];
			int currentLength = Integer.parseInt(fc.input_chr_length.get(i).split("\t")[1]);
			/**
			 * Some info printing
			 */
			FrostRunner.createLog("\n");
			FrostRunner.createLog((i+1) + ". "+ currentChr + " " + currentLength);
			FrostRunner.total_mutations = (int)((currentLength/Math.pow(10, 8) * mutRate * 5300)/(1));
			FrostRunner.total_deNovo = (int)(FrostRunner.total_mutations/5300);
			FrostRunner.createLog("Number of total mutations: " + FrostRunner.total_mutations);
			FrostRunner.createLog("Number of total de novo: " + FrostRunner.total_deNovo);
			
			int chunk = (currentLength/FrostRunner.chunk_length)+1;
			
			ids += fc.input_chr_length.get(i) + "\n";

			/**
			 * Reading reference Fasta
			 */	
//			startTime = System.currentTimeMillis();
			FastaReader fr = new FastaReader();
			fr.readSequenceFromFile(input, currentChr);	
//			System.out.println("Actual length of " + currentChr + " from file: " + fr.getLength());

//			endTime   = System.currentTimeMillis();
//			formatter = new DecimalFormat("#0.00000");
//			System.err.println("Execution time to read " + currentChr + " FASTA: " + formatter.format((endTime - startTime) / 1000d) + " seconds");
		

			/**
			 * preparing the mutations and recombination positions
			 */	

//			startTime = System.currentTimeMillis();
			InputScanner in = new InputScanner(currentChr, mutRate, recombination, generation, seed, 
					mutVary, recVary, deNovoVary, chunk);			
			in.prepare(fr.getLength());//currentLength
//			endTime   = System.currentTimeMillis();
//			formatter = new DecimalFormat("#0.00000");
//			System.err.println("Execution time to prepare " + currentChr + ": " + formatter.format((endTime - startTime) / 1000d) + " seconds");

				/**
				 * Creating the trio: invoke parental mutation, denovo for child,
				 * recombination for child
				 */	
			
			FrostRunner.createLog("CHUNK SIZE: " + chunk + "\t" + currentChr);
//			startTime = System.currentTimeMillis();
		
			
			for (int j = 0; j < chunk; j++) {

//				System.err.println("P1: " + FrostRunner.parental_chromatids[0] + "\t" + "P2: " + FrostRunner.parental_chromatids[1]);

				/**
				 * ID_List will be called in the main method to write the output files.
				 */
				FrostRunner.id_list.add(j + "_" + currentChr); //fc.input_chr_length.get(i).split("\t")[0] is the ID itself

//				startTime = System.currentTimeMillis();
				FrostRunner.createLog("CHUNKY #"+j);
				
//				System.out.println("mutation_index_parents: " + mutation_index_parent);
//				System.out.println("denovo_index_child: " + denovo_index_child);
//				System.out.println("recombination_index: " + recombination_index);

				TrioSimulator trio = new TrioSimulator(fc, fr, in);
				trio.createTrio(j, currentChr);	
				
				/**
				 * Appending into the record files
				 */
				rw = new RecordWriters(fr, in, trio);
				/**
				 * parents_run_
				 */
				rw.write_simple_string(recordFiles[0], trio.getParentInfo());
				/**
				 * child_run_
				 */
				rw.write_simple_string(recordFiles[1], trio.getChildInfo());
				/**
				 * deNovo_
				 */
				rw.write_InputData(recordFiles[2],in.getiData_deNovo_child(), j);
				/**
				 * recombination_
				 */
				rw.write_InputData(recordFiles[3],in.getiData_recombination_child(), j);
				//writing the recombined file
				String rec = "";
				for (int k = 0; k < rw.getTrio().getRecombined().size(); k++) {
					rec += rw.getTrio().getRecombined().get(k) + "\n";
				}
				/**
				 * recombined_seq_
				 */
				rw.write_simple_string(recordFiles[4], rec);


			}
			
			/**
			 * trio_
			 */
			rw.write_vcf(recordFiles[5], recordFiles[0], recordFiles[1]);
			rw.unphase(recordFiles[5]);


//			endTime   = System.currentTimeMillis();
//			formatter = new DecimalFormat("#0.00000");
//			System.err.println("Execution time to create trio of " + currentChr + ": " + formatter.format((endTime - startTime) / 1000d) + " seconds");
			
//			rw.write_InputData(INTERNAL_PATH + "deNovo_" + seed + ".txt",in.getiData_deNovo_child());			
//			mutation_index_parent = 0;
//			denovo_index_child = 0;
//			recombination_index = 0;
			System.gc();		
//			memory();
//			System.out.println();

		}
		/**
		 * ids_chunks.txt
		 */
		rw.write_simple_string(recordFiles[6], ids);
		


	}

	static void createLog(String logInfo) {
		
		try (PrintWriter pw = new PrintWriter(new FileOutputStream(FrostRunner.records[7], true))) { // new PrintWriter(new BufferedWriter(new FileWriter(fileName, true))))
			pw.write(logInfo + "\n");
		} catch (Exception e) {
			System.err.println("Error: " + e.getMessage());
			e.printStackTrace();
		}
//		Logger logger = Logger.getLogger("SimulatorLog");  
//	    FileHandler fh;  
//
//	    try {  
//
//	        // This block configure the logger with handler and formatter  
//	        fh = new FileHandler(INTERNAL_PATH + "MyLogFile.log");  
//	        logger.addHandler(fh);
//	        logger.setUseParentHandlers(false);
//	        SimpleFormatter formatter = new SimpleFormatter();  
//	        fh.setFormatter(formatter);  
//
//	        // the following statement is used to log any messages  
//	        logger.info(logInfo);  
//
//	    } catch (SecurityException e) {  
//	        e.printStackTrace();  
//	    } catch (IOException e) {  
//	        e.printStackTrace();  
//	    }  

	}
	static void memory() {
			
		int mb = 1024*1024;
        
        //Getting the runtime reference from system
        Runtime runtime = Runtime.getRuntime();
         
        FrostRunner.createLog("##### Heap utilization statistics [MB] #####");
         
        //Print used memory
        FrostRunner.createLog("Used Memory:" + (runtime.totalMemory() - runtime.freeMemory()) / mb);
 
        //Print free memory
        FrostRunner.createLog("Free Memory:" + runtime.freeMemory() / mb);
         
        //Print total available memory
        FrostRunner.createLog("Total Memory:" + runtime.totalMemory() / mb);
 
        //Print Maximum available memory
        FrostRunner.createLog("Max Memory:" + runtime.maxMemory() / mb);
		
	}
	
	public static void print_help() {
		FrostRunner.createLog("Usage:");
		FrostRunner.createLog("    -i <String> path to chromosome file in fasta format");
		FrostRunner.createLog("    -m <Double>  Mutation rate (e.g. 1.5) " + "\n" +
				                "default: 2.36 (*e-8 per bp per generation)");
		FrostRunner.createLog("    -r <Integer> number of switch positions for recombination"  + "\n" +
								"default: 1000");
		FrostRunner.createLog("    -g <Integer> number of generations since the first Homo sapiens"  + "\n" +
								"default: 5300");
		FrostRunner.createLog("    -s <Integer>  random seed for same output");
		FrostRunner.createLog("Compulsory parameter is -i");

		// -verbose:gc -Dsun.rmi.dgc.client.gcInterval=3600000

//		System.exit(0);
	}

}
