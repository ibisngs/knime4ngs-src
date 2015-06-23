package de.helmholtz_muenchen.ibis.utils.ngs.frost;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
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
//	 --bed /storageNGS/scratch/Sequenciator/secondary_files/UCSC_CodingExons.bed

	public static ArrayList<String> id_list = new ArrayList<String> ();;
	public static String[] parental_chromatids = new String[2];
	public static String varyParameter = "";
	public static String[] records = FrostNodeModel.recordFiles();
	public static int skipped_N;
	public final static int reco_gap = 18000;
	public final static double insertion_rate = 0.15;
	public final static double deletion_rate = 0.05;
	public static BufferedWriter bw6;
;


	public static void main(String[] args) throws InterruptedException, IOException {
		// TODO Auto-generated method stub

		//if -m not provided the mutation rate will be 2.36e-8;

//		for (int i = 0; i < args.length; i++) {
//			System.out.println((i+1) + "\t" + args[i]);
//		}
		ArrayList<String> args_al = new ArrayList<>(Arrays.asList(args));

//		for (String s: args)
//			System.out.println("Runner: " + s);
		String tag_input = "-i", tag_mutRate = "-m", tag_recombination = "-r", 
				tag_generation = "-g", tag_seed = "-s",
				tag_mutVary = "--mut", tag_recVary = "--reco", tag_deNovoVary = "--denovo", tag_bed = "--bed";

		/**
		 * Default values are already given in Knime
		 */
		String input = "";
		String mapFile = FrostNodeModel.DEFAULT_MAPFILE;
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
		


		if (args.length >= 3 && args.length <= 13) {
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
						FrostRunner.varyParameter = "Mutation";
						continue;
					}
					else if (args[i].equals(tag_recVary)) {
						recVary = true;
						FrostRunner.varyParameter = "Crossover points";
						continue;
					}
					else if (args[i].equals(tag_deNovoVary)) {
						deNovoVary = true;
						FrostRunner.varyParameter = "de novo";
						continue;
					}
					else if (args[i].equals(tag_bed)) {
						mapFile = args[i + 1];
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
			 * Handling the existing files first
			 */
			deleteExistingFiles();

			String using_command = "";
			for (String s : args)
				using_command += s+" ";
			
			FrostRunner.bw6 = new BufferedWriter(new FileWriter(FrostRunner.records[6], true), 10000000);
			
			FrostRunner.createLog(FrostRunner.bw6, "Using command: " + "\n" + using_command);

			run(input, mapFile, mutRate, recombination, generation, seed, mutVary, recVary, deNovoVary);
			
		}
		long endTime   = System.currentTimeMillis();
		NumberFormat formatter = new DecimalFormat("#0.00000");
		FrostRunner.createLog(FrostRunner.bw6, "Execution time is (main) " + formatter.format((endTime - startTime) / 1000d) + " seconds");
		FrostRunner.createLog(FrostRunner.bw6,"\n");
		FrostRunner.bw6.close();
		System.out.println("done");


	}
	private static void deleteExistingFiles() {
		for (String s : FrostRunner.records) {
			File f = new File (s);
			if (f.exists())
				f.delete();
		}
		File f = new File (FrostRunner.records[0].replace(".tmp", ".txt"));
		if (f.exists())
			f.delete();
		f = new File (FrostRunner.records[1].replace(".tmp", ".txt"));
		if (f.exists())
			f.delete();
	}
	/**
	 * @param input
	 * @param mutRate
	 * @param recombination
	 * @param seed
	 * @throws InterruptedException
	 * @throws IOException 
	 */
	protected static void run(String input, String mapFile, double mutRate, int recombination, int generation, int seed, 
			boolean mutVary, boolean recVary, boolean deNovoVary) throws InterruptedException, IOException {
		// TODO Auto-generated method stub

		FrostRunner.createLog(FrostRunner.bw6,"Positions to vary: " + FrostRunner.varyParameter + ", using seed: " + seed);

		/**
		 * Checking input Fasta
		 */	
		FastaCheck fc = new FastaCheck(input);
		FrostRunner.parental_chromatids = fc.getParentalChromatids();
		FrostRunner.createLog(FrostRunner.bw6,"Parental chromatids: " + FrostRunner.parental_chromatids[0] 
				+ " and " + FrostRunner.parental_chromatids[1]);
		/**
		 * get the N region map or exons at first
		 */
		
		GenomeMap gm = new GenomeMap(mapFile);

		FrostRunner.createLog(FrostRunner.bw6,"File used for positions: " + mapFile);
		
		//Writing the IDs and Chunks as file
		RecordWriters rw = new RecordWriters();
		String ids = "";
		
		System.gc();
//		memory();

		/**
		 * Staring the vcf file
		 */
		BufferedWriter bw4 = new BufferedWriter(new FileWriter(FrostRunner.records[4], true), 10000000);
		rw.write_simple_string(bw4, 
				"#CHROM" + "\t" /*0*/
				+ "#POS" + "\t" /*1*/
				+ "#ID" + "\t" /*2*/
				+ "#REF" + "\t" /*3*/
				+ "#ALT" + "\t" /*4*/
				+ "#QUAL" + "\t" /*5*/
				+ "#FILTER" + "\t" /*6*/
				+ "#INFO" + "\t" /*7*/
				+ "#FORMAT" + "\t" /*8*/
				+ "#M" + "\t" /*9*/
				+ "#F" + "\t" /*10*/
				+ "#C" + "\n"); /*11*/
		bw4.close();

		/**
		 * 
		 * 
		 * 
		 * START OF A CHROMOSOME
		 * 
		 * 
		 * 
		 * 
		 */
		for (int i = 0; i < fc.input_chr_length.size(); i++) {
			BufferedWriter bw0 = new BufferedWriter(new FileWriter(FrostRunner.records[0], true), 10000000);
			BufferedWriter bw1 = new BufferedWriter(new FileWriter(FrostRunner.records[1], true), 10000000);
			BufferedWriter bw2 = new BufferedWriter(new FileWriter(FrostRunner.records[2], true), 10000000);
			BufferedWriter bw3 = new BufferedWriter(new FileWriter(FrostRunner.records[3], true), 10000000);
			bw4 = new BufferedWriter(new FileWriter(FrostRunner.records[4], true), 10000000);

			String currentChr = fc.input_chr_length.get(i).split("\t")[0];
			int currentLength = Integer.parseInt(fc.input_chr_length.get(i).split("\t")[1]);
			/**
			 * Some info printing
			 */
			FrostRunner.createLog(FrostRunner.bw6, (i+1) + ". "+ currentChr + " " + currentLength);
						
//			int chunk = (currentLength/FrostRunner.chunk_length)+1;
			
			ids += fc.input_chr_length.get(i) + "\n";

			/**
			 * Reading reference Fasta
			 */	
			FastaReader fr = new FastaReader();
			fr.readSequenceFromFile(input, currentChr);	



			/**
			 * preparing the mutations and recombination positions
			 * and also Creating the N Region maps
			 */		

			if (recombination > currentLength) {
				FrostRunner.createLog(FrostRunner.bw6, "Check the number of crossover points. Currently it is greater than the sequence length itself.");
				System.exit(0);
			}
			InputScanner in = new InputScanner(currentChr, mutRate, recombination, generation, seed, 
					mutVary, recVary, deNovoVary, gm);
			in.prepare(fr.getLength());//currentLength

			/**
			 * Creating the trio: invoke parental mutation, denovo for child,
			 * recombination for child
			 */	
			
			
//			for (int j = 0; j < chunk; j++) {

			/**
			 * ID_List will be called in the main method to write the output files.
			 */
			FrostRunner.id_list.add(/*j + "_" + */currentChr); //fc.input_chr_length.get(i).split("\t")[0] is the ID itself

			TrioSimulator trio = new TrioSimulator(/*fc, */fr, in);
			trio.createTrio(/*j, */currentChr);	

				
			/**
			 * Appending into the record files
			 */
			rw = new RecordWriters(fr, in, trio);
			/**
			 * 	
 			 * parents_run_
			 */
			rw.write_simple_string(bw0, trio.getParentInfo());
			/**
			 * child_run_
			 */
			rw.write_simple_string(bw1, trio.getChildInfo());
//			/**
//			 * deNovo_
//			 */
//			rw.write_InputData(recordFiles[2],in.getiData_deNovo_child()/*, j*/);
			/**
			 * recombination_
			 */
			rw.write_InputData(bw2,in.getiData_recombination_child()/*, j*/);
			//writing the recombined file
			String rec = "";
			for (int k = 0; k < rw.getTrio().getRecombined().size(); k++) {
				rec += rw.getTrio().getRecombined().get(k) + "\n";
			}
			/**
			 * recombined_seq_
			 */
			rw.write_simple_string(bw3, rec);

//			}
			
			bw0.flush();
			bw1.flush();
			bw2.flush();
			bw3.flush();
			FrostRunner.bw6.flush();
			System.gc();		
//			memory();
			bw0.close();
			bw1.close();
			bw2.close();
			bw3.close();
			/**
			 * trio_phased_
			 */
			rw.write_vcf(bw4, FrostRunner.records[0], FrostRunner.records[1]);
			bw4.close();

		}

		/**
		 * 
		 * 
		 * 
		 * END OF A CHROMOSOME!
		 * 
		 * 
		 * 
		 */
		/**
		 * trio_unpahsed_
		 */
		BufferedWriter bw5 = new BufferedWriter(new FileWriter(FrostRunner.records[5], true), 10000000);
		rw.unphase(bw5, FrostRunner.records[4]); //recordfile[5]
		bw5.close();


//		/**
//		 * ids_chunks.txt
//		 */
//		rw.write_simple_string(recordFiles[7], ids);
		/**
		 * delete the recordfile[0] (child_run_) and recordfile[1](parents_run_)
		 * which are .tmp s
		 */
		for (int i = 0; i <= 1; i++) {
			File f = new File (FrostRunner.records[i]);
			if (f.exists())
			f.delete();
		}
		
	}

	static void createLog(BufferedWriter bw, String logInfo) {
		
		try { // new PrintWriter(new BufferedWriter(new FileWriter(fileName, true))))
			bw.write(logInfo + "\n");
		} catch (Exception e) {
			System.err.println("Error: " + e.getMessage());
			e.printStackTrace();
		} 

	}
	static void memory() {
			
		int mb = 1024*1024;
        
        //Getting the runtime reference from system
        Runtime runtime = Runtime.getRuntime();
         
        FrostRunner.createLog(FrostRunner.bw6,"##### Heap utilization statistics [MB] #####");
         
        //Print used memory
        FrostRunner.createLog(FrostRunner.bw6,"Used Memory:" + (runtime.totalMemory() - runtime.freeMemory()) / mb);
 
        //Print free memory
        FrostRunner.createLog(FrostRunner.bw6,"Free Memory:" + runtime.freeMemory() / mb);
         
        //Print total available memory
        FrostRunner.createLog(FrostRunner.bw6,"Total Memory:" + runtime.totalMemory() / mb);
 
        //Print Maximum available memory
        FrostRunner.createLog(FrostRunner.bw6,"Max Memory:" + runtime.maxMemory() / mb);
		
	}
	
	public static void print_help() {

		System.out.println("Usage:");
		System.out.println("    -i <String> path to chromosome file in fasta format");
		System.out.println("    -m <Double>  Mutation rate (e.g. 1.5) " + "\n" +
				                "default: 2.36 (*e-8 per bp per generation)");
		System.out.println("    -r <Integer> number of switch positions for recombination"  + "\n" +
								"default: 1000");
		System.out.println("    -g <Integer> number of generations since the first Homo sapiens"  + "\n" +
								"default: 5300");
		System.out.println("    -s <Integer>  random seed for same output with two of the three paramters" + "\n" +
								"(mut/reco/denovo) being fixed");
		System.out.println("    --mut: mutation positions of parents will vary at each run");
		System.out.println("    --reco: crossover positions of child will vary at each run ");
		System.out.println("    --denovo: denovo positions of child will vary at each run ");

		System.out.println("Compulsory parameters are -i AND (--mut OR --reco OR --donovo)");
		System.out.println("Only --mut OR --reco OR --donovo at a time");


		// -verbose:gc -Dsun.rmi.dgc.client.gcInterval=3600000

//		System.exit(0);
	}

}
