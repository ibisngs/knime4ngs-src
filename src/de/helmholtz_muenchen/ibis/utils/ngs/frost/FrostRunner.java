package de.helmholtz_muenchen.ibis.utils.ngs.frost;

import java.text.DecimalFormat;
import java.text.NumberFormat;

/**
 * @author tanzeem.haque
 *
 */
public class FrostRunner {


	public static final String INTERNAL_PATH = "/home/ibis/tanzeem.haque/Documents/Frost_outputs/";
//	public static int SEED;
	public static String[] ID_List;



	public static void main(String[] args) throws InterruptedException {
		// TODO Auto-generated method stub

		//if -m not provided the mutation rate will be 2.36e-8;

//		for (int i = 0; i < args.length; i++) {
//			System.out.println((i+1) + " HERE " + args[i]);
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
		System.out.print("Execution time is (main) " + formatter.format((endTime - startTime) / 1000d) + " seconds");


	}

	/**
	 * @param input
	 * @param mutRate
	 * @param recombination
	 * @param seed
	 * @throws InterruptedException
	 */
	protected static void run(String input, double mutRate, int recombination, int seed) throws InterruptedException {
		// TODO Auto-generated method stub

		/**
		 * Reading reference Fasta
		 */	
		FastaReader fr = new FastaReader();
		fr.readSequenceFromFile(input);		
		/**
		 * preparing the mutations and recombination positions
		 */		
		InputScanner in = new InputScanner(fr, mutRate, recombination, seed);
		in.prepare();
		/**
		 * Creating the trio: invoke parental mutation, denovo for child,
		 * recombination for child
		 */		
		Trio_Simulator trio = new Trio_Simulator(in);
		trio.createTrio();	

		FrostRunner.ID_List = new String[trio.ID_List.size()];
		FrostRunner.ID_List = trio.ID_List.toArray(FrostRunner.ID_List);
		
		write(fr, in, trio, seed);		

	}

	private static void write(FastaReader fr, InputScanner in, Trio_Simulator trio, int seed) {
		OutputWriters ow = new OutputWriters(fr, in, trio);

		//Writing the IDs as file
		String ids = "";
		for (int i = 0; i < fr.size(); i++) {
			ids += ow.getFs().getIdentifier(i)+"\n";
		}
		ow.write_simple_string(INTERNAL_PATH + "ids.txt", ids);

/*
		//Writing info outputs
		for (int i = 0; i < trio.ID_List.size(); i++) {
			String elem = trio.ID_List.get(i);
			if (elem.equals(in.getiData_arrList_deNovo_child().get(i).getId())) {
				ow.write_InputData(INTERNAL_PATH + elem + "_deNovo_" + FrostRunner.SEED + ".txt",
						in.getiData_arrList_deNovo_child().get(i));
			}
			if (elem.equals(in.getiData_arrList_recombination_child().get(i).getId())) {
				ow.write_InputData(INTERNAL_PATH + elem + "_recombination_" + FrostRunner.SEED + ".txt",
						in.getiData_arrList_recombination_child().get(i));
			}
			
			
		}
*/		
		ow.write_InputData_arrList(INTERNAL_PATH + "deNovo_" + seed + ".txt",
				in.getiData_arrList_deNovo_child());
		ow.write_InputData_arrList(INTERNAL_PATH + "recombination_" + seed + ".txt",
				in.getiData_arrList_recombination_child());
		//writing the reference vcf for parents only
		ow.write_simple_string(FrostRunner.INTERNAL_PATH + "parents_run_" + seed + ".txt", trio.getParentInfo());
		ow.write_simple_string(FrostRunner.INTERNAL_PATH + "child_run_" + seed + ".txt", trio.getChildInfo());

		//writing the recombined file
		String rec = "";
		for (int i = 0; i < ow.getTrio().getRecombined().size(); i++) {
			rec += ow.getTrio().getRecombined().get(i) + "\n";
		}
		ow.write_simple_string(FrostRunner.INTERNAL_PATH + "recombined_seq_" + seed + ".txt", rec);
	
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
