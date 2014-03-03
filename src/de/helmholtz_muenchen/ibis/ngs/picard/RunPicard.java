package de.helmholtz_muenchen.ibis.ngs.picard;

import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.nio.file.Paths;

import net.sf.picard.analysis.CollectInsertSizeMetrics;
import net.sf.picard.sam.AddOrReplaceReadGroups;
import net.sf.picard.sam.MarkDuplicates;
import net.sf.picard.sam.SortSam;
import net.sf.samtools.SAMFileWriterFactory;

public class RunPicard {
	
	//run CollectInsertSizeMetrics
	protected static void runMetrics (String input, String ref, String outputm, String outputh, String val, boolean index, double dev, double pct, String acc, boolean sorted ) throws Exception{
		
		/* command line options
		 * INPUT
		 * CREATE_INDEX
		 * VALIDATION_STRINGENCY
		 * HISTOGRAM_FILE
		 * DEVIATIONS
		 * MINIMUM_PCT
		 * METRIC_ACCUMULATION_LEVEL
		 * ASSUME_SORTED
		 * OUTPUT
		 * REFERENCE_SEQUENCE
		 */
		
		String[] args = new String[11];
		args[0]="INPUT="+input;
		args[1]="OUTPUT="+outputm;
		args[2]="HISTOGRAM_FILE="+outputh;
		args[3]="REFERENCE_SEQUENCE="+ref;
		args[4]="VALIDATION_STRINGENCY="+val;
		args[5]="CREATE_INDEX="+index;
		args[6]="DEVIATIONS="+dev;
		args[7]="MINIMUM_PCT="+pct;
		args[8]="METRIC_ACCUMULATION_LEVEL="+acc;
		args[9]="ASSUME_SORTED="+sorted;
		args[10]="TMP_DIR="+Paths.get(input).getParent().toString();
		
		PicardToolsNodeModel.logger.info("Running CollectInsertSizeMetrics "+args[0]+" "+args[1]+" "+args[2]+" "+args[3]+" "+args[4]+" "+args[5]+" "+args[6]+" "+args[7]+" "+args[8]+" "+args[9]);
		
		//redirect error stream
		redirecterr(outputm);
		System.err.println("Output of CollectInsertSizeMetrics:");
		
		Exception exception=null;
		
		int exitcode =0;
		
		try{
			//run tool
			exitcode=new CollectInsertSizeMetrics().instanceMain(args);
		}
		catch(Exception e){
			exception =e ;
		}
		//reset error stream and variable for index writing
		finally{
			SAMFileWriterFactory.setDefaultCreateIndexWhileWriting(false);
			reseterr();
		}
		
		//checks if exception has been caught
		if(exception!=null){
			throw exception;
		}
		
		if(exitcode==1){
			throw new Exception("Something went wrong while executing CollectInsertSizeMetrics, please check out the log file");
		}
		
		PicardToolsNodeModel.logger.info("CollectInsertSizeMetrics finished successfully");
	}
	
	//run AddOrReplaceReadGroups
	protected static void runRG(String input, String output, String val, boolean index, String id, String library, String sample, String unit, String platform) throws Exception{
		
		/* command line arguments
		 * INPUT
		 * OUTPUT
		 * VALIDATION_STRINGENCY
		 * CREATE_INDEX
		 * RGID id
		 * RGLB library
		 * RGPL platform
		 * RGPU platform unit
		 * RGSM sample
		 */
		
		String [] args = new String[10];
		args[0]="INPUT="+input;
		args[1]="OUTPUT="+output;
		args[2]="VALIDATION_STRINGENCY="+val;
		args[3]="CREATE_INDEX="+index;
		args[4]="RGID="+id;
		args[5]="RGLB="+library;
		args[6]="RGSM="+sample;
		args[7]="RGPU="+unit;
		args[8]="RGPL="+platform;
		args[9]="TMP_DIR="+Paths.get(input).getParent().toString();
		
		PicardToolsNodeModel.logger.info("Running AddOrReplaceReadGroups "+args[0]+" "+args[1]+" "+args[2]+" "+args[3]+" "+args[4]+" "+args[5]+" "+args[6]+" "+args[7]+" "+args[8]);
		
		//redirect error stream
		redirecterr(output);
		System.err.println("Output of AddOrReplaceReadGroups:");
		
		Exception exception=null;
		int exitcode=0;
		
		try{
			//run tool
			exitcode=new AddOrReplaceReadGroups().instanceMain(args);
		}
		catch(Exception e){
			exception =e;
		}
		//reset error stream, resert writing index variable
		finally{
			SAMFileWriterFactory.setDefaultCreateIndexWhileWriting(false);
			reseterr();
		}
		
		//checks if exception has been caught
		if(exception!=null){
			throw exception;
		}
		
		if(exitcode==1){
			throw new Exception("Something went wrong while executing MarkDuplicates, please check out the log file");
		}
		
		PicardToolsNodeModel.logger.info("AddOrReplaceReadGroups finished successfully");
	}

	
	//run MarkDuplicates
	protected static void runDupl(String input, String output, String metrics, String val, boolean index, boolean rmdupl, boolean ass_sort) throws Exception {
		
		/* command line options for MarkDuplicates
		 * INPUT
		 * OUTPUT
		 * METRICS_FILE
		 * REMOVE_DUPLICATES
		 * ASSUME_SORTED
		 * VALIDATION_STRINGENCY
		 * CREATE_INDEX
		 */
		
		String [] args = new String[8];
		args[0]="INPUT="+input;
		args[1]="OUTPUT="+output;
		args[2]="METRICS_FILE="+metrics;
		args[3]="VALIDATION_STRINGENCY="+val;
		args[4]="CREATE_INDEX="+index;
		args[5]="REMOVE_DUPLICATES="+rmdupl;
		args[6]="ASSUME_SORTED="+ass_sort;
		args[7]="TMP_DIR="+Paths.get(input).getParent().toString();
		
		PicardToolsNodeModel.logger.info("Running MarkDuplicates "+args[0]+" "+args[1]+" "+args[2]+" "+args[3]+" "+args[4]+" "+args[5]+" "+args[6]);
		
		//redirect error stream
		redirecterr(output);
		System.err.println("Output of MarkDuplicates:");
		
		Exception exception=null;
		int exitcode=0;

		try{
			//run tool
			exitcode = new MarkDuplicates().instanceMain(args);
		}
		catch(Exception e){
			exception = e;
		}
		//reset error stream, reset writing index variable
		finally{
			SAMFileWriterFactory.setDefaultCreateIndexWhileWriting(false);
			reseterr();
		}
		
		//checks if exception has been caught
		if(exception!=null){
			throw exception;
		}
		
		if(exitcode==1){
			throw new Exception("Something went wrong while executing MarkDuplicates, please check out the log file");
		}
		
		PicardToolsNodeModel.logger.info("MarkDupplicates finished successfully");
		
	}
	
	//run SortSam
	protected static void runSortSam(String input, String output, String val, boolean index, String order) throws Exception {
		
		/* command line options for SortSam
		 *INPUT
		 *OUTPUT
		 *SORT_ORDER
		 *VALIDATION_STRINGENCY
		 *CREATE_INDEX
		 */

		
		String[] args= new String[6];
		args[0]="INPUT="+input;
		args[1]="OUTPUT="+output;
		args[2]="VALIDATION_STRINGENCY="+val;
		args[3]="CREATE_INDEX="+index;
		args[4]="SORT_ORDER="+order;
		args[5]="TMP_DIR="+Paths.get(input).getParent().toString();
		
		PicardToolsNodeModel.logger.info("Running SortSam "+args[0]+" "+args[1]+" "+args[2]+" "+args[3]+" "+args[4]);
		
		
		//redirect error stream
		redirecterr(output);
		System.err.println("Output of Sortsam");
		
		Exception exception=null;
		int exitcode=0;
		
		try{
			//run tool
			exitcode=new SortSam().instanceMain(args);
		}
		catch(Exception e){
			exception = e;
		}
		//reset error stream, reset index writing variable
		finally{
			SAMFileWriterFactory.setDefaultCreateIndexWhileWriting(false);
			reseterr();
		}
		
		//checks if exception has been caught
		if(exception!=null){
			throw exception;
		}
		if(exitcode==1){
			throw new Exception("Something went wrong while executing SortSam, please check out the log file");
		}
		
		PicardToolsNodeModel.logger.info("SortSam finished successfully");

	}
	
	
	//create file that contains output of picard tools
	//code from http://www.avajava.com/tutorials/lessons/how-do-i-redirect-standard-error-to-a-file.html
	
	private static PrintStream err=null;
	
	//method to redirect error stream
	private static void redirecterr(String output) throws Exception{
		
		//makes sure err is never overridden
		if(err==null){
			err = System.err;
		}
		PrintStream errtofile = new PrintStream( new FileOutputStream(new File(output+".log")));
		System.setErr(errtofile);

		PicardToolsNodeModel.logger.info("log file can be found in "+output+".log");
	}

	//method to reset error stream
	private static void reseterr(){
		System.setErr(err);
	}
	
		
	}
	